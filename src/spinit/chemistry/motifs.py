"""Connectivity-aware local motif detection for multi-element structures.

Motif labels in this module are heuristic connectivity descriptors. Names ending
in ``_like`` indicate local-environment similarity and are not formal assignments
of bond order, oxidation state, aromaticity, or molecular charge.
"""

from __future__ import annotations

from collections import Counter
from typing import Any, Mapping, Sequence

import networkx as nx
from ase import Atoms

from ..graph.graph_ops import (
    count_neighbors_by_element,
    get_bond_length_between,
    get_neighbor_bond_length_map,
    get_neighbors,
)
from .hybridization import is_sp2_like_label


def build_component_context(atoms: Atoms, G: nx.Graph, config: Mapping[str, Any]) -> dict[str, Any]:
    """Build connectivity metadata to distinguish network-bound vs fragment atoms.

    Context: Motif meaning changes between extended carbon networks and small
    disconnected molecules, so this component map is used to separate those regimes.
    """
    carbon_symbol = str(config["carbon_symbol"])
    mcfg = config.get("motif_detection", {})
    network_min_size = int(mcfg.get("network_min_component_size", 20))
    network_requires_carbon = bool(mcfg.get("network_requires_carbon", True))
    molecular_max_size = int(mcfg.get("molecular_fragment_max_size", 16))

    components = [sorted(int(v) for v in comp) for comp in nx.connected_components(G)]
    if len(components) == 0:
        components = [[int(i)] for i in range(len(atoms))]

    atom_to_component: dict[int, int] = {}
    component_meta: dict[int, dict[str, Any]] = {}

    network_component_ids: set[int] = set()
    for cid, nodes in enumerate(components):
        counts = Counter(str(atoms[i].symbol) for i in nodes)
        size = int(len(nodes))
        contains_carbon = bool(counts.get(carbon_symbol, 0) > 0)
        is_network_like = bool(
            size >= network_min_size and (contains_carbon or (not network_requires_carbon))
        )

        if is_network_like:
            network_component_ids.add(cid)

        component_meta[cid] = {
            "size": size,
            "element_counts": {k: int(v) for k, v in sorted(counts.items())},
            "contains_carbon": contains_carbon,
            "is_network_like": is_network_like,
        }
        for i in nodes:
            atom_to_component[int(i)] = cid

    if len(network_component_ids) == 0:
        fallback = [
            cid
            for cid, meta in component_meta.items()
            if meta["contains_carbon"] and int(meta["size"]) > molecular_max_size
        ]
        if len(fallback) > 0:
            largest_fallback = max(fallback, key=lambda cid: int(component_meta[cid]["size"]))
            component_meta[largest_fallback]["is_network_like"] = True
            network_component_ids.add(int(largest_fallback))

    main_network_component_id: int | None = None
    if len(network_component_ids) > 0:
        main_network_component_id = int(
            max(network_component_ids, key=lambda cid: int(component_meta[cid]["size"]))
        )

    return {
        "atom_to_component": atom_to_component,
        "component_meta": component_meta,
        "main_network_component_id": main_network_component_id,
    }


def get_component_status_for_atom(i: int, component_context: Mapping[str, Any]) -> dict[str, Any]:
    """Return component membership/status fields for one atom."""
    atom_to_component = component_context["atom_to_component"]
    component_meta = component_context["component_meta"]

    if int(i) not in atom_to_component:
        return {
            "component_id": int(i),
            "component_size": 1,
            "component_element_counts": {},
            "is_network_bound": False,
            "is_molecular_fragment": True,
        }

    cid = int(atom_to_component[int(i)])
    meta = dict(component_meta[cid])
    is_network_bound = bool(meta.get("is_network_like", False))

    return {
        "component_id": cid,
        "component_size": int(meta.get("size", 1)),
        "component_element_counts": dict(meta.get("element_counts", {})),
        "is_network_bound": is_network_bound,
        "is_molecular_fragment": bool(not is_network_bound),
    }


def motif_effect(is_closed_shell_like: bool, is_radical_like: bool, is_network_bound: bool) -> str:
    """Map motif flags onto a compact chemistry role label."""
    if bool(is_radical_like):
        return "open_shell_candidate"
    if bool(is_closed_shell_like):
        return "quenching"
    if bool(is_network_bound):
        return "modifying"
    return "neutral"


def detect_oxygen_motif(
    atoms: Atoms,
    G: nx.Graph,
    i: int,
    neighbors: Sequence[int],
    neighbor_counts_by_element: Mapping[str, int],
    bond_length_map: Mapping[int, float],
    component_status: Mapping[str, Any],
    config: Mapping[str, Any],
) -> dict[str, Any]:
    """Detect oxygen-centered motifs (hydroxyl/ether/epoxy/carbonyl/etc.).

    Context: Oxygen can either quench spins (closed-shell functional groups) or
    occasionally host open-shell character, so explicit motif labels avoid over-seeding.
    Labels such as ``hydroxide_like_O`` denote OH-fragment-like local connectivity
    and do not imply a formal anionic charge assignment.
    """
    coord = int(len(neighbors))
    c_neighbors = [int(j) for j in neighbors if str(atoms[j].symbol) == "C"]
    h_neighbors = [int(j) for j in neighbors if str(atoms[j].symbol) == "H"]
    o_neighbors = [int(j) for j in neighbors if str(atoms[j].symbol) == "O"]

    mcfg = config.get("motif_detection", {})
    component_size = int(component_status["component_size"])
    component_elements = dict(component_status["component_element_counts"])
    is_network_bound = bool(component_status["is_network_bound"])
    is_molecular_fragment = bool(component_status["is_molecular_fragment"])

    motif_label = "oxygen_other"
    motif_tags: list[str] = []
    is_closed_shell_like = False
    is_radical_like = False
    can_host_initial_moment = False

    if coord == 0:
        motif_label = "undercoordinated_oxygen"
        is_radical_like = True
        can_host_initial_moment = True
    elif (
        coord == 2
        and len(h_neighbors) == 2
        and component_size == 3
        and int(component_elements.get("O", 0)) == 1
        and int(component_elements.get("H", 0)) == 2
    ):
        motif_label = "water_like_O"
        motif_tags.extend(["o_h", "water_like"])
        is_closed_shell_like = True
    elif (
        coord == 1
        and len(h_neighbors) == 1
        and is_molecular_fragment
        and component_size <= 2
        and int(component_elements.get("O", 0)) == 1
        and int(component_elements.get("H", 0)) == 1
    ):
        motif_label = "hydroxide_like_O"
        motif_tags.append("o_h")
        is_closed_shell_like = True
    elif (
        coord == 1
        and len(o_neighbors) == 1
        and component_size == 2
        and int(component_elements.get("O", 0)) == 2
    ):
        motif_tags.append("o_o")
        oo_length = bond_length_map.get(o_neighbors[0], get_bond_length_between(G, i, o_neighbors[0]))
        if oo_length is not None and oo_length <= float(mcfg["o2_oo_max"]):
            motif_label = "molecular_oxygen_like_O"
            is_radical_like = True
            can_host_initial_moment = True
        else:
            motif_label = "peroxide_like_O"
            is_closed_shell_like = True
    elif len(o_neighbors) > 0:
        motif_label = "peroxide_like_O"
        motif_tags.append("o_o")
        is_closed_shell_like = True
    elif coord == 2 and len(c_neighbors) == 1 and len(h_neighbors) == 1:
        motif_label = "hydroxyl_like_O"
        motif_tags.append("o_h")
        is_closed_shell_like = True
    elif coord == 2 and len(c_neighbors) == 2:
        c1, c2 = c_neighbors
        motif_label = "epoxy_like_O" if (is_network_bound and G.has_edge(c1, c2)) else "ether_like_O"
        is_closed_shell_like = True
    elif coord == 1 and len(c_neighbors) == 1:
        c = c_neighbors[0]
        co_length = bond_length_map.get(c, get_bond_length_between(G, i, c))
        c_o_neighbors = [j for j in get_neighbors(G, c) if str(atoms[j].symbol) == "O"]
        c_c_neighbors = [j for j in get_neighbors(G, c) if str(atoms[j].symbol) == "C"]
        has_oh_partner = any(
            int(count_neighbors_by_element(atoms, G, o_idx).get("H", 0)) > 0
            for o_idx in c_o_neighbors
            if int(o_idx) != int(i)
        )
        if (
            len(c_o_neighbors) >= 2
            and (len(c_c_neighbors) > 0 or has_oh_partner)
            and (co_length is None or co_length <= float(mcfg["carboxyl_co_single_max"]))
        ):
            motif_label = "carboxyl_like_O"
            is_closed_shell_like = True
        elif co_length is not None and co_length <= float(mcfg["carbonyl_co_max"]):
            motif_label = "carbonyl_like_O"
            is_closed_shell_like = True
        else:
            motif_label = "undercoordinated_oxygen"
            is_radical_like = True
            can_host_initial_moment = True
    elif coord <= 1:
        motif_label = "undercoordinated_oxygen"
        is_radical_like = True
        can_host_initial_moment = True

    if int(neighbor_counts_by_element.get("H", 0)) > 0:
        motif_tags.append("o_h")
    if int(neighbor_counts_by_element.get("O", 0)) > 0:
        motif_tags.append("o_o")

    return {
        "motif_label": motif_label,
        "motif_tags": sorted(set(motif_tags)),
        "is_closed_shell_like": bool(is_closed_shell_like),
        "is_radical_like": bool(is_radical_like),
        "can_host_initial_moment": bool(can_host_initial_moment),
    }


def detect_nitrogen_motif(
    atoms: Atoms,
    G: nx.Graph,
    i: int,
    neighbors: Sequence[int],
    neighbor_counts_by_element: Mapping[str, int],
    bond_length_map: Mapping[int, float],
    component_status: Mapping[str, Any],
    ring_sizes: Sequence[int],
    hybridization: str,
    config: Mapping[str, Any],
) -> dict[str, Any]:
    """Detect nitrogen-centered motifs (graphitic/pyridinic/pyrrolic/etc.).

    Context: Nitrogen can behave as a dopant, a closed-shell substituent, or an
    open-shell center depending on local bonding, so it needs explicit environment logic.
    ``graphitic``/``pyridinic_like``/``pyrrolic_like`` labels are topology-guided
    local categories rather than strict aromaticity proofs.
    """
    coord = int(len(neighbors))
    c_neighbors = [int(j) for j in neighbors if str(atoms[j].symbol) == "C"]
    h_neighbors = [int(j) for j in neighbors if str(atoms[j].symbol) == "H"]
    n_neighbors = [int(j) for j in neighbors if str(atoms[j].symbol) == "N"]
    o_neighbors = [int(j) for j in neighbors if str(atoms[j].symbol) == "O"]

    mcfg = config.get("motif_detection", {})
    component_size = int(component_status["component_size"])
    component_elements = dict(component_status["component_element_counts"])
    is_network_bound = bool(component_status["is_network_bound"])

    motif_label = "nitrogen_other"
    motif_tags: list[str] = []
    is_closed_shell_like = False
    is_radical_like = False
    can_host_initial_moment = False

    if coord == 0:
        motif_label = "undercoordinated_nitrogen"
        is_radical_like = True
        can_host_initial_moment = True
    elif (
        coord == 1
        and len(n_neighbors) == 1
        and component_size == 2
        and int(component_elements.get("N", 0)) == 2
    ):
        nn_length = bond_length_map.get(n_neighbors[0], get_bond_length_between(G, i, n_neighbors[0]))
        motif_label = "molecular_nitrogen_like_N"
        if nn_length is not None and nn_length <= float(mcfg["n2_nn_max"]):
            is_closed_shell_like = True
        else:
            is_radical_like = True
            can_host_initial_moment = True
    elif coord == 1 and len(c_neighbors) == 1:
        cn_length = bond_length_map.get(c_neighbors[0], get_bond_length_between(G, i, c_neighbors[0]))
        if cn_length is not None and cn_length <= float(mcfg["nitrile_cn_max"]):
            motif_label = "nitrile_like_N"
            is_closed_shell_like = True
        else:
            motif_label = "undercoordinated_nitrogen"
            is_radical_like = True
            can_host_initial_moment = True
    elif coord == 2:
        if len(c_neighbors) == 2 and is_network_bound and 5 in set(map(int, ring_sizes)):
            motif_label = "pyridinic_like_N"
            can_host_initial_moment = True
        elif len(c_neighbors) >= 1:
            motif_label = "imine_like_N"
            can_host_initial_moment = True
        elif len(h_neighbors) >= 1:
            motif_label = "amine_like_N"
            is_closed_shell_like = True
    elif coord == 3:
        if len(c_neighbors) >= 3 and is_network_bound and is_sp2_like_label(hybridization):
            motif_label = "graphitic_N"
            can_host_initial_moment = True
        elif len(c_neighbors) >= 2 and len(h_neighbors) >= 1 and 5 in set(map(int, ring_sizes)):
            motif_label = "pyrrolic_like_N"
            motif_tags.append("n_h_containing")
            can_host_initial_moment = True
        elif len(h_neighbors) >= 1:
            motif_label = "amine_like_N"
            motif_tags.append("n_h_containing")
            is_closed_shell_like = True
        else:
            motif_label = "network_nitrogen"
            can_host_initial_moment = True

    if len(o_neighbors) > 0:
        motif_tags.append("n_o_containing")
        if motif_label in {"nitrogen_other", "amine_like_N", "imine_like_N"}:
            motif_label = "n_o_containing_N"
            is_closed_shell_like = True
            can_host_initial_moment = False
    if len(h_neighbors) > 0:
        motif_tags.append("n_h_containing")

    if coord <= 1 and motif_label != "molecular_nitrogen_like_N":
        motif_label = "undercoordinated_nitrogen"
        is_closed_shell_like = False
        is_radical_like = True
        can_host_initial_moment = True

    return {
        "motif_label": motif_label,
        "motif_tags": sorted(set(motif_tags)),
        "is_closed_shell_like": bool(is_closed_shell_like),
        "is_radical_like": bool(is_radical_like),
        "can_host_initial_moment": bool(can_host_initial_moment),
    }


def detect_hydrogen_motif(
    atoms: Atoms,
    G: nx.Graph,
    i: int,
    neighbors: Sequence[int],
    bond_length_map: Mapping[int, float],
    component_status: Mapping[str, Any],
    config: Mapping[str, Any],
) -> dict[str, Any]:
    """Detect hydrogen roles (terminating H, O-H/N-H, H2-like, adsorbed H).

    Context: Hydrogen is often a terminator that suppresses radical character,
    but isolated or unusual H motifs can still indicate open-shell fragments.
    """
    coord = int(len(neighbors))
    mcfg = config.get("motif_detection", {})
    component_size = int(component_status["component_size"])
    component_elements = dict(component_status["component_element_counts"])
    is_molecular_fragment = bool(component_status["is_molecular_fragment"])

    motif_label = "hydrogen_other"
    motif_tags: list[str] = []
    is_closed_shell_like = False
    is_radical_like = False
    can_host_initial_moment = False

    if coord == 0:
        motif_label = "isolated_hydrogen"
        is_radical_like = True
        can_host_initial_moment = True
    elif coord == 1:
        j = int(neighbors[0])
        j_symbol = str(atoms[j].symbol)

        if j_symbol == "H":
            hh_length = bond_length_map.get(j, get_bond_length_between(G, i, j))
            if (
                component_size == 2
                and int(component_elements.get("H", 0)) == 2
                and hh_length is not None
                and hh_length <= float(mcfg["h2_hh_max"])
            ):
                motif_label = "h2_like_H"
                is_closed_shell_like = True
            else:
                motif_label = "isolated_hydrogen"
                is_radical_like = True
                can_host_initial_moment = True
        elif j_symbol == "O":
            o_counts = count_neighbors_by_element(atoms, G, j)
            o_coord = len(get_neighbors(G, j))
            motif_tags.append("o_h")
            if (
                o_coord == 2
                and int(o_counts.get("H", 0)) == 2
                and component_size == 3
                and int(component_elements.get("O", 0)) == 1
                and int(component_elements.get("H", 0)) == 2
            ):
                motif_label = "water_hydrogen"
            elif int(o_counts.get("C", 0)) >= 1:
                motif_label = "hydroxyl_hydrogen"
            else:
                motif_label = "o_h_hydrogen"
            is_closed_shell_like = True
        elif j_symbol == "N":
            motif_label = "n_h_hydrogen"
            motif_tags.append("n_h")
            is_closed_shell_like = True
        elif j_symbol == "C":
            c_coord = len(get_neighbors(G, j))
            if c_coord <= 2:
                motif_label = "c_dangling_terminator_H"
                is_closed_shell_like = True
            elif is_molecular_fragment:
                motif_label = "hydrocarbon_fragment_H"
                is_closed_shell_like = True
            else:
                motif_label = "adsorbed_hydrogen"
                can_host_initial_moment = bool(c_coord <= 1)
                is_radical_like = bool(can_host_initial_moment)
        else:
            motif_label = "adsorbed_hydrogen"
            is_radical_like = True
            can_host_initial_moment = True
    else:
        motif_label = "multi_bonded_hydrogen"

    return {
        "motif_label": motif_label,
        "motif_tags": sorted(set(motif_tags)),
        "is_closed_shell_like": bool(is_closed_shell_like),
        "is_radical_like": bool(is_radical_like),
        "can_host_initial_moment": bool(can_host_initial_moment),
    }


def detect_carbon_motif(
    atoms: Atoms,
    G: nx.Graph,
    i: int,
    neighbors: Sequence[int],
    neighbor_counts_by_element: Mapping[str, int],
    component_status: Mapping[str, Any],
    config: Mapping[str, Any],
) -> dict[str, Any]:
    """Detect carbon motifs, including gas-like fragments vs network-bound sites.

    Context: Carbon in the main network and carbon in detached molecules should
    not be scored the same, so this tags fragment chemistry separately from network sites.
    """
    coord = int(len(neighbors))
    c_neighbors = [int(j) for j in neighbors if str(atoms[j].symbol) == "C"]
    h_neighbors = [int(j) for j in neighbors if str(atoms[j].symbol) == "H"]
    o_neighbors = [int(j) for j in neighbors if str(atoms[j].symbol) == "O"]
    n_neighbors = [int(j) for j in neighbors if str(atoms[j].symbol) == "N"]

    mcfg = config.get("motif_detection", {})
    component_size = int(component_status["component_size"])
    component_elements = dict(component_status["component_element_counts"])
    is_network_bound = bool(component_status["is_network_bound"])
    is_molecular_fragment = bool(component_status["is_molecular_fragment"])

    motif_label = "carbon_other"
    motif_tags: list[str] = []
    is_closed_shell_like = False
    is_radical_like = False
    can_host_initial_moment = True

    if is_molecular_fragment:
        if (
            component_size == 2
            and int(component_elements.get("C", 0)) == 1
            and int(component_elements.get("O", 0)) == 1
            and coord == 1
            and len(o_neighbors) == 1
            and len(get_neighbors(G, o_neighbors[0])) == 1
        ):
            motif_label = "co_like_C"
            is_closed_shell_like = True
            can_host_initial_moment = False
        elif (
            component_size == 3
            and int(component_elements.get("C", 0)) == 1
            and int(component_elements.get("O", 0)) == 2
            and coord == 2
            and len(o_neighbors) == 2
            and all(len(get_neighbors(G, o)) == 1 for o in o_neighbors)
        ):
            motif_label = "co2_like_C"
            is_closed_shell_like = True
            can_host_initial_moment = False
        elif len(h_neighbors) == 1:
            motif_label = "ch_like_fragment_C"
            is_radical_like = bool(coord < 4)
            is_closed_shell_like = bool(coord >= 4)
            can_host_initial_moment = bool(is_radical_like)
        elif len(h_neighbors) == 2:
            motif_label = "ch2_like_fragment_C"
            is_radical_like = bool(coord < 4)
            is_closed_shell_like = bool(coord >= 4)
            can_host_initial_moment = bool(is_radical_like)
        elif len(h_neighbors) >= 3:
            motif_label = "ch3_like_fragment_C"
            is_radical_like = bool(coord < 4)
            is_closed_shell_like = bool(coord >= 4)
            can_host_initial_moment = bool(is_radical_like)
        elif (
            set(component_elements.keys()).issubset({"C", "H"})
            and component_size <= int(mcfg["small_hydrocarbon_max_size"])
        ):
            motif_label = "small_hydrocarbon_fragment_C"
            is_radical_like = bool(coord < 4)
            is_closed_shell_like = bool(coord >= 4)
            can_host_initial_moment = bool(is_radical_like)
        elif coord <= 2 and len(o_neighbors) == 0 and len(n_neighbors) == 0:
            motif_label = "undercoordinated_fragment_C"
            is_radical_like = True
            can_host_initial_moment = True
        else:
            motif_label = "molecular_fragment_C"
            is_closed_shell_like = bool(coord >= 4)
            is_radical_like = bool(coord <= 2)
            can_host_initial_moment = bool(is_radical_like)
    elif is_network_bound:
        if int(neighbor_counts_by_element.get("O", 0)) > 0 or int(neighbor_counts_by_element.get("N", 0)) > 0:
            motif_label = "heteroatom_modified_network_C"
        elif int(neighbor_counts_by_element.get("H", 0)) > 0:
            motif_label = "hydrogenated_network_C"
            is_closed_shell_like = bool(coord >= 3)
        else:
            motif_label = "network_carbon"
        can_host_initial_moment = True
    else:
        motif_label = "carbon_other"

    if len(o_neighbors) > 0:
        motif_tags.append("c_o")
    if len(n_neighbors) > 0:
        motif_tags.append("c_n")
    if len(h_neighbors) > 0:
        motif_tags.append("c_h")
    if len(c_neighbors) > 0:
        motif_tags.append("c_c")

    return {
        "motif_label": motif_label,
        "motif_tags": sorted(set(motif_tags)),
        "is_closed_shell_like": bool(is_closed_shell_like),
        "is_radical_like": bool(is_radical_like),
        "can_host_initial_moment": bool(can_host_initial_moment),
    }


def detect_generic_motif(
    atoms: Atoms,
    i: int,
    neighbors: Sequence[int],
) -> dict[str, Any]:
    """Fallback motif detection for elements beyond C/H/O/N."""
    element = str(atoms[i].symbol)
    coord = int(len(neighbors))
    if coord == 0:
        return {
            "motif_label": f"isolated_{element}",
            "motif_tags": [],
            "is_closed_shell_like": False,
            "is_radical_like": True,
            "can_host_initial_moment": True,
        }
    return {
        "motif_label": f"{element}_other",
        "motif_tags": [],
        "is_closed_shell_like": False,
        "is_radical_like": False,
        "can_host_initial_moment": False,
    }


def classify_atom_motif(
    atoms: Atoms,
    G: nx.Graph,
    i: int,
    hybridization: str,
    ring_sizes: Sequence[int],
    component_context: Mapping[str, Any],
    config: Mapping[str, Any],
) -> dict[str, Any]:
    """Classify atom motif and chemistry role (closed-shell/modifier/open-shell).

    Context: This is the single dispatch point that harmonizes element-specific
    motif detectors into shared flags used by feature building and magnetic scoring.
    """
    element = str(atoms[i].symbol)
    neighbors = get_neighbors(G, i)
    neighbor_counts = count_neighbors_by_element(atoms, G, i)
    bond_length_map = get_neighbor_bond_length_map(G, i)
    component_status = get_component_status_for_atom(i, component_context)

    if element == "O":
        motif = detect_oxygen_motif(
            atoms,
            G,
            i,
            neighbors,
            neighbor_counts,
            bond_length_map,
            component_status,
            config,
        )
    elif element == "N":
        motif = detect_nitrogen_motif(
            atoms,
            G,
            i,
            neighbors,
            neighbor_counts,
            bond_length_map,
            component_status,
            ring_sizes,
            hybridization,
            config,
        )
    elif element == "H":
        motif = detect_hydrogen_motif(
            atoms,
            G,
            i,
            neighbors,
            bond_length_map,
            component_status,
            config,
        )
    elif element == str(config["carbon_symbol"]):
        motif = detect_carbon_motif(
            atoms,
            G,
            i,
            neighbors,
            neighbor_counts,
            component_status,
            config,
        )
    else:
        motif = detect_generic_motif(atoms, i, neighbors)

    motif["is_network_bound"] = bool(component_status["is_network_bound"])
    motif["is_molecular_fragment"] = bool(component_status["is_molecular_fragment"])
    motif["component_id"] = int(component_status["component_id"])
    motif["component_size"] = int(component_status["component_size"])
    motif["component_element_counts"] = dict(component_status["component_element_counts"])
    motif["motif_effect"] = motif_effect(
        bool(motif.get("is_closed_shell_like", False)),
        bool(motif.get("is_radical_like", False)),
        bool(component_status["is_network_bound"]),
    )
    return motif
