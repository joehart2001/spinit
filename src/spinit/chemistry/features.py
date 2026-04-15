"""Atom-level feature construction for magnetic seeding."""

from __future__ import annotations

from collections import Counter
from typing import Any, Mapping, Sequence

import networkx as nx
from ase import Atoms

from ..graph.geometry import get_bond_angles
from ..graph.graph_ops import count_neighbors_by_element, get_bond_lengths, get_neighbors
from ..graph.topology import compute_topological_defect_score, get_atom_ring_sizes
from .hybridization import (
    classify_all_hybridizations,
    detect_pi_unsaturated_site,
    detect_sigma_dangling_bond,
    estimate_sp2_normal,
    evaluate_sp_carbon_environment,
    is_strained_sp2_site,
    is_undercoordinated_carbon,
    measure_sp2_pi_alignment,
    measure_sp2_pyramidalization,
)
from .motifs import build_component_context, classify_atom_motif


def compute_element_specific_features(
    atoms: Atoms,
    G: nx.Graph,
    i: int,
    hybridization: str,
    bond_angles: Sequence[float],
    neighbor_counts_by_element: Mapping[str, int],
    ring_sizes: Sequence[int],
    hybridizations: Mapping[int, str],
    motif_info: Mapping[str, Any],
    config: Mapping[str, Any],
) -> dict[str, Any]:
    """Compute element-aware site features without assuming specific heteroatoms.

    Context: This keeps carbon-centered logic while safely handling mixed-element
    systems so non-carbon atoms are described without carbon-only assumptions.
    """
    element = str(atoms[i].symbol)
    coordination = len(get_neighbors(G, i))

    carbon_symbol = str(config["carbon_symbol"])
    allow_n_weak = bool(config["elements"].get("allow_n_weak_candidate", False))
    sp_flags = evaluate_sp_carbon_environment(
        coordination=coordination,
        hybridization=hybridization,
        bond_angles=bond_angles,
        neighbor_counts_by_element=neighbor_counts_by_element,
        config=config,
    )

    if element == carbon_symbol:
        undercoordinated = is_undercoordinated_carbon(coordination, hybridization)
        sigma_dangling = detect_sigma_dangling_bond(coordination)
        pi_unsaturated = detect_pi_unsaturated_site(coordination, hybridization, bond_angles, config)
        strained_sp2 = is_strained_sp2_site(G, i, hybridizations, config)

        if bool(sp_flags["is_sp_site"]) and not bool(sp_flags["sp_radical_like"]):
            undercoordinated = False
            pi_unsaturated = False
    elif element == "N":
        undercoordinated = coordination < 3
        sigma_dangling = detect_sigma_dangling_bond(coordination)
        pi_unsaturated = detect_pi_unsaturated_site(coordination, hybridization, bond_angles, config)
        strained_sp2 = is_strained_sp2_site(G, i, hybridizations, config)
    elif element == "O":
        undercoordinated = bool(coordination <= 1 and motif_info.get("is_radical_like", False))
        sigma_dangling = bool(coordination <= 1 and motif_info.get("is_radical_like", False))
        pi_unsaturated = False
        strained_sp2 = False
    elif element == "H":
        undercoordinated = bool(coordination == 0)
        sigma_dangling = bool(coordination == 0)
        pi_unsaturated = False
        strained_sp2 = False
    else:
        undercoordinated = bool(coordination <= 1 and motif_info.get("is_radical_like", False))
        sigma_dangling = bool(coordination <= 1 and motif_info.get("is_radical_like", False))
        pi_unsaturated = False
        strained_sp2 = False

    is_candidate_element = bool(motif_info.get("can_host_initial_moment", False))
    if element == "N" and (not allow_n_weak) and (not bool(motif_info.get("is_radical_like", False))):
        is_candidate_element = False

    topological_defect_score = compute_topological_defect_score(ring_sizes, element=element, config=config)

    return {
        "is_candidate_element": bool(is_candidate_element),
        "is_undercoordinated": bool(undercoordinated),
        "sigma_dangling_bond": bool(sigma_dangling),
        "pi_unsaturated": bool(pi_unsaturated),
        "is_strained_sp2": bool(strained_sp2),
        "topological_defect_score": float(topological_defect_score),
        **sp_flags,
    }


def compute_atom_features(
    atoms: Atoms,
    G: nx.Graph,
    ring_sizes_per_atom: Sequence[Sequence[int]],
    hybridizations: Mapping[int, str],
    component_context: Mapping[str, Any],
    i: int,
    config: Mapping[str, Any],
) -> dict[str, Any]:
    """Compute plain-dictionary feature bundle for one atom.

    Context: This assembles geometry, topology, motif, and chemistry indicators
    into one inspectable record that downstream scoring can use transparently.
    """
    element = str(atoms[i].symbol)
    neighbors = get_neighbors(G, i)
    coordination = len(neighbors)
    neighbor_counts_by_element = count_neighbors_by_element(atoms, G, i)

    bond_lengths = get_bond_lengths(G, i)
    bond_angles = get_bond_angles(G, i)
    hybridization = hybridizations.get(i, "unknown")

    ring_sizes = get_atom_ring_sizes(ring_sizes_per_atom, i)
    normal = estimate_sp2_normal(G, i)
    pyramidalization = measure_sp2_pyramidalization(G, i)
    pi_alignment = measure_sp2_pi_alignment(G, i, hybridizations)
    motif_info = classify_atom_motif(
        atoms,
        G,
        i,
        hybridization,
        ring_sizes,
        component_context,
        config,
    )

    element_features = compute_element_specific_features(
        atoms,
        G,
        i,
        hybridization,
        bond_angles,
        neighbor_counts_by_element,
        ring_sizes,
        hybridizations,
        motif_info,
        config,
    )

    return {
        "atom_index": int(i),
        "element": element,
        "neighbors": neighbors,
        "coordination": int(coordination),
        "neighbor_counts_by_element": neighbor_counts_by_element,
        "bond_lengths": [float(x) for x in bond_lengths],
        "bond_angles": [float(x) for x in bond_angles],
        "hybridization": hybridization,
        "ring_sizes": ring_sizes,
        "sp2_normal": None if normal is None else [float(x) for x in normal],
        "pyramidalization": None if pyramidalization is None else float(pyramidalization),
        "pi_alignment": None if pi_alignment is None else float(pi_alignment),
        "motif_label": str(motif_info["motif_label"]),
        "motif_tags": list(motif_info.get("motif_tags", [])),
        "is_molecular_fragment": bool(motif_info["is_molecular_fragment"]),
        "is_network_bound": bool(motif_info["is_network_bound"]),
        "is_closed_shell_like": bool(motif_info["is_closed_shell_like"]),
        "is_radical_like": bool(motif_info["is_radical_like"]),
        "can_host_initial_moment": bool(motif_info["can_host_initial_moment"]),
        "motif_effect": str(motif_info.get("motif_effect", "neutral")),
        "component_id": int(motif_info["component_id"]),
        "component_size": int(motif_info["component_size"]),
        "component_element_counts": dict(motif_info.get("component_element_counts", {})),
        **element_features,
        "neighbor_motif_counts": {},
        "intrinsic_score": 0.0,
        "neighbor_modifier": 0.0,
        "magnetic_score": 0.0,
    }


def compute_all_atom_features(
    atoms: Atoms,
    G: nx.Graph,
    ring_sizes_per_atom: Sequence[Sequence[int]],
    config: Mapping[str, Any],
) -> dict[int, dict[str, Any]]:
    """Compute feature dictionary for all atoms, keyed by atom index.

    Context: This is the main feature-construction pass used before scoring, and
    it also adds neighbor motif summaries needed for local environment modifiers.
    """
    hybridizations = classify_all_hybridizations(G, config)
    component_context = build_component_context(atoms, G, config)

    feature_dict: dict[int, dict[str, Any]] = {}
    for i in range(len(atoms)):
        feature_dict[i] = compute_atom_features(
            atoms,
            G,
            ring_sizes_per_atom,
            hybridizations,
            component_context,
            i,
            config,
        )

    for i in range(len(atoms)):
        counts = Counter(
            str(feature_dict[j].get("motif_label", "unknown"))
            for j in get_neighbors(G, i)
            if j in feature_dict
        )
        feature_dict[i]["neighbor_motif_counts"] = {
            label: int(count) for label, count in sorted(counts.items())
        }

    return feature_dict
