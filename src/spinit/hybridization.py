"""Hybridization-like local chemistry and sp2 strain utilities."""

from __future__ import annotations

from itertools import combinations
from typing import Any, Mapping, Sequence

import networkx as nx
import numpy as np

from .geometry import angle_between, get_bond_angles, mean_abs_deviation, unit_vector
from .graph_ops import get_bond_vectors, get_neighbors


def classify_hybridization(
    coordination: int,
    bond_angles: Sequence[float],
    config: Mapping[str, Any],
) -> str:
    """Classify local environment into simple hybridization-like labels.

    Context: This gives a lightweight chemical-state label that downstream
    unsaturation and strain checks use when deciding where initial spin seeds
    are chemically plausible.
    """
    hcfg = config["hybridization"]

    if coordination <= 0:
        return "isolated"
    if coordination == 1:
        return "undercoordinated"

    sp_target = float(hcfg["sp_linear_target_deg"])
    sp2_target = float(hcfg["sp2_target_deg"])
    sp3_target = float(hcfg["sp3_target_deg"])

    if coordination == 2:
        if len(bond_angles) == 0:
            return "undercoordinated_distorted"
        mean_angle = float(np.mean(bond_angles))
        if abs(mean_angle - sp_target) <= float(hcfg["sp_linear_tol_deg"]):
            return "sp"
        if abs(mean_angle - sp2_target) <= float(hcfg["sp2_tol_deg"]):
            return "undercoordinated_sp2_like"
        return "undercoordinated_distorted"

    if coordination == 3:
        mad_sp2 = mean_abs_deviation(bond_angles, sp2_target)
        mad_sp3 = mean_abs_deviation(bond_angles, sp3_target)

        if mad_sp2 is None:
            return "sp2_distorted"
        if mad_sp2 <= float(hcfg["sp2_mad_max_deg"]):
            return "sp2"
        if mad_sp3 is not None and mad_sp3 <= float(hcfg["sp3_mad_max_deg"]):
            return "sp3_like"
        if mad_sp3 is None or mad_sp2 <= mad_sp3:
            return "sp2_distorted"
        return "sp3_like_distorted"

    if coordination == 4:
        return "sp3"

    return "unknown"


def classify_all_hybridizations(G: nx.Graph, config: Mapping[str, Any]) -> dict[int, str]:
    """Classify hybridization-like labels for all graph nodes.

    Context: The full-node map is reused across feature extraction so we apply
    one consistent hybridization heuristic everywhere in the pipeline.
    """
    labels: dict[int, str] = {}
    for i in G.nodes:
        i = int(i)
        coordination = len(get_neighbors(G, i))
        labels[i] = classify_hybridization(coordination, get_bond_angles(G, i), config)
    return labels


def is_undercoordinated_carbon(coordination: int, hybridization: str) -> bool:
    """Return True for carbon undercoordination patterns."""
    if coordination <= 1:
        return True
    if coordination == 2 and hybridization == "sp":
        return False
    return coordination < 3 or hybridization.startswith("undercoordinated")


def detect_sigma_dangling_bond(coordination: int) -> bool:
    """Return True for sigma dangling-bond-like patterns."""
    return coordination <= 1


def detect_pi_unsaturated_site(
    coordination: int,
    hybridization: str,
    bond_angles: Sequence[float],
    config: Mapping[str, Any],
) -> bool:
    """Return True for pi-unsaturated / localized-spin-like environments.

    Context: This flags sites where pi bonding is likely incomplete or distorted,
    which makes spin localization more likely than in well-satisfied closed-shell sites.
    """
    if coordination <= 1:
        return True
    if coordination == 2:
        if hybridization == "sp":
            if len(bond_angles) == 0:
                return False
            sp_cfg = config.get("sp_site", {})
            linear_target = float(config["hybridization"]["sp_linear_target_deg"])
            tol = float(sp_cfg.get("linearity_tolerance_deg", config["hybridization"]["sp_linear_tol_deg"]))
            mean_angle = float(np.mean(bond_angles))
            return abs(mean_angle - linear_target) > tol
        return True
    if coordination != 3:
        return False

    if hybridization in {"sp2_distorted", "undercoordinated_sp2_like"}:
        mad120 = mean_abs_deviation(bond_angles, float(config["hybridization"]["sp2_target_deg"]))
        if mad120 is None:
            return True
        return mad120 >= float(config["pi_unsaturated"]["angle_mad_min_deg"])

    return False


def evaluate_sp_carbon_environment(
    coordination: int,
    hybridization: str,
    bond_angles: Sequence[float],
    neighbor_counts_by_element: Mapping[str, int],
    config: Mapping[str, Any],
) -> dict[str, Any]:
    """Classify whether a carbon sp site looks heteroatom-modified vs radical-like.

    Context: Not all sp carbons should seed moments; this function separates
    chemically satisfied sp motifs from genuinely unsatisfied radical-like cases.
    """
    is_sp_site = bool(coordination == 2 and hybridization == "sp")
    if not is_sp_site:
        return {
            "is_sp_site": False,
            "sp_is_heteroatom_modified": False,
            "sp_is_linear": False,
            "sp_is_distorted": False,
            "sp_radical_like": False,
            "sp_mean_angle_deg": None,
        }

    modifier_elements = set(map(str, config.get("elements", {}).get("heteroatom_modifier_elements", [])))
    sp_cfg = config.get("sp_site", {})
    linear_target = float(config["hybridization"]["sp_linear_target_deg"])
    tol = float(sp_cfg.get("linearity_tolerance_deg", config["hybridization"]["sp_linear_tol_deg"]))

    mean_angle = float(np.mean(bond_angles)) if len(bond_angles) > 0 else None
    is_linear = bool(mean_angle is not None and abs(mean_angle - linear_target) <= tol)
    is_distorted = bool(mean_angle is None or not is_linear)
    is_heteroatom_modified = any(
        int(neighbor_counts_by_element.get(elem, 0)) > 0 for elem in modifier_elements
    )

    allow_linear_unmodified = bool(sp_cfg.get("allow_unmodified_linear_sp_as_radical_like", False))
    radical_like = bool((not is_heteroatom_modified) and (is_distorted or allow_linear_unmodified))

    return {
        "is_sp_site": True,
        "sp_is_heteroatom_modified": bool(is_heteroatom_modified),
        "sp_is_linear": bool(is_linear),
        "sp_is_distorted": bool(is_distorted),
        "sp_radical_like": bool(radical_like),
        "sp_mean_angle_deg": None if mean_angle is None else float(mean_angle),
    }


def is_sp2_like_label(label: str) -> bool:
    """Return True for sp2-like labels used in pi-overlap analysis."""
    return label in {"sp2", "sp2_distorted", "undercoordinated_sp2_like"}


def estimate_sp2_normal(G: nx.Graph, i: int) -> np.ndarray | None:
    """Estimate local pz-like normal for a 3-fold site.

    Context: This approximates local p-orbital orientation at sp2-like sites so
    we can detect when neighboring pz directions are misaligned and pi overlap is reduced.
    """
    vectors = get_bond_vectors(G, i)
    if len(vectors) != 3:
        return None

    unit_bonds = [unit_vector(v) for v in vectors]
    cross_normals: list[np.ndarray] = []
    for a, b in combinations(unit_bonds, 2):
        n = np.cross(a, b)
        if np.linalg.norm(n) < 1e-10:
            continue
        cross_normals.append(unit_vector(n))

    if len(cross_normals) == 0:
        return None

    reference = cross_normals[0]
    aligned: list[np.ndarray] = []
    for n in cross_normals:
        if np.dot(n, reference) < 0.0:
            n = -n
        aligned.append(n)

    normal = np.sum(aligned, axis=0)
    if np.linalg.norm(normal) < 1e-10:
        return None
    return unit_vector(normal)


def measure_sp2_pyramidalization(G: nx.Graph, i: int) -> float | None:
    """Measure deviation from ideal planar sp2 geometry in degrees.

    Context: Pyramidalization quantifies how far a nominally sp2 site bends out
    of plane, which is a practical strain signal for weakened delocalized pi bonding.
    """
    normal = estimate_sp2_normal(G, i)
    if normal is None:
        return None

    vectors = get_bond_vectors(G, i)
    if len(vectors) != 3:
        return None

    deviations: list[float] = []
    for v in vectors:
        angle_to_normal = angle_between(v, normal, degrees=True)
        if not np.isnan(angle_to_normal):
            deviations.append(abs(angle_to_normal - 90.0))

    if len(deviations) == 0:
        return None
    return float(np.mean(deviations))


def measure_sp2_pi_alignment(
    G: nx.Graph,
    i: int,
    hybridizations: Mapping[int, str],
) -> float | None:
    """Measure mean |n_i . n_j| against neighboring sp2-like normals.

    Context: Lower alignment indicates neighboring pz-like directions are not
    co-oriented, which reduces local pi overlap and supports spin-seed placement.
    """
    if not is_sp2_like_label(hybridizations.get(i, "")):
        return None

    normal_i = estimate_sp2_normal(G, i)
    if normal_i is None:
        return None

    alignments: list[float] = []
    for j in get_neighbors(G, i):
        if not is_sp2_like_label(hybridizations.get(j, "")):
            continue
        normal_j = estimate_sp2_normal(G, j)
        if normal_j is not None:
            alignments.append(abs(float(np.dot(normal_i, normal_j))))

    if len(alignments) == 0:
        return None
    return float(np.mean(alignments))


def is_strained_sp2_site(
    G: nx.Graph,
    i: int,
    hybridizations: Mapping[int, str],
    config: Mapping[str, Any],
) -> bool:
    """Return True when sp2-like local geometry indicates strain/poor pi overlap.

    Context: This combines out-of-plane bending and pi-normal misalignment into
    one conservative strained-site flag for magnetic scoring.
    """
    if not is_sp2_like_label(hybridizations.get(i, "unknown")):
        return False

    pyramidalization = measure_sp2_pyramidalization(G, i)
    pi_alignment = measure_sp2_pi_alignment(G, i, hybridizations)
    strain_cfg = config["sp2_strain"]

    pyramidalization_bad = (
        pyramidalization is not None
        and pyramidalization >= float(strain_cfg["pyramidalization_min_deg"])
    )
    alignment_bad = (
        pi_alignment is not None
        and pi_alignment < float(strain_cfg["pi_alignment_min"])
    )
    return bool(pyramidalization_bad or alignment_bad)
