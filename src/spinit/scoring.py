"""Magnetic scoring and candidate selection."""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np


def get_weights_for_element(element: str, config: Mapping[str, Any]) -> Mapping[str, float]:
    """Return intrinsic score weights for element (zero defaults when missing)."""
    weight_map = config.get("intrinsic_score_weights", {})
    defaults = {
        "undercoordination": 0.0,
        "sigma_dangling": 0.0,
        "pi_unsaturated": 0.0,
        "strained_sp2": 0.0,
        "pi_misalignment": 0.0,
        "topological_defect": 0.0,
    }
    values = dict(weight_map.get(str(element), {}))
    return {key: float(values.get(key, default)) for key, default in defaults.items()}


def compute_motif_intrinsic_modifier(features: Mapping[str, Any], config: Mapping[str, Any]) -> float:
    """Return additive score term from explicit motif labels and motif role flags.

    Context: This lets chemically meaningful motifs override generic coordination
    heuristics, which is important for closed-shell molecules and known open-shell motifs.
    """
    mcfg = config.get("motif_scoring", {})
    label_weights = dict(mcfg.get("label_weights", {}))
    motif_label = str(features.get("motif_label", ""))
    element = str(features.get("element", ""))

    term = float(label_weights.get(motif_label, 0.0))

    if bool(features.get("is_closed_shell_like", False)):
        term += float(mcfg.get("closed_shell_penalty", -1.4))
        if bool(features.get("is_molecular_fragment", False)):
            term += float(mcfg.get("molecular_closed_shell_penalty", -0.6))

    if bool(features.get("is_radical_like", False)):
        term += float(mcfg.get("radical_like_bonus", 2.3))
        if element in {"H", "O", "N"} and bool(features.get("can_host_initial_moment", False)):
            term += float(mcfg.get("heteroatom_host_bonus", 0.7))

    return float(term)


def compute_neighbor_modifier(
    features: Mapping[str, Any],
    config: Mapping[str, Any],
) -> float:
    """Return score modifier from neighboring element chemistry and motif labels.

    Context: Nearby atoms can quench or enhance local spin tendency, so this term
    captures first-neighbor chemical environment effects separately from intrinsic site terms.
    """
    element = str(features.get("element", ""))
    neighbor_counts = dict(features.get("neighbor_counts_by_element", {}))
    neighbor_motif_counts = dict(features.get("neighbor_motif_counts", {}))

    rules = config.get("neighbor_modifiers", {}).get(element, {})
    modifier = 0.0
    for neighbor_element, count in neighbor_counts.items():
        modifier += float(rules.get(str(neighbor_element), 0.0)) * float(count)

    motif_rules = config.get("motif_neighbor_modifiers", {}).get(element, {})
    for motif_label, count in neighbor_motif_counts.items():
        modifier += float(motif_rules.get(str(motif_label), 0.0)) * float(count)

    ctrl = config.get("score_controls", {})
    min_mod = float(ctrl.get("neighbor_modifier_min", -np.inf))
    max_mod = float(ctrl.get("neighbor_modifier_max", np.inf))
    return float(np.clip(modifier, min_mod, max_mod))


def compute_intrinsic_site_score(features: Mapping[str, Any], config: Mapping[str, Any]) -> float:
    """Compute intrinsic (local) site score before neighbor-element modifiers.

    Context: This is the core additive heuristic where undercoordination,
    unsaturation, strain, and motif signals are combined at the atom itself.
    """
    if not bool(features.get("is_candidate_element", False)):
        return 0.0

    element = str(features.get("element", ""))
    weights = get_weights_for_element(element, config)

    sp_cfg = config.get("sp_site", {})
    intrinsic_scale = 1.0
    if bool(features.get("is_sp_site", False)) and bool(sp_cfg.get("require_radical_like_for_intrinsic", True)):
        if not bool(features.get("sp_radical_like", False)):
            intrinsic_scale = float(sp_cfg.get("satisfied_sp_intrinsic_scale", 0.0))

    score = 0.0
    if features.get("is_undercoordinated", False):
        score += float(weights["undercoordination"])
    if features.get("sigma_dangling_bond", False):
        score += float(weights["sigma_dangling"])
    if features.get("pi_unsaturated", False):
        score += float(weights["pi_unsaturated"])
    if features.get("is_strained_sp2", False):
        score += float(weights["strained_sp2"])

    controls = config.get("score_controls", {})

    pi_alignment = features.get("pi_alignment", None)
    if pi_alignment is not None:
        alignment_ref = float(config["sp2_strain"]["pi_alignment_min"])
        misalignment = max(0.0, alignment_ref - float(pi_alignment))
        if alignment_ref > 1e-12:
            misalignment /= alignment_ref
        pi_term = float(weights["pi_misalignment"]) * misalignment
        if features.get("is_strained_sp2", False):
            pi_term *= float(controls.get("pi_misalignment_scale_when_strained_sp2", 0.5))
        score += pi_term

    topo_term = float(weights["topological_defect"]) * float(features.get("topological_defect_score", 0.0))
    if bool(controls.get("topology_requires_plausible_site", True)):
        plausible = bool(
            features.get("pi_unsaturated", False)
            or features.get("is_strained_sp2", False)
            or features.get("hybridization", "") in {"sp2", "sp2_distorted", "undercoordinated_sp2_like"}
        )
        if not plausible:
            topo_term = 0.0
        else:
            topo_term *= float(controls.get("topology_scale_for_plausible_site", 0.5))
    score += topo_term
    score += compute_motif_intrinsic_modifier(features, config)

    score *= intrinsic_scale
    return float(score)


def compute_magnetic_score(features: Mapping[str, Any], config: Mapping[str, Any]) -> float:
    """Compute total site score = intrinsic site score + neighbor modifier."""
    intrinsic = compute_intrinsic_site_score(features, config)
    modifier = compute_neighbor_modifier(features, config)
    return float(intrinsic + modifier)


def assign_all_magnetic_scores(
    feature_dict: Mapping[int, dict[str, Any]],
    config: Mapping[str, Any],
) -> dict[int, float]:
    """Assign intrinsic, modifier, and total magnetic score to all atoms.

    Context: Storing score components per atom makes the seeding decisions
    auditable and easier to tune without changing model structure.
    """
    scores: dict[int, float] = {}
    for i, features in feature_dict.items():
        intrinsic = compute_intrinsic_site_score(features, config)
        modifier = compute_neighbor_modifier(features, config)
        total = float(intrinsic + modifier)

        features["intrinsic_score"] = float(intrinsic)
        features["neighbor_modifier"] = float(modifier)
        features["magnetic_score"] = float(total)
        scores[int(i)] = float(total)
    return scores


def select_magnetic_candidate_sites(
    feature_dict: Mapping[int, Mapping[str, Any]],
    config: Mapping[str, Any],
) -> list[int]:
    """Return candidate indices with score above element-aware threshold.

    Context: Candidate filtering is kept explicit so moment assignment strategies
    (aligned, alternating, random-sign) operate on the same chemically screened set.
    """
    default_thresh = float(config["candidate_score_threshold"])
    thresh_by_elem = config.get("candidate_score_threshold_by_element", {})

    selected: list[int] = []
    for i, features in feature_dict.items():
        if not bool(features.get("is_candidate_element", False)):
            continue

        element = str(features.get("element", ""))
        threshold = float(thresh_by_elem.get(element, default_thresh))
        score = float(features.get("magnetic_score", 0.0))
        if score >= threshold:
            selected.append(int(i))

    selected.sort(key=lambda idx: float(feature_dict[idx].get("magnetic_score", 0.0)), reverse=True)
    return selected
