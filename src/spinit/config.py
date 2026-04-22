"""Configuration defaults and helpers for spinit."""

from __future__ import annotations

from copy import deepcopy
from itertools import combinations_with_replacement
from typing import Any, Mapping, Sequence

from ase import Atoms
from ase.data import atomic_numbers, covalent_radii


def pair_key(a: str, b: str) -> tuple[str, str]:
    """Return normalized element-pair key."""
    return tuple(sorted((str(a), str(b))))


def get_present_elements(atoms: Atoms) -> list[str]:
    """Return sorted list of unique element symbols present in the structure."""
    return sorted(set(map(str, atoms.get_chemical_symbols())))


def get_default_pair_cutoffs(elements_present: Sequence[str]) -> dict[tuple[str, str], float]:
    """Return default element-pair cutoffs for the present elements."""
    base: dict[tuple[str, str], float] = {
        ("C", "C"): 1.85,
        ("C", "H"): 1.25,
        ("C", "N"): 1.70,
        ("C", "O"): 1.65,
        ("H", "H"): 1.00,
        ("H", "N"): 1.25,
        ("H", "O"): 1.20,
        ("N", "N"): 1.75,
        ("N", "O"): 1.75,
        ("O", "O"): 1.70,
    }

    fallback_scale = 1.20
    cutoffs: dict[tuple[str, str], float] = {}

    for a, b in combinations_with_replacement(sorted(elements_present), 2):
        key = pair_key(a, b)
        if key in base:
            cutoffs[key] = float(base[key])
            continue

        try:
            ra = float(covalent_radii[atomic_numbers[a]])
            rb = float(covalent_radii[atomic_numbers[b]])
            cutoffs[key] = fallback_scale * (ra + rb)
        except Exception:
            cutoffs[key] = 2.0

    return cutoffs


def expand_pair_cutoffs_for_neighborlist(
    pair_cutoffs: Mapping[tuple[str, str], float],
) -> dict[tuple[str, str], float]:
    """Expand canonical pair cutoffs into both pair orders for ASE neighbor_list."""
    expanded: dict[tuple[str, str], float] = {}
    for (a, b), dist in pair_cutoffs.items():
        expanded[(a, b)] = float(dist)
        expanded[(b, a)] = float(dist)
    return expanded


def default_config() -> dict[str, Any]:
    """Return default thresholds and weights for heuristic feature detection."""
    return {
        "carbon_symbol": "C",
        "elements": {
            "allow_n_weak_candidate": False,
            "topology_active_elements": ["C"],
            "heteroatom_modifier_elements": ["H", "O"],
            "usually_zero_moment_elements": ["H", "O"],
        },
        "hybridization": {
            "sp_linear_target_deg": 180.0,
            "sp2_target_deg": 120.0,
            "sp3_target_deg": 109.5,
            "sp_linear_tol_deg": 20.0,
            "sp2_tol_deg": 18.0,
            "sp3_tol_deg": 14.0,
            "sp2_mad_max_deg": 16.0,
            "sp3_mad_max_deg": 12.0,
        },
        "pi_unsaturated": {
            "angle_mad_min_deg": 10.0,
        },
        "sp_site": {
            "linearity_tolerance_deg": 20.0,
            "allow_unmodified_linear_sp_as_radical_like": False,
            "require_radical_like_for_intrinsic": True,
            "satisfied_sp_intrinsic_scale": 0.0,
        },
        "motif_detection": {
            "network_min_component_size": 20,
            "network_requires_carbon": True,
            "molecular_fragment_max_size": 16,
            "small_hydrocarbon_max_size": 12,
            "carbonyl_co_max": 1.30,
            "carboxyl_co_single_max": 1.46,
            "peroxide_oo_max": 1.60,
            "o2_oo_max": 1.35,
            "nitrile_cn_max": 1.22,
            "n2_nn_max": 1.20,
            "h2_hh_max": 0.95,
        },
        "sp2_strain": {
            "pyramidalization_min_deg": 8.0,
            "pi_alignment_min": 0.90,
        },
        "topology_weights": {
            "size_6": 0.0,
            "size_5_7": 1.0,
            "other_non_hex": 0.7,
        },
        "intrinsic_score_weights": {
            "C": {
                "undercoordination": 2.5,
                "sigma_dangling": 3.0,
                "pi_unsaturated": 1.6,
                "strained_sp2": 1.3,
                "pi_misalignment": 2.0,
                "topological_defect": 0.35,
            },
            "N": {
                "undercoordination": 1.2,
                "sigma_dangling": 1.4,
                "pi_unsaturated": 0.8,
                "strained_sp2": 0.6,
                "pi_misalignment": 0.8,
                "topological_defect": 0.15,
            },
            "O": {
                "undercoordination": 0.8,
                "sigma_dangling": 1.0,
                "pi_unsaturated": 0.2,
                "strained_sp2": 0.0,
                "pi_misalignment": 0.0,
                "topological_defect": 0.0,
            },
            "H": {
                "undercoordination": 0.4,
                "sigma_dangling": 0.4,
                "pi_unsaturated": 0.0,
                "strained_sp2": 0.0,
                "pi_misalignment": 0.0,
                "topological_defect": 0.0,
            },
        },
        "score_controls": {
            "pi_misalignment_scale_when_strained_sp2": 0.5,
            "topology_requires_plausible_site": True,
            "topology_scale_for_plausible_site": 0.5,
            "neighbor_modifier_min": -4.0,
            "neighbor_modifier_max": 2.0,
        },
        "motif_scoring": {
            "closed_shell_penalty": -1.4,
            "molecular_closed_shell_penalty": -0.6,
            "radical_like_bonus": 2.3,
            "heteroatom_host_bonus": 0.7,
            "label_weights": {
                "graphitic_N": 0.35,
                "pyridinic_like_N": 0.55,
                "pyrrolic_like_N": 0.40,
                "undercoordinated_oxygen": 1.0,
                "undercoordinated_nitrogen": 1.0,
                "isolated_hydrogen": 1.2,
                "adsorbed_hydrogen": 0.8,
                "molecular_oxygen_like_O": 1.0,
                "molecular_nitrogen_like_N": -0.7,
                "water_like_O": -1.0,
                "hydroxyl_like_O": -0.6,
                "ether_like_O": -0.5,
                "epoxy_like_O": -0.4,
                "carbonyl_like_O": -0.6,
                "carboxyl_like_O": -0.7,
                "peroxide_like_O": -0.3,
                "co_like_C": -0.8,
                "co2_like_C": -1.0,
                "small_hydrocarbon_fragment_C": -0.6,
            },
        },
        "neighbor_modifiers": {
            "C": {"H": -1.0, "O": -1.2, "N": -0.2},
            "N": {"H": -0.5, "O": -0.6},
            "O": {"H": -0.3},
        },
        "motif_neighbor_modifiers": {
            "C": {
                "hydroxyl_like_O": -0.5,
                "ether_like_O": -0.35,
                "epoxy_like_O": -0.35,
                "carbonyl_like_O": -0.45,
                "carboxyl_like_O": -0.55,
                "water_like_O": -0.2,
                "undercoordinated_oxygen": 0.8,
                "graphitic_N": 0.3,
                "pyridinic_like_N": 0.5,
                "pyrrolic_like_N": 0.3,
                "undercoordinated_nitrogen": 0.9,
                "adsorbed_hydrogen": 0.4,
                "c_dangling_terminator_H": -0.6,
            },
            "N": {
                "hydroxyl_like_O": -0.3,
                "undercoordinated_oxygen": 0.4,
            },
        },
        "candidate_score_threshold": 2.8,
        "candidate_score_threshold_by_element": {
            "C": 2.8,
            "N": 3.4,
            "O": 3.6,
            "H": 4.0,
        },
        "moment_assignment": {
            "strong_moment": 1.0,
            "moderate_moment": 0.6,
            "weak_moment": 0.35,
            "min_candidate_moment": 0.25,
            "max_candidate_moment": 1.0,
            "score_to_moment_scale": 0.15,
            "soft_moment_scale": 0.35,
            "soft_min_candidate_moment": 0.05,
            "soft_max_candidate_moment": 0.50,
            "o2_parallel_override_in_afm": True,
            "o2_parallel_sign": 1.0,
            "element_moment_scale": {
                "C": 1.0,
                "N": 0.7,
                "O": 0.6,
                "H": 0.4,
            },
        },
    }


def merge_nested_dict(base: Mapping[str, Any], override: Mapping[str, Any]) -> dict[str, Any]:
    """Recursively merge dictionaries without mutating inputs."""
    merged: dict[str, Any] = deepcopy(dict(base))
    for key, value in override.items():
        if isinstance(value, Mapping) and isinstance(merged.get(key), Mapping):
            merged[key] = merge_nested_dict(dict(merged[key]), value)
        else:
            merged[key] = value
    return merged


def merge_config(config: Mapping[str, Any] | None) -> dict[str, Any]:
    """Return merged config (defaults overridden by user values if provided)."""
    if config is None:
        return default_config()
    return merge_nested_dict(default_config(), config)
