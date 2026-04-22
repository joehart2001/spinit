"""Magnetic scoring and moment assignment strategies."""

from .assignment import (
    assign_moments_afm_soft,
    assign_moments_afm_clusters,
    assign_moments_balanced_random_soft,
    assign_moments_fm,
    assign_moments_fm_soft,
    assign_moments_random_candidates,
    candidate_moment_amplitude,
    soft_candidate_moment_amplitude,
)
from .scoring import (
    assign_all_magnetic_scores,
    compute_intrinsic_site_score,
    compute_magnetic_score,
    compute_motif_intrinsic_modifier,
    compute_neighbor_modifier,
    get_weights_for_element,
    select_magnetic_candidate_sites,
)

__all__ = [
    "get_weights_for_element",
    "compute_motif_intrinsic_modifier",
    "compute_neighbor_modifier",
    "compute_intrinsic_site_score",
    "compute_magnetic_score",
    "assign_all_magnetic_scores",
    "select_magnetic_candidate_sites",
    "candidate_moment_amplitude",
    "soft_candidate_moment_amplitude",
    "assign_moments_fm",
    "assign_moments_fm_soft",
    "assign_moments_afm_clusters",
    "assign_moments_afm_soft",
    "assign_moments_random_candidates",
    "assign_moments_balanced_random_soft",
]
