"""Compatibility wrapper for :mod:`spinit.seeding.assignment`."""

from .seeding.assignment import (
    assign_moments_afm_clusters,
    assign_moments_fm,
    assign_moments_random_candidates,
    candidate_moment_amplitude,
)

__all__ = [
    "candidate_moment_amplitude",
    "assign_moments_fm",
    "assign_moments_afm_clusters",
    "assign_moments_random_candidates",
]
