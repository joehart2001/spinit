"""Compatibility wrapper for legacy magnetic_moments imports.

This module re-exports the new modular implementation from the spinit package.
Prefer importing from `spinit` or `spinit.api` directly in new code.
"""

from .api import assign_initial_magnetic_moments
from .config import default_config, get_default_pair_cutoffs, get_present_elements, merge_config
from .reporting import make_feature_table, print_candidate_summary, print_feature_summary

__all__ = [
    "assign_initial_magnetic_moments",
    "default_config",
    "merge_config",
    "get_present_elements",
    "get_default_pair_cutoffs",
    "make_feature_table",
    "print_feature_summary",
    "print_candidate_summary",
]
