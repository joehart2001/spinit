"""Compatibility wrapper for :mod:`spinit.chemistry.features`."""

from .chemistry.features import compute_all_atom_features, compute_atom_features, compute_element_specific_features

__all__ = ["compute_element_specific_features", "compute_atom_features", "compute_all_atom_features"]
