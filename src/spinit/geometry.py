"""Compatibility wrapper for :mod:`spinit.graph.geometry`."""

from .graph.geometry import angle_between, get_bond_angles, mean_abs_deviation, pairwise_angles, unit_vector

__all__ = ["unit_vector", "angle_between", "pairwise_angles", "get_bond_angles", "mean_abs_deviation"]
