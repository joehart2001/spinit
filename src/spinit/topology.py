"""Compatibility wrapper for :mod:`spinit.graph.topology`."""

from .graph.topology import compute_topological_defect_score, get_atom_ring_sizes, summarize_ring_statistics

__all__ = ["get_atom_ring_sizes", "compute_topological_defect_score", "summarize_ring_statistics"]
