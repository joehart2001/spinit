"""Compatibility wrapper for :mod:`spinit.graph.bond_graph`."""

from .graph.bond_graph import (
    build_bond_graph_from_atoms,
    cycle_key,
    e_smallest_local_rings_all_atoms,
    e_smallest_local_rings_for_v,
    is_primitive_ring,
    ring_memberships_by_size,
    unique_primitive_rings_by_size,
)

__all__ = [
    "build_bond_graph_from_atoms",
    "cycle_key",
    "e_smallest_local_rings_for_v",
    "e_smallest_local_rings_all_atoms",
    "ring_memberships_by_size",
    "is_primitive_ring",
    "unique_primitive_rings_by_size",
]
