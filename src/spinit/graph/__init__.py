"""Graph and geometry utilities for connectivity/topology analysis."""

from .bond_graph import (
    build_bond_graph_from_atoms,
    cycle_key,
    e_smallest_local_rings_all_atoms,
    e_smallest_local_rings_for_v,
    is_primitive_ring,
    ring_memberships_by_size,
    unique_primitive_rings_by_size,
)
from .geometry import angle_between, get_bond_angles, mean_abs_deviation, pairwise_angles, unit_vector
from .graph_ops import (
    attach_mic_geometry_to_graph,
    build_graph_with_mic_geometry,
    count_neighbors_by_element,
    get_bond_length_between,
    get_bond_lengths,
    get_bond_vectors,
    get_neighbor_bond_length_map,
    get_neighbors,
    require_bond_geometry,
)
from .topology import compute_topological_defect_score, get_atom_ring_sizes, summarize_ring_statistics

__all__ = [
    "build_bond_graph_from_atoms",
    "cycle_key",
    "e_smallest_local_rings_for_v",
    "e_smallest_local_rings_all_atoms",
    "ring_memberships_by_size",
    "is_primitive_ring",
    "unique_primitive_rings_by_size",
    "build_graph_with_mic_geometry",
    "attach_mic_geometry_to_graph",
    "require_bond_geometry",
    "get_neighbors",
    "count_neighbors_by_element",
    "get_bond_vectors",
    "get_bond_lengths",
    "get_neighbor_bond_length_map",
    "get_bond_length_between",
    "unit_vector",
    "angle_between",
    "pairwise_angles",
    "get_bond_angles",
    "mean_abs_deviation",
    "get_atom_ring_sizes",
    "compute_topological_defect_score",
    "summarize_ring_statistics",
]
