"""Compatibility wrapper for :mod:`spinit.graph.graph_ops`."""

from .graph.graph_ops import (
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

__all__ = [
    "build_graph_with_mic_geometry",
    "attach_mic_geometry_to_graph",
    "require_bond_geometry",
    "get_neighbors",
    "count_neighbors_by_element",
    "get_bond_vectors",
    "get_bond_lengths",
    "get_neighbor_bond_length_map",
    "get_bond_length_between",
]
