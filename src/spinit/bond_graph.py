"""Bond-graph construction and local ring statistics utilities."""

from __future__ import annotations

from collections import Counter
from itertools import combinations
from typing import Any, Mapping

import networkx as nx
import numpy as np
from ase.neighborlist import neighbor_list


def _prepare_neighborlist_cutoff(cutoff: float | Mapping[tuple[Any, Any], float]) -> float | dict[tuple[Any, Any], float]:
    """Return cutoff object in a form suitable for ASE neighbor_list."""
    if not isinstance(cutoff, Mapping):
        return float(cutoff)

    pair_cutoffs: dict[tuple[Any, Any], float] = {}
    for pair, value in cutoff.items():
        if not isinstance(pair, tuple) or len(pair) != 2:
            raise ValueError("Pair cutoff keys must be 2-tuples, e.g. ('C', 'H').")
        a, b = pair
        dist = float(value)
        pair_cutoffs[(a, b)] = dist
        pair_cutoffs[(b, a)] = dist
    return pair_cutoffs


def build_bond_graph_from_atoms(
    atoms,
    cutoff: float | Mapping[tuple[Any, Any], float] = 1.85,
    include_self_edges: bool = True,
) -> nx.Graph:
    """Build an undirected bond graph from atoms."""
    atoms.wrap()

    nl_cutoff = _prepare_neighborlist_cutoff(cutoff)
    i, j, shifts = neighbor_list("ijS", atoms, nl_cutoff)

    G = nx.Graph()
    G.add_nodes_from(range(len(atoms)))

    for idx in range(len(i)):
        u, v = int(i[idx]), int(j[idx])
        is_boundary_self_bond = u == v and np.any(shifts[idx] != 0)
        if u < v or (include_self_edges and is_boundary_self_bond):
            G.add_edge(u, v)

    return G


def cycle_key(cycle_nodes):
    """Return a stable key for a cycle, invariant to rotation and traversal direction.

    Context: The same ring can be discovered from different start atoms and path
    directions, so this key lets us deduplicate equivalent cycle descriptions.
    """
    if len(cycle_nodes) >= 2 and cycle_nodes[0] == cycle_nodes[-1]:
        cycle_nodes = cycle_nodes[:-1]

    cycle = list(map(int, cycle_nodes))
    if not cycle:
        return tuple()

    min_node = min(cycle)
    min_indices = [i for i, x in enumerate(cycle) if x == min_node]
    candidates = []
    for idx in min_indices:
        forward = cycle[idx:] + cycle[:idx]
        candidates.append(tuple(forward))

        reverse_cycle = cycle[::-1]
        reverse_idx = reverse_cycle.index(min_node)
        backward = reverse_cycle[reverse_idx:] + reverse_cycle[:reverse_idx]
        candidates.append(tuple(backward))

    return min(candidates)


def e_smallest_local_rings_for_v(
    G: nx.Graph,
    v: int,
    maxlength: int = 12,
    k_hops: int | None = 5,
) -> list[list[int]]:
    """For node v with degree e, find up to e smallest local rings touching v.

    Context: Local ring signatures are used as supporting topology features, so
    we keep the search near each atom instead of enumerating all global cycles.
    """
    neighbors = list(G.neighbors(v))
    degree = len(neighbors)
    if degree < 2:
        return []

    if k_hops is not None:
        neighborhood = set(nx.single_source_shortest_path_length(G, v, cutoff=k_hops).keys())
        H = G.subgraph(neighborhood).copy()
    else:
        H = G.copy()

    if not H.has_node(v):
        return []
    H.remove_node(v)

    cycles = []
    seen = set()

    for u, w in combinations(neighbors, 2):
        if not H.has_node(u) or not H.has_node(w):
            continue
        try:
            path = nx.shortest_path(H, source=u, target=w)
        except nx.NetworkXNoPath:
            continue

        cyc = [v] + path + [v]
        ring_size = len(cyc) - 1
        if ring_size > maxlength:
            continue

        key = cycle_key(cyc)
        if key in seen:
            continue
        seen.add(key)
        cycles.append(cyc)

    cycles.sort(key=lambda c: (len(c) - 1, c))
    return cycles[:degree]


def e_smallest_local_rings_all_atoms(
    G: nx.Graph,
    maxlength: int = 12,
    k_hops: int | None = 5,
) -> tuple[list[list[list[int]]], list[list[int]]]:
    """Return per-atom local rings and ring sizes.

    Context: This runs the local-ring detection once for the full structure so
    per-atom topology metrics can be reused in feature extraction and reporting.
    """
    num_nodes = G.number_of_nodes()
    rings_per_atom = []
    ring_sizes_per_atom = []

    for v in range(num_nodes):
        rings_v = e_smallest_local_rings_for_v(G, v, maxlength=maxlength, k_hops=k_hops)
        rings_per_atom.append(rings_v)
        ring_sizes_per_atom.append([len(cyc) - 1 for cyc in rings_v])

    return rings_per_atom, ring_sizes_per_atom


def ring_memberships_by_size(
    ring_sizes_per_atom: list[list[int]],
    sizes=(3, 4, 5, 6, 7, 8),
) -> dict[int, int]:
    """Count per-atom ring memberships by ring size (not unique ring count)."""
    counts = Counter()
    allowed = set(sizes)
    for atom_ring_sizes in ring_sizes_per_atom:
        for ring_size in atom_ring_sizes:
            if ring_size in allowed:
                counts[ring_size] += 1
    return {size: int(counts.get(size, 0)) for size in sizes}


def is_primitive_ring(G: nx.Graph, cycle_nodes: list[int]) -> bool:
    """Check if ring has no chords (internal edges)."""
    if cycle_nodes[0] == cycle_nodes[-1]:
        cycle_nodes = cycle_nodes[:-1]

    for i, u in enumerate(cycle_nodes):
        for j in range(i + 2, len(cycle_nodes)):
            v = cycle_nodes[j]
            if j == len(cycle_nodes) - 1 and i == 0:
                continue
            if G.has_edge(u, v):
                return False
    return True


def unique_primitive_rings_by_size(
    rings_per_atom: list[list[list[int]]],
    G: nx.Graph,
    sizes=(3, 4, 5, 6, 7, 8),
) -> dict[int, int]:
    """Deduplicate cycles across all atoms to estimate unique primitive ring counts.

    Context: Per-atom ring membership overcounts shared rings, so this provides a
    structure-level primitive-ring summary for diagnostics and sanity checks.
    """
    allowed = set(sizes)
    seen = set()
    counts = Counter()

    for rings_v in rings_per_atom:
        for cyc in rings_v:
            ring_size = len(cyc) - 1
            if ring_size not in allowed:
                continue
            if not is_primitive_ring(G, cyc):
                continue
            key = cycle_key(cyc)
            if key in seen:
                continue
            seen.add(key)
            counts[ring_size] += 1

    return {size: int(counts.get(size, 0)) for size in sizes}
