"""Graph construction and local graph access helpers."""

from __future__ import annotations

from collections import Counter
from typing import Any, Mapping

import networkx as nx
import numpy as np
from ase import Atoms
from ase.neighborlist import neighbor_list

from .bond_graph import build_bond_graph_from_atoms
from .config import (
    expand_pair_cutoffs_for_neighborlist,
    get_default_pair_cutoffs,
    get_present_elements,
    pair_key,
)


def build_graph_with_mic_geometry(
    atoms: Atoms,
    cutoff: float = 1.85,
    pair_cutoffs: Mapping[tuple[str, str], float] | None = None,
    include_self_edges: bool = False,
) -> nx.Graph:
    """Build multi-element bond graph and attach minimum-image convention (MIC) geometry.

    Context: Most downstream features need both connectivity and bond vectors
    under periodic boundary conditions, so this is the standard graph entry point.
    """
    elements_present = get_present_elements(atoms)

    if pair_cutoffs is None:
        pair_cutoffs_local = get_default_pair_cutoffs(elements_present)
        if "C" in elements_present:
            pair_cutoffs_local[pair_key("C", "C")] = float(cutoff)
        elif len(elements_present) == 1:
            elem = elements_present[0]
            pair_cutoffs_local[pair_key(elem, elem)] = float(cutoff)
    else:
        pair_cutoffs_local = {
            pair_key(a, b): float(dist)
            for (a, b), dist in pair_cutoffs.items()
        }

    G = build_bond_graph_from_atoms(
        atoms,
        cutoff=pair_cutoffs_local,
        include_self_edges=include_self_edges,
    )

    if not include_self_edges:
        G.remove_edges_from(nx.selfloop_edges(G))

    attach_mic_geometry_to_graph(
        atoms,
        G,
        cutoff=pair_cutoffs_local,
        include_self_edges=include_self_edges,
    )
    G.graph["pair_cutoffs"] = dict(pair_cutoffs_local)
    G.graph["elements_present"] = list(elements_present)
    return G


def attach_mic_geometry_to_graph(
    atoms: Atoms,
    G: nx.Graph,
    cutoff: float | Mapping[tuple[str, str], float],
    include_self_edges: bool = False,
) -> None:
    """Attach directed minimum-image convention (MIC) bond vectors and lengths.

    Context: Strain, angles, and local-orbital alignment rely on direction-aware
    bond geometry, which must be computed with periodic images handled correctly.
    """
    atoms.wrap()

    nl_cutoff: float | dict[tuple[str, str], float]
    if isinstance(cutoff, Mapping):
        canonical = {pair_key(a, b): float(d) for (a, b), d in cutoff.items()}
        nl_cutoff = expand_pair_cutoffs_for_neighborlist(canonical)
    else:
        nl_cutoff = float(cutoff)

    i_idx, j_idx, shifts = neighbor_list("ijS", atoms, nl_cutoff)

    cell = np.asarray(atoms.cell)
    bond_vectors: dict[tuple[int, int], np.ndarray] = {}
    bond_lengths: dict[tuple[int, int], float] = {}

    for idx in range(len(i_idx)):
        u, v = int(i_idx[idx]), int(j_idx[idx])
        shift = np.asarray(shifts[idx], dtype=float)

        is_boundary_self_bond = u == v and np.any(shift != 0)
        if u > v:
            continue
        if u == v and not (include_self_edges and is_boundary_self_bond):
            continue
        if not G.has_edge(u, v):
            continue

        vec = np.asarray(atoms.positions[v] + shift @ cell - atoms.positions[u], dtype=float)
        length = float(np.linalg.norm(vec))

        if u == v:
            bond_vectors[(u, u)] = vec
            bond_lengths[(u, u)] = length
            continue

        bond_vectors[(u, v)] = vec
        bond_vectors[(v, u)] = -vec
        bond_lengths[(u, v)] = length
        bond_lengths[(v, u)] = length

    G.graph["bond_vectors"] = bond_vectors
    G.graph["bond_lengths"] = bond_lengths


def require_bond_geometry(
    G: nx.Graph,
) -> tuple[dict[tuple[int, int], np.ndarray], dict[tuple[int, int], float]]:
    """Return graph minimum-image convention (MIC) geometry maps or raise a clear error.

    Context: Centralized validation here prevents silent failures when geometry-dependent
    routines are called on a graph that only has connectivity.
    """
    if "bond_vectors" not in G.graph or "bond_lengths" not in G.graph:
        raise ValueError(
            "Graph missing minimum-image convention (MIC) geometry. "
            "Use build_graph_with_mic_geometry(...) "
            "or attach_mic_geometry_to_graph(...)."
        )
    return G.graph["bond_vectors"], G.graph["bond_lengths"]


def get_neighbors(G: nx.Graph, i: int) -> list[int]:
    """Return sorted neighbor indices for node i, excluding self-loops."""
    nbrs = [int(j) for j in G.neighbors(i) if int(j) != int(i)]
    nbrs.sort()
    return nbrs


def count_neighbors_by_element(atoms: Atoms, G: nx.Graph, i: int) -> dict[str, int]:
    """Count neighbors of atom i by element symbol."""
    counts = Counter(str(atoms[j].symbol) for j in get_neighbors(G, i))
    return {elem: int(count) for elem, count in sorted(counts.items())}


def get_bond_vectors(G: nx.Graph, i: int) -> list[np.ndarray]:
    """Return directed minimum-image convention (MIC) bond vectors from i to each neighbor."""
    bond_vectors, _ = require_bond_geometry(G)
    vectors: list[np.ndarray] = []
    for j in get_neighbors(G, i):
        vec = bond_vectors.get((i, j))
        if vec is not None:
            vectors.append(np.asarray(vec, dtype=float))
    return vectors


def get_bond_lengths(G: nx.Graph, i: int) -> list[float]:
    """Return minimum-image convention (MIC) bond lengths from i to each neighbor."""
    _, bond_lengths = require_bond_geometry(G)
    lengths: list[float] = []
    for j in get_neighbors(G, i):
        length = bond_lengths.get((i, j))
        if length is not None:
            lengths.append(float(length))
    return lengths


def get_neighbor_bond_length_map(G: nx.Graph, i: int) -> dict[int, float]:
    """Return map from neighbor index to bond length for atom i."""
    _, bond_lengths = require_bond_geometry(G)
    out: dict[int, float] = {}
    for j in get_neighbors(G, i):
        length = bond_lengths.get((i, j))
        if length is not None:
            out[int(j)] = float(length)
    return out


def get_bond_length_between(G: nx.Graph, i: int, j: int) -> float | None:
    """Return bond length between i and j if present in graph geometry."""
    _, bond_lengths = require_bond_geometry(G)
    length = bond_lengths.get((int(i), int(j)))
    if length is None:
        return None
    return float(length)
