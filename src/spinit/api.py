"""High-level public API for spinit."""

from __future__ import annotations

from typing import Any, Mapping

import networkx as nx
import numpy as np
from ase import Atoms

from .chemistry.features import compute_all_atom_features
from .config import get_present_elements, merge_config
from .graph.bond_graph import (
    e_smallest_local_rings_all_atoms,
    ring_memberships_by_size,
    unique_primitive_rings_by_size,
)
from .graph.graph_ops import build_graph_with_mic_geometry
from .graph.topology import summarize_ring_statistics
from .seeding.assignment import assign_moments_afm_clusters, assign_moments_fm, assign_moments_random_candidates
from .seeding.scoring import assign_all_magnetic_scores


def assign_initial_magnetic_moments(
    atoms: Atoms,
    config: Mapping[str, Any] | None = None,
    strategy: str = "afm_clusters",
    cutoff: float = 1.85,
    pair_cutoffs: Mapping[tuple[str, str], float] | None = None,
    maxlength: int = 12,
    k_hops: int = 5,
    seed: int = 0,
) -> tuple[np.ndarray, dict[int, dict[str, Any]], nx.Graph, dict[str, Any]]:
    """End-to-end multi-element spin-seed pipeline.

    Context: This is the public entry point that connects graph building,
    motif-aware feature detection, scoring, and strategy-based moment assignment.
    """
    cfg = merge_config(config)

    G = build_graph_with_mic_geometry(
        atoms,
        cutoff=cutoff,
        pair_cutoffs=pair_cutoffs,
        include_self_edges=False,
    )

    rings_per_atom, ring_sizes_per_atom = e_smallest_local_rings_all_atoms(
        G,
        maxlength=maxlength,
        k_hops=k_hops,
    )

    feature_dict = compute_all_atom_features(atoms, G, ring_sizes_per_atom, cfg)
    assign_all_magnetic_scores(feature_dict, cfg)

    if strategy == "fm":
        magmoms = assign_moments_fm(atoms, feature_dict, cfg)
    elif strategy == "afm_clusters":
        magmoms = assign_moments_afm_clusters(atoms, G, feature_dict, cfg)
    elif strategy == "random_candidates":
        magmoms = assign_moments_random_candidates(atoms, feature_dict, cfg, seed=seed)
    else:
        valid = ["fm", "afm_clusters", "random_candidates"]
        raise ValueError(f"Unknown strategy '{strategy}'. Expected one of {valid}.")

    ring_info = {
        "rings_per_atom": rings_per_atom,
        "ring_sizes_per_atom": ring_sizes_per_atom,
        "membership_histogram": ring_memberships_by_size(ring_sizes_per_atom),
        "unique_ring_histogram": unique_primitive_rings_by_size(rings_per_atom, G),
        "summary": summarize_ring_statistics(ring_sizes_per_atom),
        "elements_present": get_present_elements(atoms),
        "pair_cutoffs": dict(G.graph.get("pair_cutoffs", {})),
    }

    return magmoms, feature_dict, G, ring_info
