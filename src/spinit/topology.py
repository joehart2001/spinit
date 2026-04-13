"""Ring and topology feature helpers."""

from __future__ import annotations

from collections import Counter
from typing import Any, Mapping, Sequence

import numpy as np


def get_atom_ring_sizes(ring_sizes_per_atom: Sequence[Sequence[int]], i: int) -> list[int]:
    """Return local ring sizes associated with atom i."""
    if i < 0 or i >= len(ring_sizes_per_atom):
        return []
    return list(map(int, ring_sizes_per_atom[i]))


def compute_topological_defect_score(
    ring_sizes_for_atom: Sequence[int],
    element: str,
    config: Mapping[str, Any],
) -> float:
    """Return topology score, active mainly for configured carbon-like elements."""
    if str(element) not in set(map(str, config["elements"]["topology_active_elements"])):
        return 0.0

    weights = config["topology_weights"]
    score = 0.0
    for ring_size in ring_sizes_for_atom:
        if ring_size == 6:
            score += float(weights["size_6"])
        elif ring_size in (5, 7):
            score += float(weights["size_5_7"])
        else:
            score += float(weights["other_non_hex"])
    return float(score)


def summarize_ring_statistics(ring_sizes_per_atom: Sequence[Sequence[int]]) -> dict[str, Any]:
    """Return compact structure-level ring-membership summary."""
    all_sizes = [int(size) for atom_sizes in ring_sizes_per_atom for size in atom_sizes]
    size_counts = dict(Counter(all_sizes))

    memberships_per_atom = [len(atom_sizes) for atom_sizes in ring_sizes_per_atom]
    mean_memberships = float(np.mean(memberships_per_atom)) if memberships_per_atom else 0.0
    max_memberships = int(np.max(memberships_per_atom)) if memberships_per_atom else 0

    non_hex_memberships = sum(1 for size in all_sizes if size != 6)
    total_memberships = len(all_sizes)
    non_hex_fraction = (float(non_hex_memberships) / float(total_memberships)) if total_memberships > 0 else 0.0

    return {
        "size_counts": size_counts,
        "mean_memberships_per_atom": mean_memberships,
        "max_memberships_per_atom": max_memberships,
        "total_memberships": int(total_memberships),
        "non_hex_membership_fraction": non_hex_fraction,
    }
