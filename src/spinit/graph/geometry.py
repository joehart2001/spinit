"""Local geometry utilities from bond vectors."""

from __future__ import annotations

from itertools import combinations
from typing import Sequence

import networkx as nx
import numpy as np

from .graph_ops import get_bond_vectors


def unit_vector(v: Sequence[float]) -> np.ndarray:
    """Return normalized vector; zeros for near-zero input."""
    arr = np.asarray(v, dtype=float)
    norm = float(np.linalg.norm(arr))
    if norm < 1e-12:
        return np.zeros_like(arr)
    return arr / norm


def angle_between(v1: Sequence[float], v2: Sequence[float], degrees: bool = True) -> float:
    """Return angle between vectors."""
    u1 = unit_vector(v1)
    u2 = unit_vector(v2)

    n1 = float(np.linalg.norm(u1))
    n2 = float(np.linalg.norm(u2))
    if n1 < 1e-12 or n2 < 1e-12:
        return float("nan")

    cos_theta = float(np.clip(np.dot(u1, u2), -1.0, 1.0))
    angle = float(np.arccos(cos_theta))
    if degrees:
        return float(np.degrees(angle))
    return angle


def pairwise_angles(vectors: Sequence[Sequence[float]]) -> list[float]:
    """Return all pairwise angles between vectors."""
    angles: list[float] = []
    for v1, v2 in combinations(vectors, 2):
        angle = angle_between(v1, v2, degrees=True)
        if not np.isnan(angle):
            angles.append(float(angle))
    return angles


def get_bond_angles(G: nx.Graph, i: int) -> list[float]:
    """Return all bond angles around node i."""
    return pairwise_angles(get_bond_vectors(G, i))


def mean_abs_deviation(values: Sequence[float], target: float) -> float | None:
    """Return mean absolute deviation from target; None for empty input."""
    if len(values) == 0:
        return None
    arr = np.asarray(values, dtype=float)
    return float(np.mean(np.abs(arr - float(target))))
