"""Magnetic moment assignment strategies."""

from __future__ import annotations

from collections import deque
from typing import Any, Mapping

import networkx as nx
import numpy as np
from ase import Atoms

from .scoring import select_magnetic_candidate_sites


def candidate_moment_amplitude(features: Mapping[str, Any], config: Mapping[str, Any]) -> float:
    """Map candidate score to moment amplitude with element-aware scaling."""
    mcfg = config["moment_assignment"]
    element = str(features.get("element", ""))

    score = max(0.0, float(features.get("magnetic_score", 0.0)))
    amp = float(mcfg["min_candidate_moment"]) + float(mcfg["score_to_moment_scale"]) * score

    if features.get("sigma_dangling_bond", False) or int(features.get("coordination", 9)) <= 1:
        amp = max(amp, float(mcfg["strong_moment"]))
    elif features.get("is_strained_sp2", False) or features.get("pi_unsaturated", False):
        amp = max(amp, float(mcfg["moderate_moment"]))
    else:
        amp = max(amp, float(mcfg["weak_moment"]))

    scale_map = mcfg.get("element_moment_scale", {})
    amp *= float(scale_map.get(element, 1.0))

    return float(np.clip(amp, float(mcfg["min_candidate_moment"]), float(mcfg["max_candidate_moment"])))


def assign_moments_fm(
    atoms: Atoms,
    feature_dict: Mapping[int, Mapping[str, Any]],
    config: Mapping[str, Any],
) -> np.ndarray:
    """Assign ferromagnetic (+) moments on selected candidates."""
    magmoms = np.zeros(len(atoms), dtype=float)
    for i in select_magnetic_candidate_sites(feature_dict, config):
        magmoms[i] = candidate_moment_amplitude(feature_dict[i], config)
    return magmoms


def assign_moments_afm_clusters(
    atoms: Atoms,
    G: nx.Graph,
    feature_dict: Mapping[int, Mapping[str, Any]],
    config: Mapping[str, Any],
) -> np.ndarray:
    """Assign AFM-like signs by alternating signs inside connected candidate clusters."""
    magmoms = np.zeros(len(atoms), dtype=float)
    candidate_sites = select_magnetic_candidate_sites(feature_dict, config)
    if len(candidate_sites) == 0:
        return magmoms

    subgraph = G.subgraph(list(candidate_sites)).copy()
    clusters = [sorted([int(i) for i in comp]) for comp in nx.connected_components(subgraph)]
    clusters.sort(key=lambda cluster: (len(cluster), cluster), reverse=True)

    for cluster in clusters:
        sign_labels: dict[int, int] = {int(cluster[0]): 1}
        queue: deque[int] = deque([int(cluster[0])])
        while queue:
            u = int(queue.popleft())
            for v_raw in subgraph.neighbors(u):
                v = int(v_raw)
                if v not in sign_labels:
                    sign_labels[v] = -sign_labels[u]
                    queue.append(v)
                elif sign_labels[v] == sign_labels[u]:
                    continue

        for i in cluster:
            sign = float(sign_labels.get(i, 1))
            amp = candidate_moment_amplitude(feature_dict[i], config)
            magmoms[i] = sign * amp

    return magmoms


def assign_moments_random_candidates(
    atoms: Atoms,
    feature_dict: Mapping[int, Mapping[str, Any]],
    config: Mapping[str, Any],
    seed: int = 0,
) -> np.ndarray:
    """Assign random +/- signs to candidate-site moments with reproducible seed."""
    rng = np.random.default_rng(seed)
    magmoms = np.zeros(len(atoms), dtype=float)

    for i in select_magnetic_candidate_sites(feature_dict, config):
        sign = float(rng.choice([-1.0, 1.0]))
        amp = candidate_moment_amplitude(feature_dict[i], config)
        magmoms[i] = sign * amp

    return magmoms
