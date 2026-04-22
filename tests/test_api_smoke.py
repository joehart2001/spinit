from __future__ import annotations

import math

from ase import Atoms

from spinit.api import assign_initial_magnetic_moments


def _active_sites(magmoms) -> set[int]:
    return {i for i, m in enumerate(magmoms) if abs(float(m)) > 1e-12}


def test_assign_initial_magnetic_moments_smoke() -> None:
    atoms = Atoms("C3", positions=[[0.0, 0.0, 0.0], [1.28, 0.0, 0.0], [2.56, 0.0, 0.0]], pbc=False)

    magmoms, feature_dict, graph, ring_info = assign_initial_magnetic_moments(
        atoms,
        strategy="fm",
        cutoff=1.85,
        maxlength=8,
        k_hops=3,
        seed=0,
    )

    assert len(magmoms) == len(atoms)
    assert len(feature_dict) == len(atoms)
    assert graph.number_of_nodes() == len(atoms)
    assert "elements_present" in ring_info


def test_o2_afm_override_parallel_triplet_like() -> None:
    atoms = Atoms("O2", positions=[[0.0, 0.0, 0.0], [1.23, 0.0, 0.0]], pbc=False)

    magmoms, feature_dict, graph, ring_info = assign_initial_magnetic_moments(
        atoms,
        strategy="afm_clusters",
        cutoff=1.85,
        maxlength=8,
        k_hops=3,
        seed=0,
    )

    assert len(magmoms) == 2
    assert abs(float(magmoms[0])) > 1e-12
    assert abs(float(magmoms[1])) > 1e-12
    assert float(magmoms[0]) * float(magmoms[1]) > 0.0
    assert feature_dict[0]["motif_label"] == "molecular_oxygen_like_O"
    assert feature_dict[1]["motif_label"] == "molecular_oxygen_like_O"


def test_fm_soft_reduces_fm_amplitudes_on_same_sites() -> None:
    atoms = Atoms("C3", positions=[[0.0, 0.0, 0.0], [1.28, 0.0, 0.0], [2.56, 0.0, 0.0]], pbc=False)

    hard_magmoms, _, _, _ = assign_initial_magnetic_moments(
        atoms,
        strategy="fm",
        cutoff=1.85,
        maxlength=8,
        k_hops=3,
        seed=0,
    )
    soft_magmoms, _, _, _ = assign_initial_magnetic_moments(
        atoms,
        strategy="fm_soft",
        cutoff=1.85,
        maxlength=8,
        k_hops=3,
        seed=0,
    )

    active = _active_sites(hard_magmoms)
    assert active == _active_sites(soft_magmoms)
    assert len(active) > 0
    for i in active:
        assert 0.0 < abs(float(soft_magmoms[i])) < abs(float(hard_magmoms[i]))


def test_afm_soft_reduces_amplitudes_and_keeps_o2_parallel_override() -> None:
    atoms = Atoms("O2", positions=[[0.0, 0.0, 0.0], [1.23, 0.0, 0.0]], pbc=False)

    hard_magmoms, _, _, _ = assign_initial_magnetic_moments(
        atoms,
        strategy="afm_clusters",
        cutoff=1.85,
        maxlength=8,
        k_hops=3,
        seed=0,
    )
    soft_magmoms, _, _, _ = assign_initial_magnetic_moments(
        atoms,
        strategy="afm_soft",
        cutoff=1.85,
        maxlength=8,
        k_hops=3,
        seed=0,
    )

    assert _active_sites(hard_magmoms) == {0, 1}
    assert _active_sites(soft_magmoms) == {0, 1}
    assert float(soft_magmoms[0]) * float(soft_magmoms[1]) > 0.0
    for i in (0, 1):
        assert abs(float(soft_magmoms[i])) < abs(float(hard_magmoms[i]))


def test_balanced_random_soft_gives_small_near_zero_seed() -> None:
    positions = [
        [1.40 * math.cos(2.0 * math.pi * k / 6.0), 1.40 * math.sin(2.0 * math.pi * k / 6.0), 0.0]
        for k in range(6)
    ]
    atoms = Atoms("C6", positions=positions, pbc=False)

    magmoms, _, _, _ = assign_initial_magnetic_moments(
        atoms,
        strategy="balanced_random_soft",
        cutoff=1.85,
        maxlength=12,
        k_hops=5,
        seed=0,
    )

    active = _active_sites(magmoms)
    assert active == set(range(6))
    assert any(float(magmoms[i]) > 0.0 for i in active)
    assert any(float(magmoms[i]) < 0.0 for i in active)
    assert abs(float(sum(magmoms))) < 1e-12
