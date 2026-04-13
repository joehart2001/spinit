from ase import Atoms

from spinit import assign_initial_magnetic_moments


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
