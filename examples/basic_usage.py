"""Minimal usage example for spinit."""

from ase.io import read

from spinit.api import assign_initial_magnetic_moments


atoms = read("structure.xyz")

magmoms, feature_dict, graph, ring_info = assign_initial_magnetic_moments(
    atoms,
    strategy="afm_clusters",
    cutoff=1.85,
    maxlength=12,
    k_hops=5,
    seed=0,
)

atoms.set_initial_magnetic_moments(magmoms)

print("Selected non-zero moments:", sum(abs(float(m)) > 1e-12 for m in magmoms))
print("Elements present:", ring_info["elements_present"])
print("Ring summary:", ring_info["summary"])
