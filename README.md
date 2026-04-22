![spinit logo](spinit_logo.svg)


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19595734.svg)](https://doi.org/10.5281/zenodo.19595734) [![PyPI](https://img.shields.io/pypi/v/spinit.svg)](https://pypi.org/project/spinit/)

## Physically informed spin initialisation for atomistic structures


`spinit` seeds **initial magnetic moments** to `ase.Atoms` using local chemistry, topology, and geometry criteria to provide physically motivated spin initialisation. This is useful for e.g. spin-polarised DFT of large structures where spins cannot (reasonably) be manually assigned.

It is designed for carbon-rich systems, including graphenic, nanoporous, and amorphous carbon, as well as mixed chemistries containing heteroatoms such as H, O, and N.

## Scientific Scope

`spinit` is a **heuristic spin seed generator**.

- It does not prove magnetic ground states.
- It provides physically motivated initial spin patterns to improve exploration of plausible spin solutions.
- It combines undercoordination, local geometry/strain, ring topology, and explicit motif detection.
- Motif labels are heuristic local-environment descriptors (`*_like`) and do not by themselves assign formal bond order, oxidation state, or molecular charge.

## Features

- PBC-aware bond graph construction (ASE + MIC vectors).
- Ring/topology analysis from local graph cycles.
- Hybridization-like local classification (`sp`, `sp2`, `sp3`, distorted variants).
- Motif-aware handling of O/N/H/C-centered local environments.
- Connectivity-aware classification of network-bound vs molecular fragments.
- Transparent additive scoring with configurable thresholds and weights.
- Multiple assignment strategies:
  - `fm`
  - `fm_soft`
  - `afm_clusters`
  - `afm_soft`
  - `random_candidates`
  - `balanced_random_soft`

## Installation

```bash
git clone https://github.com/joehart2001/spinit.git
cd spinit
pip install -e .
```

Optional reporting extras:

```bash
pip install -e .[reporting]
```

## Quickstart

```python
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
```

## SCF Restarts And Convergence Rescue

If a spin-polarised DFT calculation is struggling to converge, reinitialisaiton
of the spin with a smaller or less prescriptive magnetic seed can sometimes help 
with SCF convergence.

The restart-oriented strategies in `spinit` keep the same chemically screened
candidate sites, but soften the assigned initial moments:

- `fm_soft`: small all-positive moments on the selected candidate atoms. This is
  a good first retry when a stronger AFM/FM seed is too aggressive.
- `afm_soft`: the same cluster-based sign pattern as `afm_clusters`, but with
  reduced amplitudes. Use this when you still want local AFM-like structure
  without a large initial spin splitting.
- `balanced_random_soft`: small random `+/-` moments chosen to keep the total
  initial magnetisation close to zero. This is useful when SCF is getting stuck
  in a highly symmetric or oscillatory basin and needs a gentler symmetry break.

These soft modes are reproducible when you set `seed=...` and can be tuned
through `config["moment_assignment"]`, especially:

- `soft_moment_scale`
- `soft_min_candidate_moment`
- `soft_max_candidate_moment`

Example restart ladder for a difficult SCF:

```python
for strategy in ["afm_clusters", "afm_soft", "fm_soft", "balanced_random_soft"]:
    magmoms, feature_dict, graph, ring_info = assign_initial_magnetic_moments(
        atoms,
        strategy=strategy,
        cutoff=1.85,
        maxlength=12,
        k_hops=5,
        seed=0,
    )
    atoms.set_initial_magnetic_moments(magmoms)
    # Run your next DFT restart here.
```

## Package Layout

- `spinit.api`: high-level end-to-end API.
- `spinit.config`: defaults and config merge helpers.
- `spinit.graph`: graph/ring/geometry/topology utilities.
- `spinit.chemistry`: hybridization, motif detection, and per-atom feature construction.
- `spinit.seeding`: magnetic scoring and moment assignment strategies.
- `spinit.output`: reporting/summary utilities.


## Citation

If you use `spinit` in scientific work, please cite:

```bibtex
@software{spinit,
  author  = {Hart, Joseph},
  title   = {spinit: Physically informed spin initialisation for atomistic structures},
  year    = {2026},
  version = {0.1.0},
  url     = {https://github.com/joehart2001/spinit}
}
```

Full machine-readable metadata is available in [`CITATION.cff`](./CITATION.cff).
