![spinit logo](spinit_logo.svg)

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
  - `afm_clusters`
  - `random_candidates`

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
from spinit import assign_initial_magnetic_moments

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

## Package Layout

- `spinit.api`: high-level end-to-end API.
- `spinit.config`: defaults and config merge helpers.
- `spinit.graph`: graph/ring/geometry/topology utilities.
- `spinit.chemistry`: hybridization, motif detection, and per-atom feature construction.
- `spinit.seeding`: magnetic scoring and moment assignment strategies.
- `spinit.output`: reporting/summary utilities.
- Compatibility modules (`spinit.bond_graph`, `spinit.graph_ops`, `spinit.hybridization`, etc.) still exist and re-export the new subpackage implementations.

## Development

```bash
pip install -e .[dev]
pytest
```

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
