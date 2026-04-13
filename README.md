![spinit logo](spinit_logo.svg)

# spinit

`spinit` assigns **initial magnetic moments** for spin-polarized DFT calculations using local chemistry, topology, and geometry heuristics on `ase.Atoms` structures.

It is designed for carbon-rich systems (graphenic, nanoporous, amorphous carbon) and mixed chemistries containing heteroatoms (H/O/N/other).

## Scientific Scope

`spinit` is a **heuristic seed generator** for SCF initialization.

- It does not prove magnetic ground states.
- It provides physically motivated initial spin patterns to improve exploration of plausible spin solutions.
- It combines undercoordination, local geometry/strain, ring topology, and explicit motif detection.

## Features

- PBC-aware bond graph construction (ASE + MIC vectors).
- Ring/topology analysis from local graph cycles.
- Hybridization-like local classification (`sp`, `sp2`, `sp3`, distorted variants).
- Explicit motif detection for O/N/H/C-centered local environments.
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

- `spinit.config`: defaults and config merge helpers.
- `spinit.bond_graph`: graph/ring primitives.
- `spinit.graph_ops`: graph construction + MIC geometry + local graph accessors.
- `spinit.geometry`: vector/angle utilities.
- `spinit.hybridization`: local chemistry + sp2 strain logic.
- `spinit.motifs`: connectivity and motif detection.
- `spinit.features`: per-atom feature construction.
- `spinit.scoring`: magnetic scoring and candidate selection.
- `spinit.assignment`: magnetic moment sign/amplitude assignment strategies.
- `spinit.reporting`: summary utilities.
- `spinit.api`: high-level end-to-end API.

## Development

```bash
pip install -e .[dev]
pytest
```

## Citation

If you use `spinit` in scientific work, please cite the software metadata in [`CITATION.cff`](./CITATION.cff).
