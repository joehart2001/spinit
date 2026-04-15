"""Reporting utilities for feature and candidate summaries."""

from __future__ import annotations

from typing import Any, Mapping, Sequence

import numpy as np


def make_feature_table(feature_dict: Mapping[int, Mapping[str, Any]]) -> Any:
    """Return DataFrame when pandas exists, else list of row dictionaries."""
    rows = [dict(feature_dict[i]) for i in sorted(feature_dict)]
    try:
        import pandas as pd
    except ImportError:
        return rows
    return pd.DataFrame(rows)


def print_feature_summary(feature_dict: Mapping[int, Mapping[str, Any]], max_rows: int = 20) -> None:
    """Print compact feature summary for quick inspection."""
    table = make_feature_table(feature_dict)
    print("\n=== Atom Feature Summary ===")
    if hasattr(table, "head"):
        print(table.head(max_rows).to_string(index=False))
    else:
        for row in table[:max_rows]:
            print(row)
    print("============================\n")


def print_candidate_summary(feature_dict: Mapping[int, Mapping[str, Any]], magmoms: Sequence[float]) -> None:
    """Print candidate summary from non-zero assigned initial moments."""
    magmoms_arr = np.asarray(magmoms, dtype=float)
    active = [i for i in range(len(magmoms_arr)) if abs(magmoms_arr[i]) > 1e-12]

    print("\n=== Candidate Summary ===")
    print(f"Candidate count: {len(active)}")

    if len(active) == 0:
        print("No candidate atoms selected.")
        print("=========================\n")
        return

    active_sorted = sorted(active, key=lambda i: abs(magmoms_arr[i]), reverse=True)
    for i in active_sorted[:20]:
        f = feature_dict[i]
        print(
            f"i={i:4d} | elem={f.get('element')} | m0={magmoms_arr[i]:+6.3f} | "
            f"score={f.get('magnetic_score', 0.0):6.3f} | "
            f"intr={f.get('intrinsic_score', 0.0):6.3f} | mod={f.get('neighbor_modifier', 0.0):+6.3f} | "
            f"coord={f.get('coordination')} | hyb={f.get('hybridization')}"
        )
    print("=========================\n")
