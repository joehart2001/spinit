"""Compatibility wrapper for :mod:`spinit.chemistry.motifs`."""

from .chemistry.motifs import (
    build_component_context,
    classify_atom_motif,
    detect_carbon_motif,
    detect_generic_motif,
    detect_hydrogen_motif,
    detect_nitrogen_motif,
    detect_oxygen_motif,
    get_component_status_for_atom,
    motif_effect,
)

__all__ = [
    "build_component_context",
    "get_component_status_for_atom",
    "motif_effect",
    "detect_oxygen_motif",
    "detect_nitrogen_motif",
    "detect_hydrogen_motif",
    "detect_carbon_motif",
    "detect_generic_motif",
    "classify_atom_motif",
]
