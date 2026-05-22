"""
Aggregate interface metrics from per-domain annotations (resi, pLDDT, DeepTMHMM, DSSP).

Ported from tricomplex-design ``src.utils`` (iface_metrics / iface_sse_metrics / get_location_masks).
Expects merged CSV columns on each row: resi, plddt, deeptmhmm_annotation, mdtraj_dssp_annotation,
and interface residue list in ``matched_iface_resi``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from masif_mimicry.postprocess.metrics.secondary_structure import find_sse


def get_location_masks(row, n_res: int | None = None) -> tuple:
    """Return (intracellular, extracellular, transmembrane, unassigned) boolean arrays or Nones."""
    ann = row.get("deeptmhmm_annotation") if isinstance(row, dict) else getattr(row, "deeptmhmm_annotation", None)
    if not isinstance(ann, str):
        return None, None, None, None
    if n_res is not None and len(ann) != n_res:
        return None, None, None, None

    intracellular = np.array([label in {"S", "I"} for label in ann])
    extracellular = np.array([label in {"O"} for label in ann])
    transmembrane = np.array([label in {"M", "B"} for label in ann])
    unassigned = np.array([label in {"P", "?"} for label in ann])
    return intracellular, extracellular, transmembrane, unassigned


def iface_sse_metrics(row, iface_resi: set) -> tuple:
    """SSE counts/fractions at interface from mdtraj_dssp_annotation."""
    md = row.get("mdtraj_dssp_annotation") if isinstance(row, dict) else getattr(row, "mdtraj_dssp_annotation", None)
    if not isinstance(md, str):
        return None, None, None, None, None

    resi_raw = row.get("resi") if isinstance(row, dict) else getattr(row, "resi", None)
    if not isinstance(resi_raw, str):
        return None, None, None, None, None

    residue_ids = list(map(int, resi_raw.split(",")))
    mdtraj_dssp_labels = list(md)
    if len(residue_ids) != len(mdtraj_dssp_labels):
        return None, None, None, None, None

    iface_mdtraj_dssp_labels = [x for x, r in zip(mdtraj_dssp_labels, residue_ids) if r in iface_resi]
    segments = find_sse(residue_ids, mdtraj_dssp_labels)
    iface_segments = [seg for seg in segments if any(seg["start"] <= r <= seg["end"] for r in iface_resi)]

    total_n_sse = len(iface_segments)
    n_helix = sum(s["label"] == "H" for s in iface_segments)
    n_strand = sum(s["label"] == "E" for s in iface_segments)
    helix_frac = float(np.mean([x == "H" for x in iface_mdtraj_dssp_labels])) if iface_mdtraj_dssp_labels else None
    strand_frac = float(np.mean([x == "E" for x in iface_mdtraj_dssp_labels])) if iface_mdtraj_dssp_labels else None
    return total_n_sse, n_helix, n_strand, helix_frac, strand_frac


def iface_metrics(row, iface_col: str = "matched_iface_resi") -> tuple:
    """
    Return (
        iface_intracellular_frac, iface_extracellular_frac, iface_transmembrane_frac,
        iface_plddt, iface_n_sse, iface_n_helix, iface_n_strand, iface_helix_frac, iface_strand_frac,
    ). NaN / None when inputs are missing or invalid.
    """
    nan9 = (np.nan,) * 9

    resi_raw = row.get("resi") if isinstance(row, dict) else getattr(row, "resi", None)
    if not isinstance(resi_raw, str):
        return nan9

    resi = np.array(list(map(int, resi_raw.split(","))))
    iface_val = row.get(iface_col) if isinstance(row, dict) else getattr(row, iface_col, None)
    if not isinstance(iface_val, str):
        return nan9

    iface_resi = set(map(int, iface_val.split(",")))
    if len(iface_resi) == 0:
        return nan9

    iface_mask = np.array([r in iface_resi for r in resi])

    plddt_raw = row.get("plddt") if isinstance(row, dict) else getattr(row, "plddt", None)
    if not isinstance(plddt_raw, str):
        return nan9
    plddt = np.array(list(map(float, plddt_raw.split(","))))
    if len(plddt) != len(resi):
        return nan9

    if not iface_mask.any():
        iface_plddt = np.nan
    else:
        iface_plddt = float(np.mean(plddt[iface_mask]))

    intracellular_mask, extracellular_mask, transmembrane_mask, _ = get_location_masks(row, n_res=len(resi))

    if intracellular_mask is not None:
        if len(resi) != len(intracellular_mask):
            iface_intracellular_frac = np.nan
        else:
            intracellular_resi = set(resi[intracellular_mask])
            iface_intracellular_frac = len(iface_resi & intracellular_resi) / len(iface_resi)
    else:
        iface_intracellular_frac = np.nan

    if extracellular_mask is not None:
        if len(resi) != len(extracellular_mask):
            iface_extracellular_frac = np.nan
        else:
            extracellular_resi = set(resi[extracellular_mask])
            iface_extracellular_frac = len(iface_resi & extracellular_resi) / len(iface_resi)
    else:
        iface_extracellular_frac = np.nan

    if transmembrane_mask is not None:
        if len(resi) != len(transmembrane_mask):
            iface_transmembrane_frac = np.nan
        else:
            transmembrane_resi = set(resi[transmembrane_mask])
            iface_transmembrane_frac = len(iface_resi & transmembrane_resi) / len(iface_resi)
    else:
        iface_transmembrane_frac = np.nan

    iface_n_sse, iface_n_helix, iface_n_strand, iface_helix_frac, iface_strand_frac = iface_sse_metrics(
        row, iface_resi
    )

    return (
        iface_intracellular_frac,
        iface_extracellular_frac,
        iface_transmembrane_frac,
        iface_plddt,
        iface_n_sse,
        iface_n_helix,
        iface_n_strand,
        iface_helix_frac,
        iface_strand_frac,
    )


IFACE_AGGREGATE_COLUMNS = [
    "iface_intracellular_frac",
    "iface_extracellular_frac",
    "iface_transmembrane_frac",
    "iface_plddt",
    "iface_n_sse",
    "iface_n_helix",
    "iface_n_strand",
    "iface_helix_frac",
    "iface_strand_frac",
]


def add_iface_aggregate_metrics(match_info: dict) -> None:
    """Update match_info in place with iface aggregate columns; missing data yields NaN/None."""
    row = pd.Series(match_info)
    vals = iface_metrics(row, iface_col="matched_iface_resi")
    for col, v in zip(IFACE_AGGREGATE_COLUMNS, vals):
        match_info[col] = v
