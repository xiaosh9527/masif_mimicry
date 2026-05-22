"""Load and merge domain metadata CSV (TED / domainome) for postprocess enrichment."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd

def load_domain_metadata_filtered(
    csv_path: Path | str,
    p1_ids: Iterable[str],
    *,
    id_column: str = "id",
    chunk_size: int = 100_000,
) -> pd.DataFrame:
    """
    Stream-read a large domainome CSV and keep only rows whose id is in p1_ids.

    Keeps all metadata columns for matched ids, while avoiding loading the full CSV at once.
    """
    csv_path = Path(csv_path)
    wanted = {str(x) for x in p1_ids}
    if not wanted:
        return pd.DataFrame()

    chunks: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        csv_path,
        chunksize=chunk_size,
        low_memory=False,
    ):
        sub = chunk[chunk[id_column].astype(str).isin(wanted)]
        if not sub.empty:
            chunks.append(sub)
    if not chunks:
        return pd.DataFrame()
    return pd.concat(chunks, ignore_index=True)


def merge_metadata_columns(match_info: dict, meta_row: pd.Series | None, *, id_column: str = "id") -> None:
    """
    Add metadata fields to match_info in place. Database columns win on name clashes
    (same semantics as tricomplex gather_tables drop + merge), except id is skipped.
    """
    if meta_row is None:
        return
    for col in meta_row.index:
        if col == id_column:
            continue
        val = meta_row[col]
        if col in match_info:
            del match_info[col]
        match_info[col] = val
