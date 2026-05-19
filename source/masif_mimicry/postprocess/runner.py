"""Backward-compatible re-exports for postprocess entry points."""

from masif_mimicry.postprocess.discovery import discover_deduplicated_rows, select_representative_row
from masif_mimicry.postprocess.pipeline import process_results_mimicry

__all__ = ["discover_deduplicated_rows", "select_representative_row", "process_results_mimicry"]
