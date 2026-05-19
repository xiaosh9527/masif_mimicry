"""Post-processing metrics and search-result aggregation."""

from masif_mimicry.postprocess.discovery import discover_deduplicated_rows
from masif_mimicry.postprocess.pipeline import process_results_mimicry

__all__ = ["discover_deduplicated_rows", "process_results_mimicry"]
