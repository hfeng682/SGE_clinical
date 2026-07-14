"""
sge_pipeline — reproduce the raw-MaveDB -> annotated transform for two-score
saturation-genome-editing screens.

    from sge_pipeline.pipeline import transform
    transform("BARD1", "raw/BARD1.csv", "out/BARD1_annotated.csv")

See README.md for stages, output schema, and how to add a new screen.
"""
from .pipeline import transform          # noqa: F401

__all__ = ["transform"]
__version__ = "1.0"
