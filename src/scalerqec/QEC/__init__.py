"""QEC module for quantum error correction code construction and simulation.

This module provides the core abstractions for defining stabilizer-based quantum
error-correcting codes, attaching noise models, and compiling them into STIM
circuits suitable for logical error rate estimation.

Key components:
    - :class:`StabCode`: High-level stabilizer code representation parameterized
      by (n, k, d) that compiles stabilizer generators and logical operators into
      syndrome-extraction circuits.
    - :class:`NoiseModel`: Configurable depolarizing noise model that can be
      selectively applied to individual gate types.
    - Pre-built code instances (:func:`fivequbitCode`, :func:`steaneCode`,
      :func:`ShorCode`) for common small codes.
    - :class:`LogicalOperatorAnalyzer`: Utilities for analyzing logical operator
      structure of a stabilizer code.
"""

# Re-export high-level components for easy access

from .qeccircuit import StabCode
from .noisemodel import NoiseModel
from .small import fivequbitCode, steaneCode, ShorCode
from .analyzer import LogicalOperatorAnalyzer


__all__ = [
    "StabCode",
    "NoiseModel",
    "fivequbitCode",
    "steaneCode",
    "ShorCode",
    "LogicalOperatorAnalyzer",
]
