"""QEC module for quantum error correction code construction and simulation.

This module provides the core abstractions for defining stabilizer-based quantum
error-correcting codes, attaching noise models, and compiling them into STIM
circuits suitable for logical error rate estimation.

Key components:
    - :class:`StabCode`: High-level stabilizer code representation parameterized
      by (n, k, d) that compiles stabilizer generators and logical operators into
      syndrome-extraction circuits.
    - :class:`RepetitionCode`, :class:`SurfaceCode`: Parametric code families.
    - :class:`NoiseModel`, :class:`SD6NoiseModel`, :class:`SI1000NoiseModel`,
      :class:`SIDNoiseModel`: Noise models for circuit-level simulation.
    - Pre-built code instances (:func:`fivequbitCode`, :func:`steaneCode`,
      :func:`ShorCode`) for common small codes.
    - :class:`IRAnalyzer`, :class:`LogicalOperatorAnalyzer`: Static analysis
      and verification utilities.
"""

# Re-export high-level components for easy access

from .qeccircuit import StabCode
from .noisemodel import NoiseModel, SD6NoiseModel, SI1000NoiseModel, SIDNoiseModel
from .surface import RepetitionCode, SurfaceCode
from .small import fivequbitCode, steaneCode, ShorCode
from .analyzer import IRAnalyzer, LogicalOperatorAnalyzer, DistanceAnalyzer


__all__ = [
    "StabCode",
    "RepetitionCode",
    "SurfaceCode",
    "NoiseModel",
    "SD6NoiseModel",
    "SI1000NoiseModel",
    "SIDNoiseModel",
    "IRAnalyzer",
    "LogicalOperatorAnalyzer",
    "DistanceAnalyzer",
    "fivequbitCode",
    "steaneCode",
    "ShorCode",
]
