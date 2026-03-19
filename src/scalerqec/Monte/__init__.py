"""Monte Carlo sampling module for logical error rate estimation.

This module provides Monte Carlo methods for estimating the logical error
rate (LER) of quantum error-correcting codes. It supports multiple sampling
backends:

* **MonteLERcalc** -- Adaptive Monte Carlo sampling using stim detector
  sampling with pymatching decoder, QEPG-accelerated sampling via a compiled
  C++ backend, or sinter-based distributed sampling.
* **MonteLDPC** -- Monte Carlo sampling tailored for LDPC codes.

All estimators use adaptive batching to efficiently allocate samples: batch
sizes grow when few errors have been observed and stop once a minimum number
of logical error events has been reached or a sample/time budget is exhausted.
"""

# Re-export high-level components for easy access

from .monteLER import MonteLERcalc
from .monteLDPC import MonteLDPC
from .noise_model_parser import (
    extract_noise_model,
    NonuniformNoiseModel,
    CorrelatedNoisePair,
)


__all__ = [
    "MonteLERcalc",
    "MonteLDPC",
    "extract_noise_model",
    "NonuniformNoiseModel",
    "CorrelatedNoisePair",
]
