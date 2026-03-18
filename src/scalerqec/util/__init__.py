"""General-purpose utilities for ScaLERQEC.

This package provides helper functions used throughout the ScaLERQEC
framework:

- **Binomial utilities** (:mod:`~scalerqec.util.binomial`):
  Compute binomial weights and subspace sizes for error-weight
  analysis.
- **Pauli algebra** (:mod:`~scalerqec.util.pauli`):
  Commutation checks and multiplication of Pauli strings.
- **Formatting** (:mod:`~scalerqec.util.printer`):
  Scientific notation with uncertainty for reporting results.
"""

# Re-export high-level components for easy access

from .binomial import binomial_weight, subspace_size
from .printer import format_with_uncertainty
from .pauli import commute, anticommute, multiply_pauli_strings

__all__ = [
    "binomial_weight",
    "subspace_size",
    "format_with_uncertainty",
    "commute",
    "anticommute",
    "multiply_pauli_strings",
]
