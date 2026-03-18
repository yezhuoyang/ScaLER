"""Symbolic computation of exact logical error rates (LER).

This module computes the **exact** logical error rate as a symbolic polynomial
in the physical error probability ``p`` using dynamic programming. The result
serves as ground truth for validating Monte Carlo and other approximate methods.

The key algorithm enumerates all possible detector/observable outcomes, weights
each by its probability polynomial (expressed in SymPy), and sums those
corresponding to logical errors. A DP recurrence over noise sources builds
these probabilities efficiently, avoiding the exponential cost of brute-force
enumeration over all :math:`4^N` Pauli error patterns.

Typical usage::

    from scalerqec.Symbolic import SymbolicLERcalc

    calc = SymbolicLERcalc()
    ler = calc.calculate_LER_from_file("path/to/stim_circuit", pvalue=0.001)
"""

# Re-export high-level components for easy access

from .symbolicLER import SymbolicLERcalc


__all__ = ["SymbolicLERcalc"]
