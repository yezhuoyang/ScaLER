"""Clifford circuit representation and error propagation analysis.

This module provides the core abstractions for representing Clifford quantum
circuits and analyzing how Pauli errors propagate through them. It supports:

- Building Clifford circuits from individual gate operations or by parsing
  STIM circuit format strings.
- Representing single-qubit gates (H, S, X, Y, Z), two-qubit gates (CNOT, CZ),
  Pauli noise channels, measurements, and resets.
- Constructing Quantum Error Propagation Graphs (QEPG) that track how X, Y,
  and Z errors at each noise source affect detector and observable outcomes.

The typical workflow is:

1. Create a :class:`CliffordCircuit` and populate it with gates (manually or
   via :meth:`CliffordCircuit.compile_from_stim_circuit_str`).
2. Pass the circuit to :class:`~scalerqec.Clifford.QEPGpython.QEPGpython` and
   call :meth:`~scalerqec.Clifford.QEPGpython.QEPGpython.backword_graph_construction`
   to build the error propagation matrix.
3. Query the propagation matrix to determine which detectors and observables
   are flipped by each possible single-qubit error.
"""

# Re-export high-level components for easy access

from .clifford import CliffordCircuit
from .noiselabel import NoiseLabel, NoiseLabelMap


__all__ = ["CliffordCircuit", "NoiseLabel", "NoiseLabelMap"]
