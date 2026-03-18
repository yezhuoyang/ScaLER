"""Depolarizing noise model for quantum error correction circuits.

This module provides :class:`NoiseModel`, which injects single-qubit
depolarizing noise before each enabled gate type in a
:class:`~scalerqec.Clifford.clifford.CliffordCircuit`.  Individual error
channels can be selectively disabled via :meth:`NoiseModel.disable_error`.
"""

from enum import Enum

from ..Clifford.clifford import (
    CliffordCircuit,
    Measurement,
    Reset,
    SingleQGate,
    TwoQGate,
)


class ErrorType(Enum):
    """Enumeration of gate/operation types that can carry noise.

    Each member corresponds to a physical operation in the circuit.
    Used with :meth:`NoiseModel.disable_error` to selectively silence
    individual noise channels.

    Attributes:
        MEASUREMENT: Measurement operation.
        RESET: Qubit reset operation.
        CNOT: Two-qubit CNOT gate.
        HADAMARD: Hadamard gate.
        PHASE: Phase (S) gate.
        PAULIX: Pauli-X gate.
        PAULIY: Pauli-Y gate.
        PAULIZ: Pauli-Z gate.
        CZ: Two-qubit controlled-Z gate.
    """

    MEASUREMENT = 0
    RESET = 1
    CNOT = 2
    HADAMARD = 3
    PHASE = 4
    PAULIX = 5
    PAULIY = 6
    PAULIZ = 7
    CZ = 8


class NoiseModel:
    """Configurable depolarizing noise model.

    Wraps a single depolarizing error rate and a set of per-gate-type
    enable flags.  When applied to a
    :class:`~scalerqec.Clifford.clifford.CliffordCircuit` via
    :meth:`reconstruct_clifford_circuit`, a single-qubit depolarizing
    channel is inserted immediately before every enabled gate.

    By default all gate types are enabled.  Call
    :meth:`disable_error` to turn off noise for specific operations.

    Args:
        error_rate: Probability of depolarizing error per gate, applied
            uniformly to all enabled gate types.
    """

    def __init__(self, error_rate: float) -> None:
        self._error_rate = error_rate

        self._has_MEASUREMENT_error = True
        self._has_RESET_error = True
        self._has_CNOT_error = True
        self._has_CZ_error = True
        self._has_HADAMARD_error = True
        self._has_PHASE_error = True
        self._has_PAULIX_error = True
        self._has_PAULIY_error = True
        self._has_PAULIZ_error = True

    @property
    def error_rate(self) -> float:
        """
        Get the error rate of the noise model.

        Returns:
            The error rate as a float.
        """
        return self._error_rate

    @error_rate.setter
    def error_rate(self, value: float) -> None:
        """
        Set the error rate of the noise model.

        Args:
            value: The new error rate as a float.
        """
        self._error_rate = value

    def disable_error(self, error_type: str) -> None:
        """Disable noise injection for a specific gate type.

        After calling this method, the corresponding gate will no longer
        have a depolarizing channel inserted before it during circuit
        reconstruction.

        Args:
            error_type: String identifier of the gate type to silence.
                Accepted values: ``"MEASUREMENT"``, ``"RESET"``,
                ``"CNOT"``, ``"CZ"``, ``"H"``, ``"P"``, ``"X"``,
                ``"Y"``, ``"Z"``.
        """
        if error_type == "MEASUREMENT":
            self._has_MEASUREMENT_error = False
        elif error_type == "RESET":
            self._has_RESET_error = False
        elif error_type == "CNOT":
            self._has_CNOT_error = False
        elif error_type == "CZ":
            self._has_CZ_error = False
        elif error_type == "H":
            self._has_HADAMARD_error = False
        elif error_type == "P":
            self._has_PHASE_error = False
        elif error_type == "X":
            self._has_PAULIX_error = False
        elif error_type == "Y":
            self._has_PAULIY_error = False
        elif error_type == "Z":
            self._has_PAULIZ_error = False

    def rewrite_stim_program(self, stim_program: str) -> str:
        """Rewrite a STIM program string to incorporate depolarizing noise.

        .. note:: Not yet implemented; returns the input unchanged.

        Args:
            stim_program: The original STIM program as a string.

        Returns:
            A new STIM program string with depolarizing noise inserted
            before each enabled gate type.
        """
        # Placeholder for actual implementation
        return stim_program

    def reconstruct_clifford_circuit(
        self, clifford_circuit: CliffordCircuit
    ) -> CliffordCircuit:
        """Rebuild a Clifford circuit with depolarizing noise injected.

        Creates a new :class:`~scalerqec.Clifford.clifford.CliffordCircuit`
        by iterating over every gate in *clifford_circuit*.  For each
        gate whose type is enabled, a single-qubit depolarizing channel
        (at the stored :attr:`error_rate`) is inserted on the relevant
        qubit(s) immediately before the gate.  For two-qubit gates
        (CNOT, CZ), depolarizing noise is applied independently to both
        the control and target qubits.

        Detector and observable metadata are copied from the original
        circuit and recompiled on the new circuit.

        Args:
            clifford_circuit: The noise-free Clifford circuit to
                augment.

        Returns:
            A new :class:`~scalerqec.Clifford.clifford.CliffordCircuit`
            instance containing the same logical operations with
            depolarizing noise inserted.
        """
        # Placeholder for actual implementation

        num_qubits = clifford_circuit.qubitnum
        new_circuit = CliffordCircuit(num_qubits)

        new_circuit.error_rate = self._error_rate
        gate_list = clifford_circuit.gatelists

        for gate in gate_list:
            if isinstance(gate, TwoQGate):
                # Apply CNOT gate with noise
                if gate.name == "CNOT":
                    if self._has_CNOT_error:
                        new_circuit.add_depolarize(gate.control)
                        new_circuit.add_depolarize(gate.target)
                    new_circuit.add_cnot(gate.control, gate.target)
                elif gate.name == "CZ":
                    if self._has_CZ_error:
                        new_circuit.add_depolarize(gate.control)
                        new_circuit.add_depolarize(gate.target)
                    new_circuit.add_cz(gate.control, gate.target)

            elif isinstance(gate, SingleQGate):
                # Apply single-qubit gate with noise
                if gate.name == "H":
                    if self._has_HADAMARD_error:
                        new_circuit.add_depolarize(gate.qubitindex)
                    new_circuit.add_hadamard(gate.qubitindex)
                elif gate.name == "P":
                    if self._has_PHASE_error:
                        new_circuit.add_depolarize(gate.qubitindex)
                    new_circuit.add_phase(gate.qubitindex)
                elif gate.name == "X":
                    if self._has_PAULIX_error:
                        new_circuit.add_depolarize(gate.qubitindex)
                    new_circuit.add_paulix(gate.qubitindex)
                elif gate.name == "Y":
                    if self._has_PAULIY_error:
                        new_circuit.add_depolarize(gate.qubitindex)
                    new_circuit.add_pauliy(gate.qubitindex)
                elif gate.name == "Z":
                    if self._has_PAULIZ_error:
                        new_circuit.add_depolarize(gate.qubitindex)
                    new_circuit.add_pauliz(gate.qubitindex)

            elif isinstance(gate, Measurement):
                # Apply measurement with noise
                if self._has_MEASUREMENT_error:
                    new_circuit.add_depolarize(gate.qubitindex)
                new_circuit.add_measurement(gate.qubitindex)
                # Placeholder for adding noise after measurement

            elif isinstance(gate, Reset):
                # Apply reset with noise
                if self._has_RESET_error:
                    new_circuit.add_depolarize(gate.qubitindex)
                new_circuit.add_reset(gate.qubitindex)
                # Placeholder for adding noise after reset

        new_circuit.parityMatchGroup = clifford_circuit.parityMatchGroup
        new_circuit.observable = clifford_circuit.observable
        new_circuit.compile_detector_and_observable()

        return new_circuit

    def uniform_depolarization_single(stim_program: str) -> str:
        """Apply uniform depolarization to single-qubit gates in a STIM program.

        .. note:: Not yet implemented; returns the input unchanged.

        Args:
            stim_program: The original STIM program as a string.

        Returns:
            The modified STIM program with single-qubit depolarization
            noise applied.
        """
        return stim_program
