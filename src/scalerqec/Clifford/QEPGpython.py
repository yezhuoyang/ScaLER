"""Quantum Error Propagation Graph (QEPG) construction in pure Python.

This module implements backward propagation of Pauli errors through a
Clifford circuit to determine how each single-qubit error at every noise
source affects the circuit's detectors and logical observable. The result
is an error propagation matrix that maps noise configurations to
detector/observable outcome vectors.
"""

import numpy as np

from .clifford import CliffordCircuit, Measurement, Reset, pauliNoise


class QEPGpython:
    """Quantum Error Propagation Graph built via backward Pauli propagation.

    Given a :class:`CliffordCircuit`, this class constructs a binary
    propagation matrix that encodes, for every noise source and every Pauli
    error type (X, Y, Z), which detectors and/or the logical observable are
    flipped.

    The propagation matrix ``_propMatrix`` has shape
    ``(3 * num_noise, num_detectors + 1)`` and is organized in three
    row blocks:

    - Rows ``[0, num_noise)``: X-error propagation for each noise source.
    - Rows ``[num_noise, 2*num_noise)``: Y-error propagation.
    - Rows ``[2*num_noise, 3*num_noise)``: Z-error propagation.

    The last column of each row indicates whether the corresponding error
    flips the logical observable.

    Args:
        circuit: A fully constructed :class:`CliffordCircuit` with gates,
            noise sources, measurements, detectors, and observable defined.
    """

    def __init__(self, circuit: CliffordCircuit):
        self._circuit = circuit
        self._total_meas = self._circuit._totalMeas
        self._total_noise = self._circuit._totalnoise
        self._propMatrix = np.zeros(
            (3 * self._total_noise, len(self._circuit.parityMatchGroup) + 1),
            dtype="uint8",
        )

    def backword_graph_construction(self) -> None:
        """Build the error propagation matrix by reverse-traversing the circuit.

        Processes the gate list from last to first. Three per-qubit tracking
        arrays (``current_x_prop``, ``current_y_prop``, ``current_z_prop``)
        record how a Pauli X, Y, or Z error on each qubit would propagate
        to detectors and the observable at the current point in the reverse
        traversal.

        Gate-type rules applied during backward propagation:

        - **Measurement**: An X or Y error on the measured qubit flips the
          detectors that include this measurement (and the observable if
          applicable). Z errors commute with Z-basis measurement and have
          no effect.
        - **Reset**: Clears all propagation on the reset qubit, since the
          qubit state is discarded.
        - **CNOT**: X propagates forward from control to target; Z propagates
          backward from target to control. Y propagation combines both
          effects.
        - **Hadamard**: Swaps X and Z propagation on the gate's qubit.
        - **Noise source**: Snapshots the current propagation state for this
          qubit into the corresponding rows of ``_propMatrix``.

        After this method completes, ``_propMatrix`` is fully populated and
        ready for querying via :meth:`sample_x_error`, :meth:`sample_y_error`,
        :meth:`sample_z_error`, or :meth:`sample_noise_vector`.
        """
        nqubit = self._circuit._qubit_num

        column_size = len(self._circuit.parityMatchGroup) + 1
        current_x_prop = np.zeros((nqubit, column_size), dtype="uint8")
        current_y_prop = np.zeros((nqubit, column_size), dtype="uint8")
        current_z_prop = np.zeros((nqubit, column_size), dtype="uint8")

        current_noise_index = self._circuit._totalnoise - 1
        current_meas_index = self._total_meas - 1

        total_noise = self._total_noise
        T = len(self._circuit._gatelists)

        for t in range(T - 1, -1, -1):
            # Update current_x_prop, current_y_prop, current_z_prop based on the current gate and measurement
            gate = self._circuit._gatelists[t]
            """
            If the gate is a oiginal noise, add edges to the graph based on current propogation
            """
            if isinstance(gate, pauliNoise):
                noiseindex = current_noise_index
                # print("Noise!")
                self._propMatrix[noiseindex, :] = current_x_prop[
                    gate._qubitindex, :
                ].copy()
                self._propMatrix[total_noise + noiseindex, :] = current_y_prop[
                    gate._qubitindex, :
                ].copy()
                self._propMatrix[total_noise * 2 + noiseindex, :] = current_z_prop[
                    gate._qubitindex, :
                ].copy()
                current_noise_index -= 1
                continue
            """
            When there is a measurement, update the current propogation based on the measurement
            We just need to consider the propagation of X and Y because only 
            the X and Y error can be detected by the measurement
            """
            if isinstance(gate, Measurement):
                measureindex = current_meas_index

                qindex = gate._qubitindex
                if measureindex in self._circuit.observable:
                    current_x_prop[qindex][column_size - 1] = 1
                    current_y_prop[qindex][column_size - 1] = 1

                for parityIdx in self._circuit.get_measIdx_to_parityIdx(measureindex):
                    current_x_prop[qindex][parityIdx] = 1
                    current_y_prop[qindex][parityIdx] = 1

                current_meas_index -= 1
                continue

            if isinstance(gate, Reset):
                current_x_prop[gate._qubitindex, :] = 0
                current_y_prop[gate._qubitindex, :] = 0
                current_z_prop[gate._qubitindex, :] = 0
                continue

            """
            Deal with propagation by CNOT gate, we need to consider the propagation of X and Z
            """
            if gate._name == "CNOT":
                control = gate._control
                target = gate._target
                current_x_prop[control, :] = (
                    current_x_prop[control, :] + current_x_prop[target, :]
                ) % 2
                current_z_prop[target, :] = (
                    current_z_prop[control, :] + current_z_prop[target, :]
                ) % 2
                current_y_prop[control, :] = (
                    current_y_prop[control, :] + current_x_prop[target, :]
                ) % 2
                current_y_prop[target, :] = (
                    current_y_prop[target, :] + current_z_prop[control, :]
                ) % 2
                continue

            """
            Deal with propagation by H gate
            If there is a H gate, we need to swap the X and Z propagations
            """
            if gate._name == "H":
                qubitindex = gate._qubitindex
                tmp_row = current_x_prop[qubitindex, :].copy()
                current_x_prop[qubitindex, :] = current_z_prop[qubitindex, :]
                current_z_prop[qubitindex, :] = tmp_row
                continue

    def sample_x_error(self, noise_index: int) -> list[int]:
        """Get the detector/observable flip vector for an X error at a noise source.

        Args:
            noise_index: The noise source index (0-based).

        Returns:
            A list of length ``num_detectors + 1`` with binary entries.
            A ``1`` at position *i* means detector *i* is flipped; a ``1``
            in the last position means the observable is flipped.
        """
        return list(self._propMatrix[noise_index, :])

    def sample_y_error(self, noise_index: int) -> list[int]:
        """Get the detector/observable flip vector for a Y error at a noise source.

        Args:
            noise_index: The noise source index (0-based).

        Returns:
            A list of length ``num_detectors + 1`` with binary entries.
        """
        return list(self._propMatrix[self._total_noise + noise_index, :])

    def sample_z_error(self, noise_index: int) -> list[int]:
        """Get the detector/observable flip vector for a Z error at a noise source.

        Args:
            noise_index: The noise source index (0-based).

        Returns:
            A list of length ``num_detectors + 1`` with binary entries.
        """
        return list(self._propMatrix[2 * self._total_noise + noise_index, :])

    def sample_noise_vector(self, noise_vector: np.ndarray) -> np.ndarray:
        """Compute the detector/observable outcome for an arbitrary error pattern.

        Multiplies a binary noise vector (indicating which X/Y/Z errors are
        active at which noise sources) by the propagation matrix to obtain
        the combined detector and observable flip pattern.

        Args:
            noise_vector: A binary array of length ``3 * num_noise``. The
                first ``num_noise`` entries correspond to X errors, the next
                to Y errors, and the last to Z errors, matching the row
                layout of ``_propMatrix``.

        Returns:
            An array of length ``num_detectors + 1`` giving the XOR-combined
            detector flips and observable flip.

        Raises:
            AssertionError: If ``noise_vector`` length does not equal
                ``3 * num_noise``.
        """
        assert len(noise_vector) == 3 * self._total_noise
        return noise_vector @ self._propMatrix


if __name__ == "__main__":
    stim_str = ""
    filepath = "C:/Users/username/Documents/Sampling/stimprograms/1cnot"
    with open(filepath, "r", encoding="utf-8") as f:
        stim_str = f.read()
    circuit = CliffordCircuit(3)
    circuit.set_error_rate(0.01)
    circuit.compile_from_stim_circuit_str(stim_str)

    graph = QEPGpython(circuit)
    graph.backword_graph_construction()

    print(graph.sample_noise_vector(np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])))
