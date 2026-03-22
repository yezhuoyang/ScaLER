"""Compile ScaLERQEC IR to Qiskit QuantumCircuit with named registers.

Translates the intermediate representation produced by
:meth:`~scalerqec.QEC.qeccircuit.StabCode.construct_IR_standard_scheme`
into a ``qiskit.QuantumCircuit`` with separate ``data`` and ``anc``
quantum registers, suitable for layout-aware transpilation.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from qiskit import QuantumCircuit

from .qeccircuit import (
    DataMeasureInstruction,
    DetectorInstruction,
    IRInstruction,
    IRType,
    ObservableInstruction,
    StabPropInstruction,
)


class qiskitIRCompiler:
    """Compile ScaLERQEC IR instruction lists to Qiskit circuits.

    The output circuit uses two named ``QuantumRegister`` instances
    (``data`` and ``anc``) so that the caller can supply separate
    physical-qubit layouts for data and ancilla qubits.
    """

    def __init__(self) -> None:
        pass

    def compile_IR_to_qiskit(
        self,
        ir_list: list[IRInstruction],
        n_data: int,
        n_stabs: int,
        rounds: int,
    ) -> "QuantumCircuit":
        """Translate an IR instruction list into a Qiskit ``QuantumCircuit``.

        Args:
            ir_list: List of IR instruction objects from
                ``StabCode._IRList``.
            n_data: Number of data qubits (``StabCode._n``).
            n_stabs: Number of stabilizers (``len(StabCode._stabs)``).
            rounds: Number of syndrome extraction rounds.

        Returns:
            A ``qiskit.QuantumCircuit`` with registers ``data[0..n_data-1]``,
            ``anc[0..n_stabs-1]``, and a classical register ``c`` sized to
            hold all measurements.
        """
        from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister

        # Count total measurements to size the classical register.
        total_meas = 0
        for instr in ir_list:
            if instr._instr_type == IRType.PROP:
                total_meas += 1
            elif instr._instr_type == IRType.DATA_MEASURE:
                total_meas += n_data

        data_reg = QuantumRegister(n_data, name="data")
        anc_reg = QuantumRegister(n_stabs, name="anc")
        cl_reg = ClassicalRegister(total_meas, name="c")

        qc = QuantumCircuit(data_reg, anc_reg, cl_reg)

        meas_counter = 0
        # Map symbolic dest name -> classical bit index
        dest_to_clbit: dict[str, int] = {}
        current_round = -1

        for instr in ir_list:
            if instr._instr_type == IRType.PROP:
                assert isinstance(instr, StabPropInstruction)
                stab_idx = instr.get_stabindex()
                pauli_str = instr.stab
                dest = instr.dest

                # Barrier between rounds for readability
                if instr.round != current_round:
                    if current_round >= 0:
                        qc.barrier()
                    current_round = instr.round

                # Reset ancilla
                qc.reset(anc_reg[stab_idx])

                # Entangle data qubits with ancilla according to Pauli string
                for q, pauli in enumerate(pauli_str):
                    if pauli == "Z":
                        # Z-type: data is control, ancilla is target
                        qc.cx(data_reg[q], anc_reg[stab_idx])
                    elif pauli == "X":
                        # X-type: H-CX-H sandwich
                        qc.h(anc_reg[stab_idx])
                        qc.cx(anc_reg[stab_idx], data_reg[q])
                        qc.h(anc_reg[stab_idx])
                    elif pauli == "Y":
                        # Y = iXZ: use S-H-CX-H-Sdg sandwich
                        qc.s(anc_reg[stab_idx])
                        qc.h(anc_reg[stab_idx])
                        qc.cx(anc_reg[stab_idx], data_reg[q])
                        qc.h(anc_reg[stab_idx])
                        qc.sdg(anc_reg[stab_idx])
                    # I -> skip

                # Measure ancilla
                qc.measure(anc_reg[stab_idx], cl_reg[meas_counter])
                dest_to_clbit[dest] = meas_counter
                meas_counter += 1

            elif instr._instr_type == IRType.DATA_MEASURE:
                assert isinstance(instr, DataMeasureInstruction)
                qc.barrier()

                # Measure all data qubits
                for q in range(n_data):
                    dest = f"m{q}"
                    if instr.basis == "X":
                        qc.h(data_reg[q])
                    elif instr.basis == "Y":
                        qc.sdg(data_reg[q])
                        qc.h(data_reg[q])
                    qc.measure(data_reg[q], cl_reg[meas_counter])
                    dest_to_clbit[dest] = meas_counter
                    meas_counter += 1

            elif instr._instr_type == IRType.DETECTOR:
                assert isinstance(instr, DetectorInstruction)
                # Metadata only — add as a comment
                clbits = [
                    f"c[{dest_to_clbit[a]}]" for a in instr.args
                    if a in dest_to_clbit
                ]
                comment = f"DETECTOR {instr.dest}: XOR of {', '.join(clbits)}"
                # Attach metadata to circuit
                if not hasattr(qc, "_scalerqec_detectors"):
                    qc._scalerqec_detectors = []
                qc._scalerqec_detectors.append(
                    {"name": instr.dest, "clbits": [dest_to_clbit[a] for a in instr.args if a in dest_to_clbit]}
                )

            elif instr._instr_type == IRType.OBSERVABLE:
                assert isinstance(instr, ObservableInstruction)
                # Metadata only
                if not hasattr(qc, "_scalerqec_observables"):
                    qc._scalerqec_observables = []
                qc._scalerqec_observables.append(
                    {"name": instr.dest, "clbits": [dest_to_clbit[a] for a in instr.args if a in dest_to_clbit]}
                )

        return qc
