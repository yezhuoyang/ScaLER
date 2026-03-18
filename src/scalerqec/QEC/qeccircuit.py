"""Quantum error correction circuit construction and compilation.

This module provides the core :class:`StabCode` class for representing
stabilizer-based quantum error-correcting codes and compiling them into
STIM circuits via a two-stage pipeline:

1. **IR construction** -- stabilizer generators and logical operators are
   translated into an intermediate representation (IR) consisting of
   :class:`StabPropInstruction`, :class:`DetectorInstruction`, and
   :class:`ObservableInstruction` objects.
2. **Circuit compilation** -- the IR is lowered into a
   :class:`~scalerqec.Clifford.clifford.CliffordCircuit` and then into a
   ``stim.Circuit``.

The module also defines the :class:`SCHEME` and :class:`IRType` enumerations
and the full IR instruction hierarchy used internally by the compiler.
"""

from __future__ import annotations

from enum import Enum
from typing import Optional

import numpy as np
from numpy.typing import NDArray
import stim

from scalerqec.Clifford.clifford import CliffordCircuit
from scalerqec.QEC.noisemodel import NoiseModel
from scalerqec.util import commute


class SCHEME(Enum):
    """Syndrome-extraction schemes for stabilizer codes.

    Each member selects a different strategy for measuring stabilizer
    generators and propagating parity information.

    Attributes:
        STANDARD: Bare ancilla syndrome extraction (one ancilla per
            stabilizer generator, no verification).
        SHOR: Shor-style fault-tolerant extraction using cat states.
        KNILL: Knill-style extraction via teleportation.
        FLAG: Flag-based fault-tolerant extraction with flag qubits.
    """

    STANDARD = 0
    SHOR = 1
    KNILL = 2
    FLAG = 3


class IRType(Enum):
    """Types of intermediate-representation (IR) instructions.

    The IR captures the logical structure of a syndrome-extraction circuit
    before it is lowered to physical gates.

    Attributes:
        PROP: Stabilizer propagation (ancilla preparation, entanglement,
            and measurement).
        DETECTOR: Detector declaration comparing measurement outcomes
            across rounds.
        OBSERVABLE: Logical observable declaration.
        IF_THEN: Conditional branch (not yet implemented).
        WHILE: While-loop control flow (not yet implemented).
        REPEAT_UNTIL: Repeat-until control flow (not yet implemented).
        REPEAT: Fixed-count repetition (not yet implemented).
    """

    PROP = 0
    DETECTOR = 1
    OBSERVABLE = 2
    IF_THEN = 3
    WHILE = 4
    REPEAT_UNTIL = 5
    REPEAT = 6


class IRInstruction:
    """Base class for all intermediate-representation (IR) instructions.

    Every IR instruction carries a type tag (:class:`IRType`) that
    identifies the kind of operation it represents.  Concrete subclasses
    add the operand-specific fields.

    Args:
        instr_type: The IR instruction type tag.
    """

    def __init__(self, instr_type: IRType) -> None:
        self._instr_type = instr_type


class StabPropInstruction(IRInstruction):
    """IR instruction for stabilizer propagation (syndrome measurement).

    Represents the preparation of an ancilla qubit, entanglement with data
    qubits according to a Pauli stabilizer string, and subsequent
    measurement.  The result is stored in a symbolic destination register.

    Args:
        round: The syndrome-extraction round index (0-based).
        stabindex: Index of the stabilizer generator within the code's
            stabilizer list.
        dest: Symbolic name for the measurement result (e.g., ``"c0"``).
        stab: Pauli string describing the stabilizer generator, composed
            of characters ``I``, ``X``, ``Y``, ``Z``.
        is_observable: If ``True``, this instruction propagates a logical
            observable rather than a stabilizer generator.
        observable_index: Index of the logical observable when
            *is_observable* is ``True``; ``-1`` otherwise.
    """

    def __init__(
        self,
        round: int,
        stabindex: int,
        dest: str,
        stab: str,
        is_observable: bool = False,
        observable_index: int = -1,
    ) -> None:
        super().__init__(IRType.PROP)
        self._round = round
        self._stabindex = stabindex
        self._dest = dest
        self._stab = stab
        self._is_observable = is_observable
        self._observable_index = observable_index

    @property
    def round(self) -> int:
        """
        Get the round number of the stabilizer propagation.

        Returns:
            int: The round number.
        """
        return self._round

    def is_observable(self) -> bool:
        """
        Check if the stabilizer propagation is for an observable.

        Returns:
            bool: True if it is an observable, False otherwise.
        """
        return self._is_observable

    def get_observable_index(self) -> int:
        """
        Get the index of the observable if applicable.

        Returns:
            int: The observable index.
        """
        return self._observable_index

    def get_stabindex(self) -> int:
        """
        Get the stabilizer index of the stabilizer propagation.

        Returns:
            int: The stabilizer index.
        """
        return self._stabindex

    def __str__(self) -> str:
        if self._is_observable:
            return f"{self._dest} = Prop {self._stab}"
        else:
            return f"{self._dest} = Prop[r={self._round}, s={self._stabindex}] {self._stab}"

    @property
    def dest(self) -> str:
        """
        Get the destination qubit/observable/detector from the instruction.

        Returns:
            str: The destination string.
        """
        return self._dest

    @property
    def stab(self) -> str:
        """
        Get the stabilizer from the instruction.

        Returns:
            str: The stabilizer string.
        """
        return self._stab


class ParityInstruction(IRInstruction):
    """Base class for parity-check IR instructions.

    A parity instruction computes the XOR (parity) of a set of
    measurement results identified by their symbolic names.  It is the
    common base for :class:`DetectorInstruction` and
    :class:`ObservableInstruction`.

    Args:
        ir_type: The IR instruction type (``DETECTOR`` or ``OBSERVABLE``).
        dest: Symbolic name for the parity result (e.g., ``"d0"`` or
            ``"o0"``).
        args: List of symbolic measurement-result names whose parity is
            computed.
    """

    def __init__(self, ir_type: IRType, dest: str, args: list[str]) -> None:
        super().__init__(ir_type)
        self._dest = dest
        self._args = args

    def __str__(self) -> str:
        return f"{self._dest} = Parity {' '.join(self._args)}"

    @property
    def dest(self) -> str:
        return self._dest

    @property
    def args(self) -> list[str]:
        return self._args


class DetectorInstruction(ParityInstruction):
    """IR instruction declaring a detector.

    A detector compares stabilizer measurement results across consecutive
    syndrome-extraction rounds.  In the absence of errors the parity of
    the referenced measurements is deterministic (always 0).

    Args:
        dest: Symbolic detector name (e.g., ``"d0"``).
        args: Symbolic names of the measurement results whose parity
            defines this detector.
    """

    def __init__(self, dest: str, args: list[str]) -> None:
        super().__init__(IRType.DETECTOR, dest, args)


class ObservableInstruction(ParityInstruction):
    """IR instruction declaring a logical observable.

    A logical observable ties measurement results to the logical-level
    outcome that the decoder must predict.  Typically it references the
    measurement of a logical Z (or X) operator.

    Args:
        dest: Symbolic observable name (e.g., ``"o0"``).
        args: Symbolic names of the measurement results whose parity
            defines this observable.
    """

    def __init__(self, dest: str, args: list[str]) -> None:
        super().__init__(IRType.OBSERVABLE, dest, args)


class IF_THENInstruction(IRInstruction):
    """IR instruction for conditional (if-then) control flow.

    .. note:: Not yet implemented; reserved for future use.

    Args:
        condition: Boolean expression string evaluated at runtime.
        then_branch: List of IR instructions to execute when
            *condition* is true.
    """

    def __init__(self, condition: str, then_branch: list[IRInstruction]) -> None:
        super().__init__(IRType.IF_THEN)
        self._condition = condition
        self._then_branch = then_branch


class WHILEInstruction(IRInstruction):
    """IR instruction for while-loop control flow.

    .. note:: Not yet implemented; reserved for future use.

    Args:
        condition: Boolean expression string re-evaluated each iteration.
        body: List of IR instructions forming the loop body.
    """

    def __init__(self, condition: str, body: list[IRInstruction]) -> None:
        super().__init__(IRType.WHILE)
        self._condition = condition
        self._body = body


class REPEAT_UNTILInstruction(IRInstruction):
    """IR instruction for repeat-until control flow.

    .. note:: Not yet implemented; reserved for future use.

    Args:
        body: List of IR instructions executed at least once per
            iteration.
        until_condition: Boolean expression string; the loop exits when
            this evaluates to true.
    """

    def __init__(self, body: list[IRInstruction], until_condition: str) -> None:
        super().__init__(IRType.REPEAT_UNTIL)
        self._body = body
        self._until_condition = until_condition


class REPEATInstruction(IRInstruction):
    """IR instruction for fixed-count repetition.

    .. note:: Not yet implemented; reserved for future use.

    Args:
        body: List of IR instructions forming the loop body.
        times: Number of times the body is repeated.
    """

    def __init__(self, body: list[IRInstruction], times: int) -> None:
        super().__init__(IRType.REPEAT)
        self._body = body
        self._times = times


class StabCode:
    """Stabilizer quantum error-correcting code.

    ``StabCode`` is the primary entry point for defining a stabilizer code,
    attaching a noise model, and compiling the resulting syndrome-extraction
    circuit into a ``stim.Circuit``.

    The compilation is a two-stage pipeline:

    1. **IR construction** -- call :meth:`construct_IR_standard_scheme` (or
       the scheme-specific variant) to build an intermediate representation
       from the stabilizer generators and logical operators.
    2. **Circuit compilation** -- call
       :meth:`compile_stim_circuit_from_IR_standard` to lower the IR into a
       :class:`~scalerqec.Clifford.clifford.CliffordCircuit` and a
       ``stim.Circuit``.

    Both stages are executed automatically by :meth:`construct_circuit`.

    Example::

        code = StabCode(n=5, k=1, d=3)
        code.add_stab("XZZXI")
        code.add_stab("IXZZX")
        code.add_stab("XIXZZ")
        code.add_stab("ZXIXZ")
        code.set_logical_Z(0, "ZZZZZ")
        code.noisemodel = NoiseModel(error_rate=0.001)
        code.scheme = "Standard"
        code.rounds = 3
        code.construct_circuit()
        print(code.stimcirc)

    Args:
        n: Number of physical (data) qubits.
        k: Number of logical qubits encoded by the code.
        d: Code distance (minimum weight of a non-trivial logical
            operator).
    """

    def __init__(self, n: int, k: int, d: int) -> None:
        self._n: int = n
        self._k: int = k
        self._d: int = d
        self._stabs: list[str] = []
        self._scheme: SCHEME = SCHEME.STANDARD
        self._circuit: CliffordCircuit | None = None
        self._stimcirc: stim.Circuit | None = None
        self._IRList: list[IRInstruction] = []
        self._rounds: int = 3 * d
        # Define the k different logical Z operators
        self._logicalZ: dict[int, str] = {}
        self._paritymatrix: Optional[NDArray[np.int_]] = None
        self._noisemodel: Optional[NoiseModel] = None
        self._IR_compiled = False
        self._circuit_compiled = False

    def is_IR_compiled(self) -> bool:
        """
        Check if the intermediate representation (IR) has been compiled.

        Returns:
            bool: True if the IR is compiled, False otherwise.
        """
        return self._IR_compiled

    def is_circuit_compiled(self) -> bool:
        """
        Check if the quantum error-correcting circuit has been compiled.

        Returns:
            bool: True if the circuit is compiled, False otherwise.
        """
        return self._circuit_compiled

    @property
    def n(self) -> int:
        """
        Get the number of physical qubits in the QECC.

        Returns:
            int: The number of physical qubits.
        """
        return self._n

    @property
    def k(self) -> int:
        """
        Get the number of logical qubits in the QECC.

        Returns:
            int: The number of logical qubits.
        """
        return self._k

    @property
    def d(self) -> int:
        """
        Get the distance of the QECC.

        Returns:
            int: The distance.
        """
        return self._d

    @property
    def noisemodel(self) -> Optional[NoiseModel]:
        """
        Get the noise model associated with the QECC.

        Returns:
            NoiseModel: The noise model.
        """
        return self._noisemodel

    @noisemodel.setter
    def noisemodel(self, noisemodel: NoiseModel) -> None:
        """
        Set the noise model associated with the QECC.

        Args:
            noisemodel (NoiseModel): The noise model to set.
        """
        self._noisemodel = noisemodel

    def init_by_parity_check_matrix(self, paritymatrix: NDArray[np.int_]) -> None:
        """Initialize the code from a binary parity-check matrix.

        Sets the number of physical qubits ``n`` and logical qubits ``k``
        from the matrix dimensions and replaces the current stabilizer
        list.  The stabilizer generators themselves are not yet extracted
        from the matrix (to be implemented).

        Args:
            paritymatrix: Binary parity-check matrix of shape
                ``(n - k, n)`` where rows correspond to stabilizer
                generators and columns to physical qubits.
        """
        self._paritymatrix = paritymatrix
        self._n = paritymatrix.shape[1]
        self._k = self._n - paritymatrix.shape[0]
        self._stabs = []
        pass

    def construct_parity_check_matrix(self) -> None:
        """Construct the binary symplectic parity-check matrix from stabilizers.

        Converts the stored Pauli-string stabilizer generators into a
        binary matrix in the standard ``[X | Z]`` symplectic format and
        stores it internally.

        .. note:: Not yet implemented.
        """
        pass

    def get_parity_check_matrix(self) -> Optional[NDArray[np.int_]]:
        """Return the binary symplectic parity-check matrix, if available.

        Returns:
            The parity-check matrix as a NumPy integer array of shape
            ``(n - k, n)``, or ``None`` if it has not been set or
            constructed.
        """
        return self._paritymatrix

    @property
    def circuit(self) -> Optional[CliffordCircuit]:
        """
        Get the Clifford circuit for the quantum error-correcting code.

        Returns:
            The Clifford circuit, or None if not yet compiled.
        """
        return self._circuit

    @property
    def stimcirc(self):
        """
        Get the STIM circuit for the quantum error-correcting code.

        Returns:
            The STIM circuit, or None if not yet compiled.
        """
        return self._stimcirc

    def set_logical_Z(self, index: int, logicalZ: str) -> None:
        """Set the logical Z operator for a logical qubit.

        The logical Z operator must be specified for every logical qubit
        (indices ``0`` through ``k - 1``) before calling
        :meth:`construct_circuit`.

        Args:
            index: Zero-based index of the logical qubit.
            logicalZ: Pauli string of length ``n`` (characters from
                ``I``, ``X``, ``Y``, ``Z``) representing the logical Z
                operator.

        Raises:
            AssertionError: If *logicalZ* has incorrect length or contains
                invalid characters.
        """
        assert len(logicalZ) == self._n, "Logical Z length must match number of qubits."
        assert all(c in "IXYZ" for c in logicalZ), (
            "Logical Z must only contain I, X, Y, and Z."
        )

        self._logicalZ[index] = logicalZ

    @property
    def rounds(self) -> int:
        """
        Get the number of error correction rounds.

        Returns:
            int: The number of rounds.
        """
        return self._rounds

    @rounds.setter
    def rounds(self, rounds: int) -> None:
        """
        Set the number of error correction rounds.

        Args:
            rounds (int): The number of rounds to set.
        """
        self._rounds = rounds

    def add_stab(self, stab: str) -> None:
        """Add a stabilizer generator to the code.

        Args:
            stab: Pauli string of length ``n`` (characters from ``I``,
                ``X``, ``Y``, ``Z``) representing the stabilizer
                generator.

        Raises:
            AssertionError: If *stab* has incorrect length or contains
                invalid characters.
        """
        assert len(stab) == self._n, "Stabilizer length must match number of qubits."
        assert all(c in "IXYZ" for c in stab), (
            "Stabilizer must only contain I, X, Y, Z."
        )

        self._stabs.append(stab)

    @property
    def scheme(self) -> SCHEME:
        """
        Get the error correction scheme for the code.

        Returns:
            SCHEME: The error correction scheme.
        """
        return self._scheme

    @property
    def stabilizers(self) -> list[str]:
        """
        Get the list of stabilizer generators for this code.

        Returns:
            list[str]: The stabilizer generators as Pauli strings.
        """
        return self._stabs

    @scheme.setter  # type: ignore[no-redef, attr-defined]
    def scheme(self, scheme: str) -> None:
        """
        Set the error correction scheme for the code.

        Args:
            scheme (SCHEME): The error correction scheme to use.
        """
        if scheme == "Standard":
            self._scheme = SCHEME.STANDARD
        elif scheme == "Shor":
            self._scheme = SCHEME.SHOR
        elif scheme == "Knill":
            self._scheme = SCHEME.KNILL
        elif scheme == "Flag":
            self._scheme = SCHEME.FLAG
        else:
            raise ValueError(f"Unknown scheme: {scheme}")

    def construct_circuit(self) -> None:
        """Run the full compilation pipeline for the selected scheme.

        Performs two stages:

        1. **IR construction** -- translates stabilizer generators and
           logical operators into an intermediate representation.
        2. **Circuit compilation** -- lowers the IR into a
           :class:`~scalerqec.Clifford.clifford.CliffordCircuit` and a
           ``stim.Circuit``.

        If a :attr:`noisemodel` has been attached, depolarizing noise is
        injected into the compiled circuit after stage 2.

        The IR uses a simple text-like notation internally::

            c0 = Prop XYZIX
            c1 = Prop IXYZI
            d0 = Parity c0 c1
            o0 = Parity c0

        Raises:
            NotImplementedError: If a scheme other than ``STANDARD`` is
                selected.
            ValueError: If any logical Z operator required by the code
                has not been set.
        """
        if self._scheme == SCHEME.STANDARD:
            self.construct_IR_standard_scheme()
            self.compile_stim_circuit_from_IR_standard()
            if self._noisemodel is not None:
                self._circuit = self._noisemodel.reconstruct_clifford_circuit(
                    self._circuit
                )
                self._stimcirc = self._circuit._stimcircuit
        else:
            raise NotImplementedError(f"Scheme {self._scheme} not implemented yet.")

    def construct_IR_shor_scheme(self) -> None:
        """Build the IR for Shor-style fault-tolerant syndrome extraction.

        .. note:: Not yet implemented.
        """
        pass

    def compile_stim_circuit_from_shor_standard(self) -> Optional[str]:
        """Compile a STIM circuit from the Shor-scheme IR.

        Returns:
            The compiled STIM circuit as a string, or ``None``.

        .. note:: Not yet implemented.
        """
        pass

    def construct_IR_knill_scheme(self) -> None:
        """Build the IR for Knill-style fault-tolerant syndrome extraction.

        .. note:: Not yet implemented.
        """
        pass

    def compile_stim_circuit_from_knill(self) -> Optional[str]:
        """Compile a STIM circuit from the Knill-scheme IR.

        Returns:
            The compiled STIM circuit as a string, or ``None``.

        .. note:: Not yet implemented.
        """
        pass

    def construct_IR_standard_scheme(self) -> None:
        """Build the IR for standard (bare-ancilla) syndrome extraction.

        For each of :attr:`rounds` syndrome-extraction rounds, a
        :class:`StabPropInstruction` is emitted for every stabilizer
        generator.  Starting from the second round, a
        :class:`DetectorInstruction` is added comparing the current
        round's measurement to the previous round's measurement for the
        same stabilizer.  Finally, :class:`ObservableInstruction` entries
        are emitted for each logical Z operator.

        This method is idempotent: calling it more than once has no
        effect.

        Raises:
            ValueError: If a logical Z operator has not been set for
                every logical qubit index ``0 .. k-1``.
        """
        if self._IR_compiled:
            return
        current_measurement_idx = 0
        current_detector_idx = 0
        prev_stab_meas_addr: dict[str, int] = {}
        for r in range(self._rounds):
            stabidx = 0
            for stab in self._stabs:
                dest = f"c{current_measurement_idx}"
                instr = StabPropInstruction(r, stabidx, dest, stab)
                self._IRList.append(instr)
                current_measurement_idx += 1
                stabidx += 1
                """
                Since the second round, add detectors comparing with previous round
                """
                if r > 0:
                    prev_dest = prev_stab_meas_addr[stab]
                    detector_dest = f"d{current_detector_idx}"
                    detector_instr = DetectorInstruction(
                        detector_dest, [prev_dest, dest]
                    )
                    self._IRList.append(detector_instr)
                    current_detector_idx += 1
                prev_stab_meas_addr[stab] = dest
        # Logical observables
        for logical_idx in range(self._k):
            if logical_idx not in self._logicalZ:
                raise ValueError(
                    f"Logical Z operator for logical qubit {logical_idx} not defined. "
                    f"Use set_logical_Z({logical_idx}, '<pauli_string>') before constructing the circuit."
                )
            logicalZ = self._logicalZ[logical_idx]

            dest = f"c{current_measurement_idx}"
            instr = StabPropInstruction(
                0, 0, dest, logicalZ, is_observable=True, observable_index=logical_idx
            )

            self._IRList.append(instr)
            current_measurement_idx += 1
            observable_dest = f"o{logical_idx}"
            observable_instr = ObservableInstruction(observable_dest, [dest])
            self._IRList.append(observable_instr)
        self._IR_compiled = True

    def show_IR(self) -> None:
        """Print the intermediate representation to stdout.

        Each IR instruction is printed on its own line using its
        ``__str__`` representation.  Useful for debugging the compilation
        pipeline.
        """
        for irinst in self._IRList:
            print(irinst)

    def compile_stim_circuit_from_IR_standard(self) -> str | None:
        """Lower the standard-scheme IR into a STIM circuit.

        Iterates over the IR instruction list and emits physical-level
        gates into a
        :class:`~scalerqec.Clifford.clifford.CliffordCircuit`:

        * **StabPropInstruction** -- resets the ancilla, applies
          controlled-Pauli entangling gates for each non-identity Pauli
          in the stabilizer, and measures the ancilla.
        * **DetectorInstruction** -- records a parity check across
          measurement results.
        * **ObservableInstruction** -- records a logical observable.

        Qubit layout convention:

        * Data qubits: indices ``0 .. n-1``.
        * Syndrome ancillas: indices ``n .. n + num_stabilizers - 1``.
        * Observable ancillas: indices
          ``n + num_stabilizers .. n + num_stabilizers + k - 1``.

        Returns:
            The compiled STIM circuit as a string if the circuit was
            already compiled, or ``None`` on first compilation (the
            result is stored in :attr:`stimcirc`).

        Raises:
            RuntimeError: If the IR has not been compiled yet.
        """
        # Convension: Stabilizer k stored in qubit n+k-1
        # Observable k stored in qubit n+num_syndromes+k-1

        if not self._IR_compiled:
            raise RuntimeError("IR not compiled yet.")
        if self._circuit_compiled:
            return str(self._stimcirc)
        self._circuit = CliffordCircuit(self._n + len(self._stabs) + self._k)
        parity_match_group = []
        observable_parity_group = []

        dest_to_measure_index = {}
        current_measure_index = 0
        for irinst in self._IRList:
            if isinstance(irinst, StabPropInstruction):
                stab = irinst.stab
                dest_index = int(irinst.dest[1:])
                if irinst.is_observable():
                    helper_qubit_index = (
                        self._n + len(self._stabs) + irinst.get_observable_index()
                    )
                else:
                    helper_qubit_index = self._n + irinst.get_stabindex()

                self._circuit.add_reset(helper_qubit_index)
                for qubit_index, pauli in enumerate(stab):
                    if pauli == "X":
                        self._circuit.add_hadamard(qubit_index)
                        self._circuit.add_cnot(
                            control=qubit_index, target=helper_qubit_index
                        )
                        self._circuit.add_hadamard(qubit_index)
                    elif pauli == "Z":
                        self._circuit.add_cnot(
                            control=qubit_index, target=helper_qubit_index
                        )
                    elif pauli == "I":
                        continue
                    elif pauli == "Y":
                        # Y = iXZ, so we need to apply both X and Z parity propagation
                        # First apply X part: H-CNOT-H
                        self._circuit.add_hadamard(qubit_index)
                        self._circuit.add_cnot(
                            control=qubit_index, target=helper_qubit_index
                        )
                        self._circuit.add_hadamard(qubit_index)
                        # Then apply Z part: CNOT
                        self._circuit.add_cnot(
                            control=qubit_index, target=helper_qubit_index
                        )

                self._circuit.add_measurement(helper_qubit_index)
                dest_to_measure_index[irinst.dest] = current_measure_index
                current_measure_index += 1

            elif isinstance(irinst, DetectorInstruction):
                args = irinst.args
                args_measure_indices = [dest_to_measure_index[arg] for arg in args]
                parity_match_group.append(args_measure_indices)

            elif isinstance(irinst, ObservableInstruction):
                args = irinst.args
                args_indices = [dest_to_measure_index[arg] for arg in args]
                observable_parity_group.append(args_indices)

        self._circuit.parityMatchGroup = parity_match_group
        self._circuit.observable = observable_parity_group[0]
        self._circuit.compile_detector_and_observable()
        self._stimcirc = self._circuit._stimcircuit
        self._circuit_compiled = True


def test_commute():
    assert commute("IXYZ", "IYZX") == False
    assert commute("XZZI", "ZXXI") == False
    assert commute("IIII", "ZZZZ") == True
    assert commute("XIZY", "YZXI") == True


if __name__ == "__main__":
    noise_model = NoiseModel(error_rate=0.001)
    qeccirc = StabCode(n=5, k=1, d=3)
    qeccirc.noisemodel = noise_model
    # Specify your stabilizers
    # Stabilizer generators
    qeccirc.add_stab("XZZXI")
    qeccirc.add_stab("IXZZX")
    qeccirc.add_stab("XIXZZ")
    qeccirc.add_stab("ZXIXZ")
    qeccirc.set_logical_Z(0, "ZZZZZ")
    # Set stabilizer parity measurement scheme, round of repetition
    qeccirc.scheme = "Standard"  # type: ignore[misc, assignment]
    qeccirc.rounds = 3
    qeccirc.construct_circuit()
    stim_circuit = qeccirc.stimcirc
    print(stim_circuit)
