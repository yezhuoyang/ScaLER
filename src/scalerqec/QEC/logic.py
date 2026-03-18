"""Logical-level circuit abstraction and parser.

This module provides a high-level representation of logical quantum
circuits operating on QEC code blocks.  It is designed for expressing
protocols such as magic state distillation that are naturally described
at the logical level rather than the physical level.

Key components:

* :class:`CodeBlock` -- a named QEC code block parameterized by
  ``(n, k, d)`` that holds logical qubits.
* Logical gate classes (:class:`LogicalH`, :class:`LogicalCNOT`,
  :class:`LogicalT`, :class:`InjectT`, :class:`LogicalMeasure`,
  :class:`LogicalReset`) -- operations on logical qubits within code
  blocks.
* :class:`LogicalCircuit` -- ordered container of code blocks and
  logical gates.
* :class:`LogicalParser` -- text parser that converts a human-readable
  circuit description string into a :class:`LogicalCircuit`.
"""

from __future__ import annotations

import re
from .qeccircuit import StabCode


"""

Synatatic design:

First, user has to allocate QEC code blocks to hold logical qubits.
All operations are defined on the allocated blocks


Type surface:
    The stabilizers should be a function of d
    This type describe the stabilizer structure.


surface q1 [n1,k1,d1]
surface q2 [n2,k2,d2]

surface t0 [n3,k3,d3]   # magic T state block

q1[0] = LogicH q1[0]

t0 = Distill15to1_T[d=25]     # returns a magic_T handle
InjectT q1[0], t0

q2[1] = LogicCNOT q1[0], q2[1]

c1 = LogicMeasure q1[0]
c2 = LogicMeasure q2[1]
"""


"""
A compiler translate LogicQ to a slightly lower-level IR which is more specific how to implement logical operations on specific QEC codes.


TRANSVERSAL H q1[0]            # Transversal H on logical qubit 0 of block q1
QECCycle q1                    # Perform one QEC cycle on block q1, including stabilizer measurements and corrections
Transversal CNOT q1[0], q2[1]  # Transversal CNOT between logical qubit 0 of block q1 and logical qubit 1 of block q2
Transversal H q2[1]
XXZZXIXII   q1[0]                # Also support direct Pauli operations, if needed

QECCycle q1                    # Perform one QEC cycle on block q1, including stabilizer measurements and corrections
QECCycle q2                    # Perform one QEC cycle on block q2, including stabilizer measurements and corrections

c1 = Prop LogicZ q1[0]          # Propagate logical Z operator to classical bit c1
c2 = Prop LogicZ q2[1]          # Propagate logical Z operator to classical bit c2


"""


"""
Magic state distillation protocol:

Now we have a language which works on the logical level. We can define the MGD protocol as follows:
15-to-1 Reed-Muller magic state distillation (logical-level).

Precondition:
  - surface f has k >= 15 logical qubits.
  - On each attempt, f[0..14] are initialized to 15 noisy |T> states
    (supplied externally; this IR does not prepare |T> states).

Postcondition:
  - returns a distilled magic_T handle backed by the encoded logical qubit.


  protocol Distill15to1_T(surface f, int d) -> magic_T:
    Repeat:

        # ---- X-type stabilizer checks ----
        c_x1 = LogicProp IIIIIIIXXXXXXXX
        c_x2 = LogicProp IIIXXXXIIIIXXXX
        c_x3 = LogicProp IXXIIXXIIXXIIXX
        c_x4 = LogicProp XIXIXIXIXIXIXIX

        # ---- Z-type stabilizer checks ----
        c_z1  = LogicProp IIIIIIIIZZZZZZZZ
        c_z2  = LogicProp IIIZZZZIIIIZZZZ
        c_z3  = LogicProp IZZIIZZIIZZIIZZ
        c_z4  = LogicProp ZIZIZIZIZIZIZIZ
        c_z12 = LogicProp IIIIIIIIIIZZZZ
        c_z13 = LogicProp IIIIIIIIZZIIIZZ
        c_z14 = LogicProp IIIIIIIIZIZIZIZ
        c_z23 = LogicProp IIIIIZZIIIIIIZZ
        c_z24 = LogicProp IIIIZIZIIIIIZIZ
        c_z34 = LogicProp IIZIIIZIIIZIIIZ

        Success = c_x1 == 0 && c_x2 == 0 && c_x3 == 0 && c_x4 == 0 &&
                  c_z1 == 0 && c_z2 == 0 && c_z3 == 0 && c_z4 == 0 &&
                  c_z12 == 0 && c_z13 == 0 && c_z14 == 0 &&
                  c_z23 == 0 && c_z24 == 0 && c_z34 == 0
        Until Success

        return 


"""


class CodeBlock:
    """A named QEC code block that holds logical qubits.

    Represents a single instance of a quantum error-correcting code
    (e.g., a surface code patch) identified by a user-chosen name and
    parameterized by ``(n, k, d)``.

    Args:
        type: Code family name (e.g., ``"surface"``, ``"color"``,
            ``"LDPC"``).
        name: Unique identifier for this block within a
            :class:`LogicalCircuit`.
        n: Number of physical qubits.
        k: Number of logical qubits encoded in this block.
        d: Code distance.
    """

    def __init__(self, type: str, name: str, n: int, k: int, d: int):
        self._name = name  # Name of the code block
        self._n = n  # Number of physical qubits
        self._k = k  # Number of logical qubits
        self._d = d  # Code distance
        self._type = type  # Type of the code, e.g., 'surface', 'color', 'LDPC'

    def __repr__(self):
        return f"{self._name}[{self._n},{self._k},{self._d}]"


class LogicalGate:
    """Abstract base class for logical-level gate operations.

    Subclasses must implement :meth:`__repr__` to provide a human-readable
    representation of the gate.

    Args:
        type: Gate type identifier (e.g., ``"H"``, ``"CNOT"``, ``"T"``).
    """

    def __init__(self, type: str):
        self._type = type  # Type of the logical gate, e.g., 'CNOT', 'H', 'T'

    @property
    def type(self) -> str:
        """str: The gate type identifier."""
        return self._type

    def __repr__(self):
        raise NotImplementedError("Subclasses must implement __repr__ method.")


class LogicalH(LogicalGate):
    """Logical Hadamard gate applied to a single logical qubit.

    Args:
        block: The code block containing the target logical qubit.
        index: Zero-based index of the logical qubit within *block*.
    """

    def __init__(self, block: CodeBlock, index: int):
        super().__init__("H")
        self._block = block
        self._index = index

    def __repr__(self):
        return f"LogicalH {self._block._name}[{self._index}]"


class LogicalCNOT(LogicalGate):
    """Logical CNOT gate between two logical qubits.

    The control and target qubits may reside in different code blocks.

    Args:
        control_block: Code block containing the control logical qubit.
        control_index: Zero-based index of the control logical qubit.
        target_block: Code block containing the target logical qubit.
        target_index: Zero-based index of the target logical qubit.
    """

    def __init__(
        self,
        control_block: CodeBlock,
        control_index: int,
        target_block: CodeBlock,
        target_index: int,
    ):
        super().__init__("CNOT")
        self._control_block = control_block
        self._target_block = target_block
        self._control_index = control_index
        self._target_index = target_index

    def __repr__(self):
        return f"LogicalCNOT {self._control_block._name}[{self._control_index}], {self._target_block._name}[{self._target_index}]"


class LogicalT(LogicalGate):
    """Logical T gate applied to a single logical qubit.

    The T gate is a non-Clifford gate typically realized via magic state
    injection at the logical level.

    Args:
        block: The code block containing the target logical qubit.
        index: Zero-based index of the logical qubit within *block*.
    """

    def __init__(self, block: CodeBlock, index: int):
        super().__init__("T")
        self._block = block
        self._index = index

    def __repr__(self):
        return f"LogicalT {self._block._name}[{self._index}]"


class InjectT(LogicalGate):
    """Inject a distilled magic T state into a logical qubit.

    Consumes a previously prepared magic-T handle and teleports the
    encoded T state into the destination logical qubit.

    Args:
        dest_block: Code block containing the destination logical qubit.
        dest_index: Zero-based index of the destination logical qubit.
        magic_T_handle: Identifier of the magic-T resource produced by
            a distillation protocol.
    """

    def __init__(self, dest_block: CodeBlock, dest_index: int, magic_T_handle: str):
        super().__init__("InjectT")
        self._dest_block = dest_block
        self._dest_index = dest_index
        self._magic_T_handle = magic_T_handle

    def __repr__(self):
        return f"InjectT {self._dest_block._name}[{self._dest_index}], {self._magic_T_handle}"


class LogicalMeasure(LogicalGate):
    """Destructive measurement of a logical qubit in the Z basis.

    The measurement outcome is stored in a classical register addressed
    by *cindex*.

    .. note::
        Currently only logical-Z measurement is supported.  More general
        logical measurements (e.g., MXX, MZX) are planned.

    Args:
        block: Code block containing the logical qubit to measure.
        cindex: Index of the classical bit that receives the outcome.
        index: Zero-based index of the logical qubit within *block*.
    """

    def __init__(self, block: CodeBlock, cindex: int, index: int):
        super().__init__("Measure")
        self._block = block
        self._cindex = cindex
        self._index = index

    @property
    def index(self) -> int:
        """int: Zero-based logical qubit index within the code block."""
        return self._index

    @property
    def block(self) -> CodeBlock:
        """CodeBlock: The code block containing the measured qubit."""
        return self._block

    @property
    def cindex(self) -> int:
        """int: Classical register index for the measurement outcome."""
        return self._cindex

    def __repr__(self):
        return f"c[{self._cindex}] = LogicalMeasure {self._block._name}[{self._index}]"


class LogicalReset(LogicalGate):
    """Reset a logical qubit to the logical |0> state.

    Args:
        block: Code block containing the logical qubit to reset.
        index: Zero-based index of the logical qubit within *block*.
    """

    def __init__(self, block: CodeBlock, index: int):
        super().__init__("Reset")
        self._block = block
        self._index = index

    @property
    def index(self) -> int:
        """int: Zero-based logical qubit index within the code block."""
        return self._index

    @property
    def block(self) -> CodeBlock:
        """CodeBlock: The code block containing the qubit to reset."""
        return self._block

    def __repr__(self):
        return f"LogicalReset {self._block._name}[{self._index}]"


class LogicalCircuit:
    """Ordered container of code blocks and logical gates.

    A ``LogicalCircuit`` collects :class:`CodeBlock` declarations and
    a sequence of :class:`LogicalGate` operations that act on logical
    qubits within those blocks.  Gate compatibility (qubit index within
    range) is validated on insertion.
    """

    def __init__(self):
        self.gates: list[LogicalGate] = []
        self._blocks: list[CodeBlock] = []
        self._qec_types: list[str] = []
        self._qec_type_names: list[str] = []
        self._MGT_handles: list[str] = []

    def add_qec_type(self, qec_type: str) -> None:
        """Register a QEC code family name (e.g., ``"surface"``).

        Args:
            qec_type: Code family identifier used during parsing.
        """
        self._qec_types.append(qec_type)
        self._qec_type_names.append(qec_type)

    def add_MGT_handle(self, handle: str) -> None:
        """Register a magic-T state handle for use with :class:`InjectT`.

        Args:
            handle: Unique identifier for the magic-T resource.
        """
        self._MGT_handles.append(handle)

    def add_block(self, block: CodeBlock) -> None:
        """Add a code block to the circuit.

        Args:
            block: The :class:`CodeBlock` instance to register.
        """
        self._blocks.append(block)

    def add_gate(self, gate: LogicalGate) -> None:
        """Validate and append a logical gate to the circuit.

        Checks that all qubit indices referenced by *gate* are within the
        valid range of their respective code blocks before appending.

        Args:
            gate: The :class:`LogicalGate` to add.

        Raises:
            ValueError: If a qubit index exceeds the number of logical
                qubits in its code block, or if the gate type is
                unsupported, or if an ``InjectT`` references an
                unregistered magic-T handle.
        """
        if gate.type in ("H", "T", "Measure", "Reset"):
            # Single-qubit gates
            if gate._index >= gate._block._k:  # type: ignore[attr-defined]
                raise ValueError(
                    f"Index {gate._index} out of range for block {gate._block._name} with {gate._block._k} logical qubits."
                )  # type: ignore[attr-defined]
            self.gates.append(gate)
        elif gate.type == "CNOT":
            # Two-qubit gates
            if gate._control_index >= gate._control_block._k:  # type: ignore[attr-defined]
                raise ValueError(
                    f"Control index {gate._control_index} out of range for block {gate._control_block._name} with {gate._control_block._k} logical qubits."
                )  # type: ignore[attr-defined]
            if gate._target_index >= gate._target_block._k:  # type: ignore[attr-defined]
                raise ValueError(
                    f"Target index {gate._target_index} out of range for block {gate._target_block._name} with {gate._target_block._k} logical qubits."
                )  # type: ignore[attr-defined]
            self.gates.append(gate)
        elif gate.type == "InjectT":
            if gate._dest_index >= gate._dest_block._k:  # type: ignore[attr-defined]
                raise ValueError(
                    f"Index {gate._dest_index} out of range for block {gate._dest_block._name} with {gate._dest_block._k} logical qubits."
                )  # type: ignore[attr-defined]
            if gate._MGT_handle not in self._MGT_handles:  # type: ignore[attr-defined]
                raise ValueError(f"Magic T handle {gate._MGT_handle} not recognized.")  # type: ignore[attr-defined]
            self.gates.append(gate)
        else:
            raise ValueError(f"Unsupported gate type: {gate.type}")

    def __repr__(self):
        output = ""
        for gate in self.gates:
            output += repr(gate) + "\n"
        return output


class LogicalParser:
    """Parser that converts a text description into a :class:`LogicalCircuit`.

    The input language supports three kinds of statements:

    * **Type declarations** -- ``Type surface``
    * **Block definitions** -- ``surface q1 [13,1,3]``
    * **Gate instructions** -- ``LogicalH q1[0]``,
      ``LogicalCNOT q1[0], q2[0]``, ``InjectT q1[0], t0``,
      ``c[0] = LogicalMeasure q1[0]``

    Lines may contain ``#``-style comments.  Blank lines are ignored.

    Example::

        Type surface
        surface q1 [13,1,3]
        surface q2 [13,1,3]
        LogicalH q1[0]
        LogicalCNOT q1[0], q2[0]
        c[1] = LogicalMeasure q1[0]
        c[2] = LogicalMeasure q2[0]
    """

    # Patterns defined at class level to avoid AttributeError
    _index_re = re.compile(r"(\w+)\[(\d+)\]")
    _triplet_re = re.compile(r"\[(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\]")
    _param_re = re.compile(r"\[(\w+)\s*=\s*(\d+)\]")  # For [d=25]

    def __init__(self):
        self._qec_types: list[str] = []
        self._qec_type_names: list[str] = []
        self._blocksmap: dict[str, CodeBlock] = {}

    def _parse_indexed(self, text: str) -> tuple[str, int]:
        """Extract a name and integer index from an indexed reference.

        Args:
            text: String of the form ``"name[index]"`` (e.g., ``"q1[0]"``).

        Returns:
            A ``(name, index)`` tuple.

        Raises:
            ValueError: If *text* does not match the expected pattern.
        """
        match = self._index_re.search(text)
        if not match:
            raise ValueError(f"Could not parse indexed reference: {text}")
        return match.group(1), int(match.group(2))

    def parse(self, logical_circ_string: str) -> LogicalCircuit:
        """Parse a circuit description string into a :class:`LogicalCircuit`.

        Performs three sequential passes over the input:

        1. **Type declarations** -- registers QEC code family names.
        2. **Block definitions** -- creates :class:`CodeBlock` instances.
        3. **Gate instructions** -- creates and validates logical gates.

        Args:
            logical_circ_string: Multi-line circuit description in the
                logical circuit language.

        Returns:
            A fully populated :class:`LogicalCircuit` instance.

        Raises:
            ValueError: If an indexed reference cannot be parsed, or if
                an unsupported gate operation is encountered.
        """
        circuit = LogicalCircuit()
        # Remove comments and empty lines
        lines = []
        for raw_line in logical_circ_string.strip().split("\n"):
            clean_line = raw_line.split("#")[0].strip()
            if clean_line:
                lines.append(clean_line)

        # Pass 1: Types
        for line in lines:
            parts = line.split()
            if parts[0] == "Type":
                circuit.add_qec_type(parts[1])

        # Pass 2: Block Definitions
        for line in lines:
            parts = line.split()
            if parts[0] in circuit._qec_type_names:
                qec_type = parts[0]
                block_name = parts[1]
                # Use regex to find [n,k,d]
                triplet_match = self._triplet_re.search(line)
                if triplet_match:
                    n, k, d = map(int, triplet_match.groups())
                    block = CodeBlock(qec_type, block_name, n, k, d)
                    self._blocksmap[block_name] = block
                    circuit.add_block(block)

        # Pass 3: Gates and Assignments
        for line in lines:
            parts = line.replace(
                ",", " "
            ).split()  # Replace comma with space for easy splitting
            if parts[0] in circuit._qec_type_names or parts[0] == "Type":
                continue

            # Case 1: Assignment (e.g., c[0] = LogicalMeasure q1[0])
            if "=" in line:
                lhs, rhs_full = line.split("=", 1)
                lhs_parts = lhs.strip().split()
                rhs_parts = rhs_full.strip().split()

                # Check for classical assignment
                if lhs_parts[0].startswith("c["):
                    c_name, c_idx = self._parse_indexed(lhs_parts[0])
                    gate_op = rhs_parts[0]

                    if gate_op == "LogicalMeasure":
                        b_name, b_idx = self._parse_indexed(rhs_parts[1])
                        circuit.add_gate(
                            LogicalMeasure(self._blocksmap[b_name], c_idx, b_idx)
                        )

                    else:
                        raise ValueError(
                            f"Unsupported gate operation in assignment: {gate_op}"
                        )

            # Case 2: In-place Gates (e.g., InjectT q1[0], t0)
            else:
                gate_type = parts[0]
                if gate_type == "LogicalH":
                    name, idx = self._parse_indexed(parts[1])
                    circuit.add_gate(LogicalH(self._blocksmap[name], idx))

                elif gate_type == "LogicalCNOT":
                    ctrl_name, ctrl_idx = self._parse_indexed(parts[1])
                    tgt_name, tgt_idx = self._parse_indexed(parts[2])
                    circuit.add_gate(
                        LogicalCNOT(
                            self._blocksmap[ctrl_name],
                            ctrl_idx,
                            self._blocksmap[tgt_name],
                            tgt_idx,
                        )
                    )

                elif gate_type == "InjectT":
                    dest_name, dest_idx = self._parse_indexed(parts[1])
                    handle = parts[2]
                    circuit.add_gate(
                        InjectT(self._blocksmap[dest_name], dest_idx, handle)
                    )

        return circuit


if __name__ == "__main__":
    parser = LogicalParser()

    logical_circ_str = """
    Type surface
    surface q1 [13,1,3]
    surface q2 [13,1,3]
    LogicalH q1[0]
    LogicalCNOT q1[0], q2[0]
    c[1] = LogicalMeasure q1[0]
    c[2] = LogicalMeasure q2[0]
    """

    logical_circuit = parser.parse(logical_circ_str)
    print(logical_circuit)
