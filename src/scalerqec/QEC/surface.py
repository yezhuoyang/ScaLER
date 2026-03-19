"""Standard quantum error correction code families.

Provides parametric code classes that automatically generate stabilizer
generators and logical operators:

* :class:`RepetitionCode` -- distance-d Z-repetition code.
* :class:`SurfaceCode` -- distance-d rotated surface code.

Both inherit from :class:`~scalerqec.QEC.qeccircuit.StabCode` and can be
compiled to Stim circuits via ``construct_circuit()``.
"""

from __future__ import annotations

from .qeccircuit import StabCode


class RepetitionCode(StabCode):
    """Distance-d Z-repetition code.

    The repetition code encodes 1 logical qubit into *d* physical qubits
    using *d - 1* weight-2 ZZ stabilizers.  It can detect and correct
    up to ``(d - 1) / 2`` single-qubit Z (phase-flip) errors.

    Stabilizers: ``ZZ I...I``, ``I ZZ I...I``, ..., ``I...I ZZ``
    Logical Z: ``Z Z ... Z`` (all qubits)

    Args:
        distance: Code distance (number of physical qubits).
        rounds: Number of syndrome extraction rounds. Defaults to *d*.

    Example::

        code = RepetitionCode(distance=3)
        code.rounds = 3
        code.construct_circuit()
        print(code.stimcirc)
    """

    def __init__(self, distance: int, rounds: int | None = None) -> None:
        super().__init__(n=distance, k=1, d=distance)
        self._construct_stabilizers()
        self.set_logical_Z(0, "Z" * distance)
        if rounds is not None:
            self._rounds = rounds
        else:
            self._rounds = distance

    def _construct_stabilizers(self) -> None:
        for i in range(self._n - 1):
            stab = "I" * i + "ZZ" + "I" * (self._n - i - 2)
            self.add_stab(stab)


class SurfaceCode(StabCode):
    """Distance-d rotated surface code.

    The rotated surface code encodes 1 logical qubit into ``d * d``
    physical qubits.  Data qubits sit on a d x d grid; X-type and
    Z-type stabilizers occupy alternating plaquettes.

    Data qubit at grid position ``(row, col)`` has index ``row * d + col``.

    * Z stabilizers: plaquettes with top-left corner at even row+col
    * X stabilizers: plaquettes with top-left corner at odd row+col
    * Logical Z: leftmost column of Z operators
    * Logical X: top row of X operators

    Args:
        distance: Code distance (odd integer >= 3).
        rounds: Number of syndrome extraction rounds. Defaults to *d*.

    Example::

        code = SurfaceCode(distance=3)
        code.construct_circuit()
        print(code.stimcirc)
    """

    def __init__(self, distance: int, rounds: int | None = None) -> None:
        n = distance * distance
        super().__init__(n=n, k=1, d=distance)
        self._distance = distance

        # Assign qubit coordinates
        self._qubit_coords: dict[int, tuple[float, float]] = {}
        for row in range(distance):
            for col in range(distance):
                q = row * distance + col
                self._qubit_coords[q] = (2.0 * col, 2.0 * row)

        self._construct_stabilizers()
        self._set_logical_operators()

        if rounds is not None:
            self._rounds = rounds
        else:
            self._rounds = distance

    def _qubit_index(self, row: int, col: int) -> int:
        return row * self._distance + col

    def _construct_stabilizers(self) -> None:
        """Generate X-type and Z-type stabilizers for the rotated surface code.

        We use the standard rotated layout where plaquettes tile the
        ``(d-1) x (d-1)`` grid of face centers.  Each face at position
        ``(fr, fc)`` (0-indexed) touches up to 4 data qubits:
        ``(fr, fc), (fr, fc+1), (fr+1, fc), (fr+1, fc+1)``.

        The stabilizer type alternates in a checkerboard pattern:
        - ``(fr + fc) % 2 == 0`` → Z-type stabilizer
        - ``(fr + fc) % 2 == 1`` → X-type stabilizer

        Boundary stabilizers (weight-2) come from virtual faces beyond
        the grid:
        - Top boundary (fr=-1): Z-type when fc is odd
        - Bottom boundary (fr=d-1): Z-type when fc is even
        - Left boundary (fc=-1): X-type when fr is even
        - Right boundary (fc=d-1): X-type when fr is odd
        """
        d = self._distance

        # Bulk stabilizers (weight-4)
        for fr in range(d - 1):
            for fc in range(d - 1):
                qubits = [
                    self._qubit_index(fr, fc),
                    self._qubit_index(fr, fc + 1),
                    self._qubit_index(fr + 1, fc),
                    self._qubit_index(fr + 1, fc + 1),
                ]
                pauli_type = "Z" if (fr + fc) % 2 == 0 else "X"
                stab = ["I"] * self._n
                for q in qubits:
                    stab[q] = pauli_type
                self.add_stab("".join(stab))

        # Boundary stabilizers (weight-2)
        # Top boundary (virtual face at fr=-1):
        # Type is Z when (-1 + fc) % 2 == 0 → fc is odd
        for fc in range(d - 1):
            if fc % 2 == 1:  # fc odd → Z-type
                qubits = [self._qubit_index(0, fc), self._qubit_index(0, fc + 1)]
                stab = ["I"] * self._n
                for q in qubits:
                    stab[q] = "Z"
                self.add_stab("".join(stab))

        # Bottom boundary (virtual face at fr=d-1):
        # Type is Z when (d-1 + fc) % 2 == 0 → fc even (for odd d)
        for fc in range(d - 1):
            if (d - 1 + fc) % 2 == 0:  # Z-type
                qubits = [
                    self._qubit_index(d - 1, fc),
                    self._qubit_index(d - 1, fc + 1),
                ]
                stab = ["I"] * self._n
                for q in qubits:
                    stab[q] = "Z"
                self.add_stab("".join(stab))

        # Left boundary (virtual face at fc=-1):
        # Type is X when (fr + (-1)) % 2 == 1 → fr is even
        for fr in range(d - 1):
            if fr % 2 == 0:  # fr even → X-type
                qubits = [
                    self._qubit_index(fr, 0),
                    self._qubit_index(fr + 1, 0),
                ]
                stab = ["I"] * self._n
                for q in qubits:
                    stab[q] = "X"
                self.add_stab("".join(stab))

        # Right boundary (virtual face at fc=d-1):
        # Type is X when (fr + d-1) % 2 == 1 → fr is odd (for odd d)
        for fr in range(d - 1):
            if (fr + d - 1) % 2 == 1:  # X-type
                qubits = [
                    self._qubit_index(fr, d - 1),
                    self._qubit_index(fr + 1, d - 1),
                ]
                stab = ["I"] * self._n
                for q in qubits:
                    stab[q] = "X"
                self.add_stab("".join(stab))

    def _set_logical_operators(self) -> None:
        """Set logical Z = leftmost column, logical X = top row."""
        d = self._distance

        # Logical Z: Z on all qubits in column 0
        logZ = ["I"] * self._n
        for row in range(d):
            logZ[self._qubit_index(row, 0)] = "Z"
        self.set_logical_Z(0, "".join(logZ))
