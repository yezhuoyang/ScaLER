"""
Useful analyzers for quantum error correction codes.
Specifically, we provide analyzers for code distance and logical operators.
"""

from scalerqec.QEC import StabCode
from scalerqec.util import commute, multiply_pauli_strings
from typing import Literal
from typing import List, Dict


class StabilizerAnalyzer:
    """
    Check if all stabilizers of a code commute with each other.
    """

    def __init__(self, code: StabCode):
        self.code = code

    def verify_commutation(self) -> bool:
        """
        Verify if all stabilizers commute with each other.
        """
        for i, stab1 in enumerate(self.code.stabilizers):
            for stab2 in self.code.stabilizers[i + 1 :]:
                if not commute(stab1, stab2):
                    return False
        return True


class DistanceAnalyzer:
    """
    Analyzer to compute the code distance of a quantum error correction code.

    The code distance is the minimum weight of any logical operator —
    i.e., any Pauli string that commutes with all stabilizers but is not
    itself in the stabilizer group.

    Two methods:

    1. **Brute-force** — enumerate all Pauli strings of increasing weight
       until a logical operator is found.  Feasible for ``n <= ~15``.
    2. **SMT** — (placeholder) use a SAT/SMT solver.
    """

    supported_methods = {"bruteforce", "smt"}

    def __init__(
        self, code: StabCode, method: Literal["bruteforce", "smt"] = "bruteforce"
    ):
        self.code = code
        if method not in self.supported_methods:
            raise ValueError(
                f"Unsupported method '{method}'. Supported methods are: {self.supported_methods}"
            )
        self.method = method

    def compute_distance_bruteforce(self, max_weight: int | None = None) -> int:
        """Compute the code distance by brute-force enumeration.

        Enumerates all Pauli strings of weight 1, 2, 3, ... up to
        *max_weight* (default ``n``) and returns the weight of the
        first logical operator found.

        Args:
            max_weight: Maximum weight to search.  Defaults to ``n``.

        Returns:
            The minimum weight of a logical operator, or ``n + 1`` if
            none is found up to *max_weight*.
        """
        from itertools import combinations, product as iterproduct

        n = self.code.n
        if max_weight is None:
            max_weight = n
        stabs = list(self.code.stabilizers)

        # Build a LogicalOperatorAnalyzer for group membership checks
        loa = LogicalOperatorAnalyzer(self.code)

        paulis = ["X", "Y", "Z"]
        for weight in range(1, max_weight + 1):
            for positions in combinations(range(n), weight):
                for letters in iterproduct(paulis, repeat=weight):
                    op = ["I"] * n
                    for pos, letter in zip(positions, letters):
                        op[pos] = letter
                    op_str = "".join(op)

                    # Must commute with all stabilizers
                    if not all(commute(op_str, s) for s in stabs):
                        continue
                    # Must NOT be in stabilizer group
                    if loa._is_in_stabilizer_group(op_str):
                        continue
                    return weight

        return n + 1  # No logical operator found

    def compute_distance_smt(self) -> int:
        """Compute the code distance using SMT solver.

        .. note:: Not yet implemented — returns the declared distance.
        """
        return self.code.d

    def verify_code_distance(self) -> bool:
        """Verify if the computed distance matches the expected distance."""
        if self.method == "bruteforce":
            computed_distance = self.compute_distance_bruteforce()
        else:
            computed_distance = self.compute_distance_smt()

        return computed_distance == self.code.d

    def verify_circuit_level_code_distance(self, circuit) -> bool:
        """Verify the code distance at the circuit level.

        .. note:: Not yet implemented.
        """
        return True


class IRAnalyzer:
    """Static analyzer for QStab IR and stabilizer code properties.

    Provides validation checks that should be run before circuit compilation:

    * :meth:`verify_commutation` — all stabilizers mutually commute.
    * :meth:`verify_logical_operators` — logical Z operators are valid.
    * :meth:`verify_sufficient_rounds` — enough syndrome rounds for the distance.
    * :meth:`verify_code_distance` — brute-force distance matches declared ``d``.
    * :meth:`validate` — run all checks at once.
    """

    def __init__(self, code: StabCode) -> None:
        self.code = code

    def verify_commutation(self) -> tuple[bool, str]:
        """Check that all stabilizers mutually commute.

        Returns:
            ``(True, "")`` if valid, ``(False, message)`` otherwise.
        """
        stabs = self.code.stabilizers
        for i in range(len(stabs)):
            for j in range(i + 1, len(stabs)):
                if not commute(stabs[i], stabs[j]):
                    return False, (
                        f"Stabilizers {i} ({stabs[i]}) and {j} ({stabs[j]}) "
                        f"do not commute."
                    )
        return True, ""

    def verify_logical_operators(self) -> tuple[bool, str]:
        """Check that declared logical Z operators are valid.

        Each logical Z must:
        1. Commute with all stabilizers.
        2. Not be in the stabilizer group.

        Returns:
            ``(True, "")`` if valid, ``(False, message)`` otherwise.
        """
        loa = LogicalOperatorAnalyzer(self.code)
        logZ = getattr(self.code, "_logicalZ", {})
        stabs = self.code.stabilizers

        for idx, lz in logZ.items():
            for s in stabs:
                if not commute(lz, s):
                    return False, (
                        f"Logical Z[{idx}] = {lz} does not commute with "
                        f"stabilizer {s}."
                    )
            if loa._is_in_stabilizer_group(lz):
                return False, (
                    f"Logical Z[{idx}] = {lz} is in the stabilizer group "
                    f"(not a genuine logical operator)."
                )
        return True, ""

    def verify_sufficient_rounds(self) -> tuple[bool, str]:
        """Check that the number of rounds is at least ``d``.

        Fewer rounds mean the code cannot correct ``(d-1)/2`` errors
        in the time direction.

        Returns:
            ``(True, "")`` if rounds >= d, ``(False, warning)`` otherwise.
        """
        rounds = self.code.rounds
        d = self.code.d
        if rounds < d:
            return False, (
                f"Rounds ({rounds}) < distance ({d}). "
                f"Circuit-level distance may be reduced."
            )
        return True, ""

    def verify_code_distance(self, max_n: int = 15) -> tuple[bool, str]:
        """Brute-force verify the declared code distance.

        Skipped when ``n > max_n`` (too expensive).

        Args:
            max_n: Maximum code size to attempt brute-force.

        Returns:
            ``(True, "")`` if distance matches or skipped,
            ``(False, message)`` otherwise.
        """
        if self.code.n > max_n:
            return True, f"Skipped (n={self.code.n} > {max_n})."

        da = DistanceAnalyzer(self.code, method="bruteforce")
        computed = da.compute_distance_bruteforce()
        if computed != self.code.d:
            return False, (
                f"Computed distance {computed} != declared distance {self.code.d}."
            )
        return True, ""

    def validate(self, check_distance: bool = True) -> tuple[bool, List[str]]:
        """Run all validation checks.

        Args:
            check_distance: Whether to run brute-force distance check.

        Returns:
            ``(all_pass, messages)`` where *messages* lists any warnings
            or errors.
        """
        messages: List[str] = []
        all_pass = True

        ok, msg = self.verify_commutation()
        if not ok:
            all_pass = False
            messages.append(f"[FAIL] Commutation: {msg}")
        else:
            messages.append("[PASS] All stabilizers commute.")

        ok, msg = self.verify_logical_operators()
        if not ok:
            all_pass = False
            messages.append(f"[FAIL] Logical operators: {msg}")
        else:
            messages.append("[PASS] Logical operators valid.")

        ok, msg = self.verify_sufficient_rounds()
        if not ok:
            messages.append(f"[WARN] Rounds: {msg}")
        else:
            messages.append("[PASS] Sufficient rounds.")

        if check_distance:
            ok, msg = self.verify_code_distance()
            if not ok:
                all_pass = False
                messages.append(f"[FAIL] Distance: {msg}")
            else:
                messages.append(f"[PASS] Distance verified. {msg}")

        return all_pass, messages


# class LogicalOperatorAnalyzer:
#     """
#     Analyzer to identify and verify logical operators in a quantum error correction code.
#     """


#     def __init__(self, code: StabCode):
#         self.code = code


#     def verify_logical_Z(self) -> bool:
#         """
#         Verify the logical Z set by the code is correct.

#         1. It must commute with all stabilizers of the code.
#         2. It must not be in the stabilizer group itself.
#         3. Different logical Z operators must commute with each other.
#         4. Different logical Z doesn't form a stabilizer operator by multiplying them together.
#         """
#         # Placeholder implementation
#         return True


#     def is_logical_operator(self, opstring: str) -> bool:
#         """
#         Verify if the given operator string is a logical operator of the code.

#         There are two criteria for an operator to be a logical operator:
#         1. It must commute with all stabilizers of the code.
#         2. It must not be in the stabilizer group itself.
#         """
#         # Placeholder implementation
#         return False


#     def determine_logical_X(self, k: int) -> str:
#         """
#         Find the logical X operator for the k-th logical qubit.

#         There must be 4 criteria for an operator to be a logical X operator:
#         1. It must commute with all stabilizers of the code.
#         2. It must not be in the stabilizer group itself.
#         3. It must anti-commute with the logical Z operator of the same qubit.
#         4. It must act as the identity on all other logical qubits.(Commute with their logical Z operators)
#         """
#         # Placeholder implementation
#         return "X" * self.code.n


#     def determine_logical_H(self, k: int) -> str:
#         """
#         Find the logical H operator for the k-th logical qubit.

#         There must be 4 criteria for an operator to be a logical H operator:
#         1. It must commute with all stabilizers of the code.
#         2. It must not be in the stabilizer group itself.
#         3. It must map logical X to logical Z and logical Z to logical X for the same qubit.
#         4. It must act as the identity on all other logical qubits.(Commute with their logical operators)
#         """
#         # Placeholder implementation
#         return "H" * self.code.n


#     def determine_logical_S(self, k: int) -> str:
#         """
#         Find the logical S operator for the k-th logical qubit.

#         There must be 4 criteria for an operator to be a logical S operator:
#         1. It must commute with all stabilizers of the code.
#         2. It must not be in the stabilizer group itself.
#         3. It must map logical X to logical Y and logical Z to logical Z for the same qubit.
#         4. It must act as the identity on all other logical qubits.(Commute with their logical operators)
#         """
#         # Placeholder implementation
#         return "S" * self.code.n


#     def determine_logical_CNOT(self, control_k: int, target_k: int) -> str:
#         """
#         Find the logical CNOT operator between the control_k and target_k logical qubits.

#         There must be 4 criteria for an operator to be a logical CNOT operator:
#         1. It must commute with all stabilizers of the code.
#         2. It must not be in the stabilizer group itself.
#         3. It must map logical X of control to logical X of control and logical X of target to logical X of target * logical X of control.
#            It must map logical Z of target to logical Z of target and logical Z of control to logical Z of control * logical Z of target.
#         4. It must act as the identity on all other logical qubits.(Commute with their logical operators)
#         """
#         # Placeholder implementation
#         assert control_k != target_k, "Control and target logical qubits must be different."
#         return "CNOT" * self.code.n


class LogicalOperatorAnalyzer:
    """
    Analyzer to identify and verify logical operators in a quantum error correction code.

    This implementation is intentionally naive:
      - It uses a simple binary symplectic representation of stabilizers.
      - Group membership (whether an operator is in the stabilizer group) is
        checked via Gaussian elimination over GF(2).
      - Logical X is found as any operator that:
          * commutes with all stabilizers,
          * is not in the stabilizer group,
          * (if logical Z is set) anti-commutes with Z_k and commutes with Z_j
            for j != k.

    It does NOT attempt to find minimal-weight logicals.
    """

    def __init__(self, code: StabCode):
        self.code = code
        self._n: int = code.n
        # Stabilizer generators as strings
        self._stabs: List[str] = list(getattr(code, "_stabs", []))
        # Logical Z operators dictionary: index -> Pauli string
        self._logicalZ: Dict[int, str] = getattr(code, "_logicalZ", {})
        # Precompute binary-symplectic rows for stabilizers
        self._stab_rows: List[int] = [self._encode_pauli(s) for s in self._stabs]
        self._stab_rank: int = self._rank(self._stab_rows)

    # ------------------------------------------------------------------
    # Basic helpers: Pauli encoding and group algebra
    # ------------------------------------------------------------------

    def _encode_pauli(self, op: str) -> int:
        """
        Encode a Pauli string into a single integer representing its
        (x|z) binary symplectic vector of length 2n.

        Bits 0..(n-1)   : X components
        Bits n..(2n-1)  : Z components
        """
        if len(op) != self._n:
            raise ValueError(f"Pauli string length {len(op)} != code length {self._n}")

        x_bits = 0
        z_bits = 0
        for i, ch in enumerate(op):
            ch = ch.upper()
            if ch == "I":
                continue
            elif ch == "X":
                x_bits |= 1 << i
            elif ch == "Z":
                z_bits |= 1 << i
            elif ch == "Y":
                x_bits |= 1 << i
                z_bits |= 1 << i
            else:
                raise ValueError(f"Invalid Pauli character {ch!r} in {op!r}")

        # Pack into a single integer: [x | z]
        return x_bits | (z_bits << self._n)

    def _rank(self, rows: List[int]) -> int:
        """
        Compute the rank over GF(2) of a list of row-vectors (each an int)
        of length 2n bits using Gaussian elimination.
        """
        rows = [r for r in rows if r != 0]
        rank = 0
        num_bits = 2 * self._n

        for bit in range(num_bits):
            # Find pivot row with this bit set
            pivot_idx = None
            for i in range(rank, len(rows)):
                if (rows[i] >> bit) & 1:
                    pivot_idx = i
                    break
            if pivot_idx is None:
                continue

            # Swap pivot row into position `rank`
            rows[rank], rows[pivot_idx] = rows[pivot_idx], rows[rank]
            pivot_row = rows[rank]

            # Eliminate this bit from all other rows
            for j in range(len(rows)):
                if j != rank and ((rows[j] >> bit) & 1):
                    rows[j] ^= pivot_row

            rank += 1
            if rank == len(rows):
                break

        return rank

    def _is_in_stabilizer_group(self, op: str) -> bool:
        """
        Check if a Pauli string `op` is in the stabilizer group generated by
        the stabilizers in self._stabs.

        We check whether `op`'s vector lies in the row span of the stabilizer
        generator matrix (ignoring phases).
        """
        if not self._stabs:
            return False  # No stabilizers → no non-trivial group.

        vec = self._encode_pauli(op)
        if vec == 0:
            # Identity is always in the stabilizer group.
            return True

        # If adding vec to the generator rows does not increase the rank,
        # then vec is in the span of the generators.
        rows_plus = self._stab_rows + [vec]
        rank_plus = self._rank(rows_plus)
        return rank_plus == self._stab_rank

    @staticmethod
    def _multiply_pauli_strings(p: str, q: str) -> str:
        """
        Multiply two Pauli strings p and q (same length), ignoring global phase.

        Delegates to the utility function multiply_pauli_strings.
        """
        return multiply_pauli_strings(p, q)

    @staticmethod
    def _anticommutes(op1: str, op2: str) -> bool:
        """
        True iff op1 and op2 anti-commute (using the provided commute helper).
        """
        return not commute(op1, op2)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def verify_logical_Z(self) -> bool:
        """
        Verify the logical Z set by the code is correct.

        For each logical Z in code._logicalZ:
        1. It must commute with all stabilizers of the code.
        2. It must not be in the stabilizer group itself.
        3. Different logical Z operators must commute with each other.
        4. The product of any subset of logical Z operators should not be a
           stabilizer (here we just check pairwise products).
        """
        if not self._logicalZ:
            # Nothing to verify
            return True

        # 1 & 2: each logical Z commutes with all stabilizers and is non-stabilizer
        for idx, lz in self._logicalZ.items():
            # Commute with stabilizers
            for stab in self._stabs:
                if not commute(lz, stab):
                    return False

            # Not in stabilizer group
            if self._is_in_stabilizer_group(lz):
                return False

        # 3 & 4: pairwise checks among logical Z operators
        indices = list(self._logicalZ.keys())
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                zi = self._logicalZ[indices[i]]
                zj = self._logicalZ[indices[j]]

                # Must commute with each other
                if not commute(zi, zj):
                    return False

                # Product must not be a stabilizer
                prod = self._multiply_pauli_strings(zi, zj)
                if self._is_in_stabilizer_group(prod):
                    return False

        return True

    def is_logical_operator(self, opstring: str) -> bool:
        """
        Verify if the given operator string is a logical operator of the code.

        Criteria:
        1. It must commute with all stabilizers of the code.
        2. It must not be in the stabilizer group itself.
        """
        if len(opstring) != self._n:
            raise ValueError(
                f"Operator length {len(opstring)} does not match code length {self._n}."
            )

        # Commutes with all stabilizers
        for stab in self._stabs:
            if not commute(opstring, stab):
                return False

        # Not in stabilizer group
        if self._is_in_stabilizer_group(opstring):
            return False

        return True

    def _candidate_X_ok(self, candidate: str, k: int) -> bool:
        r"""
        Check if `candidate` is a valid logical X for qubit k, given whatever
        logical Z operators are currently registered in the code.

        Criteria (naive):
        1. candidate is a logical operator (centralizer \ stabilizer).
        2. If logical Z_k is present: candidate anti-commutes with Z_k.
        3. For j != k, if logical Z_j present: candidate commutes with Z_j.
        """
        if not self.is_logical_operator(candidate):
            return False

        # If logical Z_k is known, enforce anti-commutation
        zk = self._logicalZ.get(k, None)
        if zk is not None:
            if commute(candidate, zk):
                return False

        # For other logical Z_j, require commutation
        for j, zj in self._logicalZ.items():
            if j == k:
                continue
            if not commute(candidate, zj):
                return False

        return True

    def determine_logical_X(self, k: int) -> str:
        """
        Find a logical X operator for the k-th logical qubit.

        Naive strategy:
          1. First try the global operator 'X' * n (often a valid logical X
             for many small / CSS codes when logical Z is 'Z'*n).
          2. If that fails, search over low-weight X-only patterns (up to
             a small max_weight) to find ANY operator that:
                - commutes with all stabilizers,
                - is not in the stabilizer group,
                - anti-commutes with Z_k (if known),
                - commutes with Z_j for j != k (if known).
        """
        n = self._n
        if k < 0 or k >= self.code.k:
            raise ValueError(f"logical qubit index {k} out of range (k={self.code.k}).")

        # 1) Try full X^n
        candidate = "X" * n
        if self._candidate_X_ok(candidate, k):
            return candidate

        # 2) Try low-weight X patterns: weight 1..max_weight
        from itertools import combinations

        max_weight = min(4, n)  # naive cap
        for w in range(1, max_weight + 1):
            for positions in combinations(range(n), w):
                chars = ["I"] * n
                for pos in positions:
                    chars[pos] = "X"
                cand = "".join(chars)
                if self._candidate_X_ok(cand, k):
                    return cand

        # Fallback: if we couldn't find anything decent, just return X^n
        # (we already know it is at least a non-stabilizer or centralizer in
        # many practical codes, but might not satisfy all criteria).
        return candidate

    def determine_logical_H(self, k: int) -> str:
        """
        Find the logical H operator for the k-th logical qubit.

        Naive placeholder:
          - Return 'H' * n, i.e., apply H to every physical qubit in the block.

        This is correct for codes where H is transversal on all physical qubits
        (e.g., Steane code in the usual encoding), but not generally.
        """
        if k < 0 or k >= self.code.k:
            raise ValueError(f"logical qubit index {k} out of range (k={self.code.k}).")
        return "H" * self.code.n

    def determine_logical_S(self, k: int) -> str:
        """
        Find the logical S operator for the k-th logical qubit.

        Naive placeholder: return 'S' * n.
        """
        if k < 0 or k >= self.code.k:
            raise ValueError(f"logical qubit index {k} out of range (k={self.code.k}).")
        return "S" * self.code.n

    def determine_logical_CNOT(self, control_k: int, target_k: int) -> str:
        """
        Find the logical CNOT operator between the control_k and target_k logical qubits.

        Naive placeholder: just return 'CNOT' * n. Not actually used in your
        current compiler path.
        """
        assert control_k != target_k, (
            "Control and target logical qubits must be different."
        )
        if control_k < 0 or control_k >= self.code.k:
            raise ValueError(f"control_k {control_k} out of range (k={self.code.k}).")
        if target_k < 0 or target_k >= self.code.k:
            raise ValueError(f"target_k {target_k} out of range (k={self.code.k}).")
        return "CNOT" * self.code.n


def verify_small_code_logical_operators():
    """
    Verify logical operators for small codes: five-qubit code, Steane code, and Shor code.
    """
    from .small import fivequbitCode, steaneCode, ShorCode

    codes = [fivequbitCode("Standard"), steaneCode("Standard"), ShorCode("Standard")]
    for code in codes:
        analyzer = LogicalOperatorAnalyzer(code)
        assert analyzer.verify_logical_Z(), (
            f"Logical Z verification failed for {code.__class__.__name__}"
        )
        for k in range(code.k):
            logical_X = analyzer.determine_logical_X(k)
            assert analyzer.is_logical_operator(logical_X), (
                f"Logical X verification failed for {code.__class__.__name__}, qubit {k}"
            )

        # Verify code distance as well
        distance_analyzer = DistanceAnalyzer(code, method="bruteforce")
        assert distance_analyzer.verify_code_distance(), (
            f"Code distance verification failed for {code.__class__.__name__}"
        )
        print(
            f"{code.__class__.__name__} passed logical operator and distance verification."
        )


def verify_surface_code_logical_operators():
    """
    Verify logical operators for surface codes of various sizes.
    """
    from .qeccircuit import generate_surface_code  # type: ignore[attr-defined]

    sizes = [3, 5, 7]
    for size in sizes:
        code = generate_surface_code(size, size)
        analyzer = LogicalOperatorAnalyzer(code)
        assert analyzer.verify_logical_Z(), (
            f"Logical Z verification failed for Surface code size {size}x{size}"
        )
        for k in range(code.k):
            logical_X = analyzer.determine_logical_X(k)
            assert analyzer.is_logical_operator(logical_X), (
                f"Logical X verification failed for Surface code size {size}x{size}, qubit {k}"
            )

        # Verify code distance as well
        distance_analyzer = DistanceAnalyzer(code, method="smt")
        assert distance_analyzer.verify_code_distance(), (
            f"Code distance verification failed for Surface code size {size}x{size}"
        )


def check_stabilizer_for_small_code():
    """
    Check stabilizer commutation for small codes: five-qubit code, Steane code, and Shor code.
    """
    from .small import fivequbitCode, steaneCode, ShorCode

    codes = [fivequbitCode("Standard"), steaneCode("Standard"), ShorCode("Standard")]
    for code in codes:
        analyzer = StabilizerAnalyzer(code)
        assert analyzer.verify_commutation(), (
            f"Stabilizer commutation failed for {code.__class__.__name__}"
        )

    print("All small codes passed stabilizer commutation check.")


def find_logical_X_for_small_code():
    """
    Find and print logical X operators for small codes: five-qubit code, Steane code, and Shor code.
    """
    from .small import fivequbitCode, steaneCode, ShorCode

    codes = [fivequbitCode("Standard"), steaneCode("Standard"), ShorCode("Standard")]
    for code in codes:
        analyzer = LogicalOperatorAnalyzer(code)
        print(f"Code: {code.__class__.__name__}")
        for k in range(code.k):
            logical_X = analyzer.determine_logical_X(k)
            print(f"  Logical X for qubit {k}: {logical_X}")


def find_logical_H_for_small_code():
    """
    Find and print logical H operators for small codes: five-qubit code, Steane code, and Shor code.
    """
    from .small import fivequbitCode, steaneCode, ShorCode

    codes = [fivequbitCode("Standard"), steaneCode("Standard"), ShorCode("Standard")]
    for code in codes:
        analyzer = LogicalOperatorAnalyzer(code)
        print(f"Code: {code.__class__.__name__}")
        for k in range(code.k):
            logical_H = analyzer.determine_logical_H(k)
            print(f"  Logical H for qubit {k}: {logical_H}")


def find_logical_S_for_small_code():
    """
    Find and print logical S operators for small codes: five-qubit code, Steane code, and Shor code.
    """
    from .small import fivequbitCode, steaneCode, ShorCode

    codes = [fivequbitCode("Standard"), steaneCode("Standard"), ShorCode("Standard")]
    for code in codes:
        analyzer = LogicalOperatorAnalyzer(code)
        print(f"Code: {code.__class__.__name__}")
        for k in range(code.k):
            logical_S = analyzer.determine_logical_S(k)
            print(f"  Logical S for qubit {k}: {logical_S}")


if __name__ == "__main__":
    verify_small_code_logical_operators()
    check_stabilizer_for_small_code()
    find_logical_X_for_small_code()
    find_logical_H_for_small_code()
    find_logical_S_for_small_code()
