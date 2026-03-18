"""Hardcoded symbolic LER example for a specific 2-qubit, 4-noise circuit.

This module is a **reference / teaching implementation** that demonstrates
the symbolic DP algorithm on a small, fixed circuit with:

- 1 detector (D0)
- 1 observable (O0)
- 4 independent noise sources

The detector/observable outcome is encoded as a 2-bit vector ``(D0, O0)``
mapped to an integer index ``D0*2 + O0`` (range 0..3).  The propagation
vectors ``PROP_X``, ``PROP_Y``, ``PROP_Z`` are hardcoded for this
specific circuit.

The DP table ``dp[i][j][v]`` is filled at module load time and the LER
can be computed by calling :func:`calculate_LER` with the list of
outcome indices that correspond to logical errors.

This module is useful for:

- Understanding the DP recurrence before reading the general
  :class:`~scalerqec.Symbolic.symbolicLER.SymbolicLERcalc`.
- Manually verifying results for a known small circuit.
- Generating LaTeX tables of the DP entries via :func:`print_table`.
"""

from sympy import symbols, binomial, Rational, simplify, latex

# ----------------------------------------------------------------------
# Physical-error model
p = symbols("p")
q = 1 - p  #   probability of "no error"
Px = Py = Pz = p / 3  #   probabilities of  X  Y  Z


# Syndrome-propagation vectors (two bits, XOR addition)
# index encoding  D0 + 2*O0   (00->0, 01->1, 10->2, 11->3)
def vec_to_idx(vec):
    """Convert a 2-element binary vector to an integer index.

    Uses the encoding ``vec[0]*2 + vec[1]`` (i.e., ``(D0, O0) -> D0*2 + O0``).

    Args:
        vec (tuple[int, int]): A pair ``(D0, O0)`` of binary digits.

    Returns:
        int: Integer index in the range 0..3.
    """
    return vec[1] + 2 * vec[0]


def idx_to_vec(idx):
    """Convert an integer index (0..3) to a 2-element binary vector.

    Inverse of :func:`vec_to_idx`.

    Args:
        idx (int): Integer in range 0..3.

    Returns:
        tuple[int, int]: Binary vector ``(D0, O0)``.
    """
    return ((idx >> 1) & 1, idx & 1)


PROP_X = [(1, 1), (0, 1), (1, 0), (0, 1)]
"""list[tuple[int, int]]: X-error propagation vectors for each of the 4 noise sources."""

PROP_Y = [(1, 1), (0, 1), (1, 0), (0, 1)]
"""list[tuple[int, int]]: Y-error propagation vectors for each of the 4 noise sources."""

PROP_Z = [(0, 0), (0, 0), (0, 0), (0, 0)]
"""list[tuple[int, int]]: Z-error propagation vectors (all zero for this circuit)."""

PROPAGATORS = (PROP_X, PROP_Y, PROP_Z)
"""tuple: Combined propagation vectors for X, Y, Z errors."""


MAX_degree = 200
"""int: Maximum polynomial degree retained in series expansions."""


def xor_vec(a, b):
    """Element-wise XOR of two 2-element binary vectors.

    Args:
        a (tuple[int, int]): First binary vector.
        b (tuple[int, int]): Second binary vector.

    Returns:
        tuple[int, int]: Element-wise XOR result.
    """
    return (a[0] ^ b[0], a[1] ^ b[1])


def verity_table(i):
    """Verify that DP row ``i`` sums to 1 (valid probability distribution).

    Sums ``dp[i][j][v]`` over all weights ``j`` and outcome indices
    ``v``, simplifies symbolically, and asserts the result equals 1.

    Args:
        i (int): The noise-source row index to verify.

    Raises:
        AssertionError: If the probabilities do not sum to 1.
    """
    sum = 0
    for j in range(0, i + 1):
        for vec_index in range(4):
            sum += dp[i][j][vec_index]
        # print(dp[i][j][vec_index])
    sum = simplify(sum)
    print(i, sum)
    assert sum == 1


# ----------------------------------------------------------------------
# DP tables
MAX_I = 4  # change here if you need longer arrays
dp = [
    [[0] * 4 for _ in range(MAX_I + 1)]  # dp[i][j][vecIdx]
    for _ in range(MAX_I + 1)
]

# Base:  i = 0, j = 0 ⇒ probability 1 at vec = (0,0)
dp[0][0][vec_to_idx((0, 0))] = 1


# ----------------------------------------------------------------------
# Fill   dp[i][j][·]   using the recurrence in Eq. (1)
for i in range(1, MAX_I + 1):
    dp[i][0][0] = (1 - p) ** i

    for j in range(1, i + 1):  # j ≤ i
        for vec_idx in range(4):
            vec = idx_to_vec(vec_idx)

            # 1) “no error’’ branch
            acc = q * dp[i - 1][j][vec_idx]

            # 2) X, Y, Z branches  (need j-1 ≥ 0)
            if j >= 1:
                for prob, prop in (
                    (Px, PROP_X[i - 1]),
                    (Py, PROP_Y[i - 1]),
                    (Pz, PROP_Z[i - 1]),
                ):
                    prev_vec = xor_vec(prop, vec)
                    acc += prob * dp[i - 1][j - 1][vec_to_idx(prev_vec)]

            dp[i][j][vec_idx] = simplify(acc)

            # print(f"dp[{i}][{j}][{vec_idx}] = {dp[i][j][vec_idx].series(p, 0, MAX_degree).removeO()}")
    verity_table(i)


def binom(k):
    """Compute the binomial probability of exactly ``k`` errors out of 4.

    Returns :math:`\\binom{4}{k} p^k (1-p)^{4-k}` as a simplified
    SymPy expression.

    Args:
        k (int): Number of errors (0 to 4).

    Returns:
        sympy.Expr: The binomial probability polynomial.
    """
    binomweight = binomial(4, k) * p**k * q ** (4 - k)
    return simplify(binomweight)


def calculate_LER(error_row_indices):
    """Compute the LER polynomial by summing DP entries at error rows.

    For each error weight (1..4), sums the final DP entries
    ``dp[4][weight][row]`` over all rows that correspond to logical
    errors, then combines the per-weight contributions into the total
    LER polynomial.

    Args:
        error_row_indices (list[int]): Outcome indices (0..3) where
            the observable bit indicates a logical error.

    Returns:
        sympy.Expr: The LER as a simplified, expanded polynomial in ``p``.
    """
    LER = 0
    for weight in range(1, 5):
        subLER = 0
        for rowindex in error_row_indices:
            subLER += dp[MAX_I][weight][rowindex]
        print("Weight: {}".format(weight))
        print(subLER.expand())
        LER += simplify(subLER)

    # LER=LER.series(p, 0, MAX_degree).removeO()    # no .expand()
    return simplify(LER).expand()


def print_table(i, j):
    """Print a LaTeX-formatted table of DP entries ``dp[i][j][*]``.

    Generates a LaTeX ``tabular`` environment showing the detector
    outcome, observable outcome, and probability polynomial for each
    of the 4 outcome indices.  Outcomes with ``O0 = 1`` are highlighted
    in red as errors.

    Args:
        i (int): Noise-source index (row in the DP table).
        j (int): Error weight (column in the DP table).
    """
    header = rf"\begin{{table}}[h!]\centering\caption{{Table of $dp[{i}][{j}]$ for all $\vec{{D}}$}}\n"
    header += r"\resizebox{\columnwidth}{!}{%" "\n"
    header += (
        r"\begin{tabular}{|l|l|l|l|l|}"
        "\n   \\hline "
        r"\n Index  & $D_0$ & $O_0$  & Probability\\\n  \hline"
    )
    print(header)

    for idx in range(4):
        D0, O0 = idx_to_vec(idx)
        prob_latex = latex(dp[i][j][idx])
        error_tag = " ({\\color{red} Error})" if O0 == 1 else ""

        print(f" {idx} {error_tag}  & {D0} & {O0}  & ${prob_latex}$\\\\")
        print("  \\hline")

    footer = r"\end{tabular}" "\n}\n" + r"\end{table}" + "\n"
    print(footer)


# ----------------------------------------------------------------------
if __name__ == "__main__":
    """
    for i in range(1, MAX_I+1):
        for j in range(1, i+1):
            print_table(i, j)
    """

    print("------------------LER----------------------")
    sum = calculate_LER([1, 3])
    print(sum)
    print(sum.evalf(subs={p: 0.1}))
