"""Brute-force (naive) symbolic LER computation for validation.

This module computes the exact logical error rate by explicitly enumerating
all :math:`4^N` possible Pauli error patterns (I, X, Y, Z at each of
``N`` noise sources).  For the hardcoded 4-noise-source circuit, this means
:math:`4^4 = 256` patterns.

For each pattern the module:

1. Determines the combined detector/observable outcome by XOR-ing the
   propagation vectors of the active errors.
2. Checks whether the observable bit is 1 (indicating a logical error).
3. Accumulates the probability :math:`(p/3)^w (1-p)^{N-w}` where ``w``
   is the number of non-identity errors.

The result is an exact polynomial in ``p`` that should match the output
of the DP-based :mod:`~scalerqec.Symbolic.symbolic` module, serving as
an independent cross-check.

Warning:
    The :math:`4^N` enumeration is exponential and only feasible for very
    small circuits (N < ~15).  For production use, prefer
    :class:`~scalerqec.Symbolic.symbolicLER.SymbolicLERcalc`.
"""

from sympy import simplify, symbols

# ----------------------------------------------------------------------
# Physical-error model
p = symbols("p")
q = 1 - p  #   probability of "no error"
Px = Py = Pz = p / 3  #   probabilities of  X  Y  Z


# Syndrome-propagation vectors (two bits, XOR addition)
# index encoding  D0 + 2*O0   (00->0, 01->1, 10->2, 11->3)
def vec_to_idx(vec):
    """Convert a 2-element binary vector to an integer index.

    Args:
        vec (tuple[int, int]): Binary vector ``(D0, O0)``.

    Returns:
        int: Integer index ``D0*2 + O0``.
    """
    return vec[1] + 2 * vec[0]


def idx_to_vec(idx):
    """Convert an integer index (0..3) to a 2-element binary vector.

    Args:
        idx (int): Integer in range 0..3.

    Returns:
        tuple[int, int]: Binary vector ``(D0, O0)``.
    """
    return ((idx >> 1) & 1, idx & 1)


def xor_vec(a, b):
    """Element-wise XOR of two 2-element binary vectors.

    Args:
        a (tuple[int, int]): First binary vector.
        b (tuple[int, int]): Second binary vector.

    Returns:
        tuple[int, int]: Element-wise XOR result.
    """
    return (a[0] ^ b[0], a[1] ^ b[1])


def pos_int_to_vec(vecidx):
    """Convert a 4-bit integer to a list of its individual bits.

    Returns bits in MSB-first order.

    Args:
        vecidx (int): Integer in range 0..15.

    Returns:
        list[int]: List of 4 binary digits.
    """
    return [vecidx >> (3 - i) & 1 for i in range(4)]


PROP_X = [(1, 1), (0, 1), (1, 0), (0, 1)]
"""list[tuple[int, int]]: X-error propagation vectors for each noise source."""

PROP_Y = [(1, 1), (0, 1), (1, 0), (0, 1)]
"""list[tuple[int, int]]: Y-error propagation vectors for each noise source."""

PROP_Z = [(0, 0), (0, 0), (0, 0), (0, 0)]
"""list[tuple[int, int]]: Z-error propagation vectors (all zero for this circuit)."""

PROPAGATORS = (PROP_X, PROP_Y, PROP_Z)
"""tuple: Combined propagation vectors for X, Y, Z errors."""


MAX_degree = 100
"""int: Maximum polynomial degree retained in series truncation."""


LER = p * 0
# We use an 4bit integer to represent the position of the pauli noise
for n0 in range(0, 4):
    n0count = n0 != 0
    for n1 in range(0, 4):
        n1count = n1 != 0
        for n2 in range(0, 4):
            n2count = n2 != 0
            for n3 in range(0, 4):
                n3count = n3 != 0
                count = n0count + n1count + n2count + n3count

                print("{},{},{},{}".format(n0, n1, n2, n3))

                init_vec = (0, 0)
                if n0 == 1:
                    init_vec = xor_vec(init_vec, PROP_X[0])
                elif n0 == 2:
                    init_vec = xor_vec(init_vec, PROP_Y[0])
                elif n0 == 3:
                    init_vec = xor_vec(init_vec, PROP_Z[0])

                if n1 == 1:
                    init_vec = xor_vec(init_vec, PROP_X[1])
                elif n1 == 2:
                    init_vec = xor_vec(init_vec, PROP_Y[1])
                elif n1 == 3:
                    init_vec = xor_vec(init_vec, PROP_Z[1])

                if n2 == 1:
                    init_vec = xor_vec(init_vec, PROP_X[2])
                elif n2 == 2:
                    init_vec = xor_vec(init_vec, PROP_Y[2])
                elif n2 == 3:
                    init_vec = xor_vec(init_vec, PROP_Z[2])

                if n3 == 1:
                    init_vec = xor_vec(init_vec, PROP_X[3])
                elif n3 == 2:
                    init_vec = xor_vec(init_vec, PROP_Y[3])
                elif n3 == 3:
                    init_vec = xor_vec(init_vec, PROP_Z[3])

                print(init_vec)

                if init_vec[1] == 1:
                    LER += simplify((p / 3) ** count * (1 - p) ** (4 - count))
            # keep only terms with deg(p) < MAX_degree
    LER = LER.series(p, 0, MAX_degree).removeO()  # no .expand()


print(LER)
