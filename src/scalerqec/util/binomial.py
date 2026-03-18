"""Binomial probability and combinatorial utilities.

Provides functions for computing binomial probability mass values and
subspace sizes, used primarily in the ScaLER algorithm for weighting
error-weight subspaces.
"""

import math
from scipy.stats import binom


def subspace_size(num_noise: int, weight: int) -> int:
    """Compute the number of error patterns of a given weight.

    Returns the binomial coefficient :math:`\\binom{N}{w}`, which is the
    number of ways to choose ``weight`` noise sources out of
    ``num_noise`` total sources.

    Args:
        num_noise (int): Total number of independent noise sources (N).
        weight (int): Number of noise sources that fire (w).

    Returns:
        int: The binomial coefficient :math:`\\binom{N}{w}`.
    """
    return math.comb(num_noise, weight)


def binomial_weight(N: int, W: int, p: float) -> float:
    """Compute the binomial probability mass function value.

    Returns the probability that exactly ``W`` out of ``N`` independent
    Bernoulli trials (each with success probability ``p``) succeed:

    .. math::

        P(X = W) = \\binom{N}{W} p^W (1-p)^{N-W}

    Uses ``scipy.stats.binom.pmf`` for numerical stability with large
    ``N``.

    Args:
        N (int): Total number of trials (noise sources).
        W (int): Number of successes (errors).
        p (float): Probability of success (physical error rate) per
            trial, in the range [0, 1].

    Returns:
        float: The binomial PMF value :math:`P(X = W)`.
    """
    return float(binom.pmf(W, N, p))
