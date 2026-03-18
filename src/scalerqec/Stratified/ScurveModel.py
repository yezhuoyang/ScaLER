from __future__ import annotations

"""Legacy S-curve model functions for stratified LER estimation.

This module contains the original standalone S-curve functions used by
:class:`~scalerqec.Stratified.stratifiedScurveLER.StratifiedScurveLERcalc`.
They implement the log-logit transformed linear model::

    y(w) = log(0.5 / P_L(w) - 1) = a * w + b + c / sqrt(w - t)

which corresponds to the S-curve in probability space::

    P_L(w) = 0.5 / (1 + exp(a * w + b + c / sqrt(w - t)))

.. deprecated::
    For new code, prefer the model classes in
    :mod:`scalerqec.Stratified.models` (``OurScurveModel``,
    ``IBMScurveModel``), which provide an abstract base class, factory
    pattern, and tunable gamma parameter for sweet-spot calculation.
"""

import warnings
import numpy as np
from typing import Callable
from numpy.typing import NDArray


def scurve_function(x: float, center: float, sigma: float) -> float:
    """Evaluate the basic S-curve (logistic) model.

    Computes ``P_L(x) = 0.5 / (1 + exp(-(x - center) / sigma))``.

    Args:
        x: Weight value (number of faults).
        center: Midpoint of the S-curve (weight at which P_L = 0.25).
        sigma: Scale parameter controlling the steepness of the curve.

    Returns:
        Conditional logical error probability in (0, 0.5).
    """
    return float(0.5 / (1 + np.exp(-(x - center) / sigma)))


def inv_logit_half(y: float) -> float:
    """Inverse log-logit transform: map y back to P_L.

    Computes ``P_L = 0.5 / (1 + exp(y))``.

    Args:
        y: Value in the log-logit transformed space.

    Returns:
        Conditional logical error probability in (0, 0.5).
    """
    return float(0.5 / (1 + np.exp(y)))


def linear_function(
    x: float | NDArray[np.floating], a: float, b: float
) -> float | NDArray[np.floating]:
    """Linear model ``y = a * x + b`` used as an initial fit for the S-curve.

    In the log-logit space the S-curve is approximately linear far from
    the pole at *t*, so this function is used for an initial parameter
    estimate before the full modified-linear fit.

    Args:
        x: Weight value(s).
        a: Slope (related to -1/alpha).
        b: Intercept (related to mu/alpha).

    Returns:
        Predicted y value(s) in log-logit space.
    """
    return a * x + b  # type: ignore[return-value]


def modified_linear_function_with_d(
    x: float | NDArray[np.floating], a: float, b: float, c: float, d: float
) -> float | NDArray[np.floating]:
    """Modified linear model with a pole term: ``y = a*x + b + c/sqrt(x-d)``.

    This is the log-logit transform of the Definition 1 S-curve.  The
    ``1/sqrt(x - d)`` term creates a pole at *d = (code_distance-1)/2*,
    ensuring P_L -> 0 for weights at or below the fault-tolerant threshold.

    Args:
        x: Weight value(s). Must satisfy ``x > d``.
        a: Slope parameter (= -1/alpha).
        b: Intercept parameter (= mu/alpha).
        c: Pole strength (= beta).
        d: Pole location, typically ``(code_distance - 1) / 2``.

    Returns:
        Predicted y value(s) in log-logit space.
    """
    eps = 1e-12
    delta = (x - d) ** 0.5
    delta = np.where(np.abs(delta) < eps, np.sign(delta) * eps, delta)
    return a * x + b + c / delta


def modified_linear_function(
    d: float,
) -> Callable[[float, float, float, float, float], float]:
    """Return a curried version of :func:`modified_linear_function_with_d` with *d* fixed.

    This is used with ``scipy.optimize.curve_fit`` where the pole
    location *d* (fault-tolerant threshold) is known a priori and should
    not be optimized.

    Args:
        d: Fixed pole location, typically ``(code_distance - 1) / 2``.

    Returns:
        A callable ``f(x, a, b, c)`` suitable for ``curve_fit``.
    """

    def tempfunc(x: float, a: float, b: float, c: float, d: float = d) -> float:
        return modified_linear_function_with_d(x, a, b, c, d)

    return tempfunc


def modified_sigmoid_function(
    x: float | NDArray[np.floating], a: float, b: float, c: float, d: float
) -> float | NDArray[np.floating]:
    """Evaluate the Definition 1 S-curve in probability space.

    Computes ``P_L(x) = 0.5 / (1 + exp(a*x + b + c/sqrt(x-d)))``.

    This is the inverse log-logit of
    :func:`modified_linear_function_with_d`.

    Args:
        x: Weight value(s). Must satisfy ``x > d``.
        a: Slope parameter (= -1/alpha).
        b: Intercept parameter (= mu/alpha).
        c: Pole strength (= beta).
        d: Pole location, typically ``(code_distance - 1) / 2``.

    Returns:
        Conditional logical error probability P_L in (0, 0.5).
    """
    z = a * x + b + c / ((x - d) ** 0.5)
    # ignore overflows in exp → exp(z) becomes np.inf, so 0.5/(1+inf) = 0.0
    with np.errstate(over="ignore"):
        y = 0.5 / (1 + np.exp(z))
    return y


def quadratic_function(x: float, a: float, b: float, c: float) -> float:
    """Quadratic model ``y = a*x^2 + b*x + c`` for curve fitting.

    Args:
        x: Weight value.
        a: Quadratic coefficient.
        b: Linear coefficient.
        c: Constant term.

    Returns:
        Predicted y value.
    """
    return a * x**2 + b * x + c


def poly_function(x: float, a: float, b: float, c: float, d: float) -> float:
    """Cubic polynomial ``y = a*x^3 + b*x^2 + c*x + d`` for curve fitting.

    Args:
        x: Weight value.
        a: Cubic coefficient.
        b: Quadratic coefficient.
        c: Linear coefficient.
        d: Constant term.

    Returns:
        Predicted y value.
    """
    return a * x**3 + b * x**2 + c * x + d


def refined_sweet_spot(
    alpha: float, beta: float, t: float, ratio: float = 0.05
) -> float:
    """Compute the sweet-spot weight where the pole term becomes negligible.

    The sweet spot is defined as the weight *w* at which the second
    derivative of ``y(w)`` (from the ``1/sqrt(w-t)`` term) drops to
    *ratio* times the first derivative (from the linear term).  Solving
    ``1/alpha = ratio * (1/2) * beta / (w - t)^{3/2}`` gives::

        w = t + [(ratio * beta * alpha / 2)]^{2/3}

    Below the sweet spot, the pole term dominates and the S-curve is
    highly nonlinear, making it the most informative weight for fitting.

    Args:
        alpha: Scale parameter of the S-curve (= -1/a).
        beta: Pole-strength parameter (= c).
        t: Fault-tolerant threshold ``(code_distance - 1) / 2``.
        ratio: Fraction of the linear-term derivative at which the pole
            term is considered negligible. Defaults to 0.05.

    Returns:
        Sweet-spot weight as a float.
    """
    return float(t + ((ratio * beta * alpha / 2) ** (2 / 3)))


def sigma_estimator(N: int, M: int) -> float:
    """Estimate the standard deviation of the log-logit statistic ``y(w)``.

    Given *N* samples of which *M* are logical errors (with
    ``P_w = M/N``), this returns the delta-method standard error of
    ``y = log(0.5/P_w - 1)`` used as weights in the weighted
    least-squares S-curve fit.

    Args:
        N: Total number of samples at a given weight.
        M: Number of observed logical errors.

    Returns:
        Estimated standard deviation of ``y(w)``.
    """
    return float(np.sqrt(N**2 * (N - M) / (M * (N - 1) * (N - 2 * M) ** 2)))


def subspace_sigma_estimator(N: int, M: int) -> float:
    """Estimate the standard deviation of the subspace LER ``P_w = M/N``.

    Uses the finite-population corrected binomial standard error.

    Args:
        N: Total number of samples at a given weight.
        M: Number of observed logical errors.

    Returns:
        Estimated standard deviation of ``P_w``.
    """
    return float(np.sqrt(M * (N - M) / (N - 1)) / N)


def bias_estimator(N: int, M: int) -> float:
    """Estimate the second-order bias of the log-logit statistic ``y(w)``.

    Because ``y = f(P_w) = log(1/(2*P_w) - 1)`` is a nonlinear
    transformation of the sample proportion ``P_w = M/N``, its
    expectation has a bias of approximately
    ``(1/2) * f''(P_w) * Var(P_w)``.  This function returns that
    bias estimate so it can be subtracted during curve fitting.

    Args:
        N: Total number of samples at a given weight.
        M: Number of observed logical errors.

    Returns:
        Estimated bias of ``y(w)``.
    """
    return 1 / 2 * (N / M) * (N - 4 * M) / (N - 2 * M) ** 2 * (N - M) / (N - 1)


def show_bias_estimator(N: int, M: int) -> float:
    """Estimate the bias of the observed subspace LER ``P_w = M/N``.

    A simpler, first-order bias estimate of P_w itself (rather than
    the log-logit transform).

    Args:
        N: Total number of samples at a given weight.
        M: Number of observed logical errors.

    Returns:
        Estimated bias of ``P_w``.
    """
    Pw = M / N
    return (1 - Pw) / (2 * Pw * N)


def evenly_spaced_ints(minw: int, maxw: int, N: int) -> list[int]:
    """Return *N* approximately evenly spaced integers in [minw, maxw].

    Used to choose which error weights to sample when a uniform spread
    of measurement points is desired across the S-curve.

    Args:
        minw: Lower bound of the weight range (inclusive).
        maxw: Upper bound of the weight range (inclusive).
        N: Desired number of points.

    Returns:
        Sorted list of *N* distinct integers spanning [minw, maxw].
        If ``N`` exceeds the number of available integers, all integers
        in the range are returned.
    """
    if N == 1:
        return [minw]
    if N > (maxw - minw + 1):
        return list(range(minw, maxw + 1))

    # Use high-resolution linspace, round, then deduplicate
    raw = np.linspace(minw, maxw, num=10 * N)
    rounded = sorted(set(map(int, raw)))

    # Pick N evenly spaced indices from the unique set
    indices = np.linspace(0, len(rounded) - 1, num=N, dtype=int)
    return [rounded[i] for i in indices]
