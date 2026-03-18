from __future__ import annotations

"""Definition 2 (IBM) S-curve model implementation.

This module implements the "min-fail enclosure" model proposed by IBM
for approximating the logical error rate spectrum::

    P_L(w) = 0.5 * [1 - exp(-2 * f0 * (w / w0) ^ gamma_ibm)]

The model is characterised by a power-law onset at the onset weight
*w0* (typically ``t + 1``), and saturates toward 0.5 for large *w*.

Parameters:
    f0: Failure probability at the onset weight (controls amplitude).
    w0: Onset weight (first weight with non-zero logical errors).
    gamma_ibm: Power exponent controlling the steepness of the curve.

.. note::
    ``gamma_ibm`` is a *model* parameter distinct from the sweet-spot
    tuning parameter ``gamma`` inherited from :class:`ScurveModelBase`.
"""

from typing import Any, Callable, Dict, List, Tuple
import numpy as np
from numpy.typing import NDArray

from scalerqec.Stratified.models.base import ScurveModelBase


class IBMScurveModel(ScurveModelBase):
    """IBM "min-fail enclosure" S-curve model (Definition 2).

    Uses the formula::

        P_L(w) = 0.5 * [1 - exp(-2 * f0 * (w / w0) ^ gamma_ibm)]

    Key differences from :class:`OurScurveModel`:

    * Exponential saturation rather than logistic.
    * Characterised by an onset weight *w0* instead of a pole at *t*.
    * Different asymptotic behaviour near saturation -- approaches 0.5
      as ``1 - exp(-z)`` rather than ``1 / (1 + exp(-z))``.

    Fitted parameters:
        f0 (float): Failure probability at the onset weight, controlling
            how quickly saturation is reached.
        w0 (float): Onset weight (typically ``t + 1``), the first weight
            at which logical errors appear.
        gamma_ibm (float): Power exponent controlling curve steepness.

    Args:
        t: Fault-tolerant threshold ``(d - 1) / 2``.
        gamma: Sweet-spot tuning parameter (base-class parameter, not
            the IBM model exponent).
    """

    def __init__(self, t: int = 0, gamma: float = 0.05):
        super().__init__(t=t, gamma=gamma)
        # Initialize default parameters
        self._params = {
            "f0": 0.5,
            "w0": float(t + 1),
            "gamma_ibm": 1.0,
        }

    @property
    def name(self) -> str:
        return "IBMModel"

    @property
    def param_names(self) -> Tuple[str, ...]:
        return ("f0", "w0", "gamma_ibm")

    def predict(
        self, w: float | NDArray[np.floating[Any]]
    ) -> float | NDArray[np.floating[Any]]:
        """
        Compute P_L(w) = 0.5 * [1 - exp(-2*f0*(w/w0)^γ)].

        For w ≤ t, returns 0 (fault-tolerant region).
        """
        w_arr = np.atleast_1d(np.asarray(w, dtype=float))
        result = np.zeros_like(w_arr)

        f0 = self._params.get("f0", 0.5)
        w0 = self._params.get("w0", float(self._t + 1))
        gamma_ibm = self._params.get("gamma_ibm", 1.0)

        # Avoid division by zero
        if w0 <= 0:
            w0 = 1.0

        # Only compute for w > t (fault-tolerant threshold)
        valid = w_arr > self._t
        if np.any(valid):
            w_valid = w_arr[valid]
            ratio = w_valid / w0

            with np.errstate(over="ignore", under="ignore"):
                exp_term = np.exp(-2 * f0 * np.power(ratio, gamma_ibm))
                result[valid] = 0.5 * (1 - exp_term)

        # Return scalar if input was scalar
        if np.isscalar(w):
            return float(result[0])
        return result

    def transform(
        self, p_w: float | NDArray[np.floating[Any]]
    ) -> float | NDArray[np.floating[Any]]:
        """
        Transform P_L to linear space.

        For IBM model: P = 0.5*(1 - exp(-x)) where x = 2*f0*(w/w0)^γ
        Solving: exp(-x) = 1 - 2*P
        So: x = -log(1 - 2*P)

        We use y = log(0.5/P - 1) for consistency with the fitting framework.
        """
        p_w_arr = np.atleast_1d(np.asarray(p_w, dtype=float))
        # Clip to valid range
        p_clipped = np.clip(p_w_arr, 1e-12, 0.5 - 1e-12)
        result = np.log(0.5 / p_clipped - 1)

        if np.isscalar(p_w):
            return float(result[0])
        return result

    def inverse_transform(
        self, y: float | NDArray[np.floating[Any]]
    ) -> float | NDArray[np.floating[Any]]:
        """
        Transform y back to P_L.
        """
        y_arr = np.atleast_1d(np.asarray(y, dtype=float))

        with np.errstate(over="ignore"):
            result = 0.5 / (1 + np.exp(y_arr))

        if np.isscalar(y):
            return float(result[0])
        return result

    def linear_prediction(
        self, w: float | NDArray[np.floating[Any]]
    ) -> float | NDArray[np.floating[Any]]:
        """
        Compute y(w) in transformed linear space for IBM model.

        y(w) = log(0.5/P(w) - 1) = log((1 + exp(-x))/(1 - exp(-x)))
        where x = 2*f0*(w/w0)^γ

        For w ≤ t, returns inf (fault-tolerant region where P_L = 0).
        """
        w_arr = np.atleast_1d(np.asarray(w, dtype=float))
        result = np.full_like(w_arr, np.inf)

        f0 = self._params.get("f0", 0.5)
        w0 = self._params.get("w0", float(self._t + 1))
        gamma_ibm = self._params.get("gamma_ibm", 1.0)

        if w0 <= 0:
            w0 = 1.0

        # Only compute for w > t (fault-tolerant threshold)
        valid = w_arr > self._t
        if np.any(valid):
            w_valid = w_arr[valid]
            ratio = w_valid / w0

            with np.errstate(over="ignore", under="ignore", divide="ignore"):
                x = 2 * f0 * np.power(ratio, gamma_ibm)
                exp_neg_x = np.exp(-x)
                # Avoid log(0) and log(negative)
                exp_neg_x = np.clip(exp_neg_x, 1e-12, 1 - 1e-12)
                result[valid] = np.log((1 + exp_neg_x) / (1 - exp_neg_x))

        if np.isscalar(w):
            return float(result[0])
        return result

    def _get_fit_function(self) -> Callable[..., Any]:
        """
        Return the function for curve fitting.

        We fit in the transformed space: y = log((1 + exp(-x))/(1 - exp(-x)))
        where x = 2*f0*(w/w0)^γ
        """

        def ibm_linear_function(
            w: float | NDArray[np.floating[Any]], f0: float, w0: float, gamma_ibm: float
        ) -> float | NDArray[np.floating[Any]]:
            w_arr = np.atleast_1d(np.asarray(w, dtype=float))

            if w0 <= 0:
                w0 = 1.0

            ratio = w_arr / w0

            with np.errstate(over="ignore", under="ignore", divide="ignore"):
                x = 2 * f0 * np.power(ratio, gamma_ibm)
                exp_neg_x = np.exp(-x)
                exp_neg_x = np.clip(exp_neg_x, 1e-12, 1 - 1e-12)
                result = np.log((1 + exp_neg_x) / (1 - exp_neg_x))

            if np.isscalar(w):
                return float(result[0])
            return result

        return ibm_linear_function

    def _get_initial_guess(
        self, x_list: List[float], y_list: List[float]
    ) -> Tuple[float, ...]:
        """
        Get initial parameter guess for IBM model.
        """
        # Initial w0 is the minimum weight in the data (approximately t+1)
        w0 = min(x_list) if x_list else float(self._t + 1)
        if w0 <= 0:
            w0 = 1.0

        # Initial f0: use the P value at the smallest weight
        # If y is large, P is small, so f0 should be small
        if y_list:
            y_min = min(y_list)
            # Rough estimate: at w=w0, x = 2*f0, so exp(-2*f0) ~ (1-2*P)
            # If P is very small (y very large), f0 ~ P
            p_at_min = 0.5 / (1 + np.exp(y_min)) if y_min < 100 else 1e-10
            f0 = max(0.01, min(1.0, p_at_min * 10))
        else:
            f0 = 0.5

        # Initial gamma_ibm
        gamma_ibm = 1.0

        return (f0, w0, gamma_ibm)

    def _get_bounds(
        self, initial_guess: Tuple[float, ...]
    ) -> Tuple[List[float], List[float]]:
        """
        Get parameter bounds for IBM model fitting.
        """
        f0, w0, gamma_ibm = initial_guess

        # f0: must be positive, typically < 1
        # w0: must be positive, around the onset weight
        # gamma_ibm: typically between 0.5 and 5

        lower = [
            max(1e-6, f0 * 0.01),
            max(1.0, w0 * 0.5),
            0.1,
        ]
        upper = [
            min(10.0, f0 * 100),
            w0 * 2.0,
            10.0,
        ]

        # Ensure lower < upper
        for i in range(3):
            if lower[i] >= upper[i]:
                lower[i], upper[i] = upper[i], lower[i]
            if lower[i] == upper[i]:
                lower[i] = lower[i] * 0.5
                upper[i] = upper[i] * 2.0

        return lower, upper

    def _set_params_from_fit(self, popt: Tuple[float, ...]) -> None:
        """
        Set model parameters from fitted values.
        """
        self._params = {
            "f0": popt[0],
            "w0": popt[1],
            "gamma_ibm": popt[2],
        }

    def calculate_sweet_spot(self) -> int:
        """
        Calculate the sweet spot for IBM model using d²y/dw² = γ · |dy/dw|.

        For IBM model: P_L(w) = 0.5 * [1 - exp(-z)] where z = 2μ(w/β)^α

        In transformed space: y = log(0.5/P_L - 1) = -z - log(1 - e^(-z))

        For the small-z regime (where P_L is small, which is our regime of interest):
            y ≈ -z - ln(z)
            dy/dw ≈ -α/w
            d²y/dw² ≈ α/w²

        The sweet spot condition d²y/dw² = γ·|dy/dw| gives:
            α/w² = γ · α/w
            w_sweet = 1/γ

        This is clamped to [t+1, w0] to stay in the relevant region.

        Parameters:
            μ (f0): amplitude parameter
            β (w0): onset weight
            α (gamma_ibm): power exponent
            γ (self._gamma): sweet spot tuning parameter
        """
        f0 = self._params.get("f0", 0.5)
        w0 = self._params.get("w0", float(self._t + 1))
        alpha = self._params.get("gamma_ibm", 1.0)

        if self._gamma <= 0:
            return max(1, self._t + 1)

        # In small-z regime: w_sweet = 1/γ
        w_sweet_base = 1.0 / self._gamma

        # Adjust based on model parameters for the transition region
        # When z is moderate, use: w_sweet = w0 * (1/γ / w0)^(1/α) if that's smaller
        # This accounts for the power-law scaling
        if alpha > 1 and w0 > 0:
            # Blend between small-z formula and onset weight
            # Use the smaller of (1/γ) and w0 as the upper bound
            w_sweet = min(w_sweet_base, w0 * 0.8)
        else:
            w_sweet = w_sweet_base

        # Ensure sweet spot is in valid range
        w_sweet = max(self._t + 1, min(w_sweet, w0))

        return int(w_sweet)

    def get_onset_weight(self) -> float:
        """Return the fitted onset weight *w0*.

        Returns:
            The onset weight, or ``t + 1`` if the model has not been
            fitted yet.
        """
        return self._params.get("w0", float(self._t + 1))
