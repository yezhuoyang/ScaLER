"""
IBM's S-curve model implementation (Definition 2 from the paper).

Model:
    f[f0, w0, γ](w) = 0.5 * [1 - exp(-2*f0*(w/w0)^γ)]

This is the "min-fail enclosure" model proposed by IBM for approximating
the logical error rate spectrum in quantum error correction.

Parameters:
    f0: Failure probability at the onset weight (controls amplitude)
    w0: Onset weight (minimum weight with logical errors, typically t+1)
    gamma_ibm: Power exponent (controls the shape of the curve)

Note: gamma_ibm is the model parameter, distinct from the sweet spot
tuning parameter gamma in the base class.
"""

from typing import Any, Callable, Dict, List, Tuple
import numpy as np
from numpy.typing import NDArray

from scalerqec.Stratified.models.base import ScurveModelBase


class IBMScurveModel(ScurveModelBase):
    """
    IBM's min-fail enclosure S-curve model from Definition 2 in the paper.

    This model uses the formula:
        P_L(w) = 0.5 * [1 - exp(-2*f0*(w/w0)^γ)]

    Key differences from OurModel:
    - Uses exponential decay instead of logistic
    - Characterized by onset weight w0 rather than pole at t
    - Different asymptotic behavior near saturation

    Parameters:
        f0: Failure probability at onset (controls how fast saturation occurs)
        w0: Onset weight (typically t+1, the first weight with errors)
        gamma_ibm: Power exponent (controls curve steepness)
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

    def predict(self, w: float | NDArray[np.floating[Any]]) -> float | NDArray[np.floating[Any]]:
        """
        Compute P_L(w) = 0.5 * [1 - exp(-2*f0*(w/w0)^γ)].

        For w ≤ 0, returns 0.
        """
        w_arr = np.atleast_1d(np.asarray(w, dtype=float))
        result = np.zeros_like(w_arr)

        f0 = self._params.get("f0", 0.5)
        w0 = self._params.get("w0", float(self._t + 1))
        gamma_ibm = self._params.get("gamma_ibm", 1.0)

        # Avoid division by zero
        if w0 <= 0:
            w0 = 1.0

        # Only compute for w > 0
        valid = w_arr > 0
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

    def transform(self, p_w: float | NDArray[np.floating[Any]]) -> float | NDArray[np.floating[Any]]:
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

    def inverse_transform(self, y: float | NDArray[np.floating[Any]]) -> float | NDArray[np.floating[Any]]:
        """
        Transform y back to P_L.
        """
        y_arr = np.atleast_1d(np.asarray(y, dtype=float))

        with np.errstate(over="ignore"):
            result = 0.5 / (1 + np.exp(y_arr))

        if np.isscalar(y):
            return float(result[0])
        return result

    def linear_prediction(self, w: float | NDArray[np.floating[Any]]) -> float | NDArray[np.floating[Any]]:
        """
        Compute y(w) in transformed linear space for IBM model.

        y(w) = log(0.5/P(w) - 1) = log((1 + exp(-x))/(1 - exp(-x)))
        where x = 2*f0*(w/w0)^γ
        """
        w_arr = np.atleast_1d(np.asarray(w, dtype=float))
        result = np.full_like(w_arr, np.inf)

        f0 = self._params.get("f0", 0.5)
        w0 = self._params.get("w0", float(self._t + 1))
        gamma_ibm = self._params.get("gamma_ibm", 1.0)

        if w0 <= 0:
            w0 = 1.0

        valid = w_arr > 0
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

    def _get_initial_guess(self, x_list: List[float], y_list: List[float]) -> Tuple[float, ...]:
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

    def _get_bounds(self, initial_guess: Tuple[float, ...]) -> Tuple[List[float], List[float]]:
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
        Calculate the sweet spot for IBM model.

        For IBM model: P = 0.5*(1 - exp(-2*f0*(w/w0)^γ))

        The sweet spot is where the curvature is maximal.
        Using d²y/dw² = γ * dy/dw (general definition):

        For this model, we approximate numerically or use a heuristic:
        The inflection point is approximately at:
            w ≈ w0 * (log(2) / (2*f0))^(1/γ)
        """
        f0 = self._params.get("f0", 0.5)
        w0 = self._params.get("w0", float(self._t + 1))
        gamma_ibm = self._params.get("gamma_ibm", 1.0)

        if f0 <= 0 or w0 <= 0 or gamma_ibm <= 0:
            return max(1, self._t + 1)

        # Heuristic: sweet spot is where x = 2*f0*(w/w0)^γ ≈ 0.5
        # => (w/w0)^γ = 0.5/(2*f0) = 0.25/f0
        # => w/w0 = (0.25/f0)^(1/γ)
        ratio = (0.25 / f0) ** (1 / gamma_ibm)
        w_sweet = w0 * ratio

        # Alternatively, using the gamma tuning parameter:
        # Find where second derivative condition is satisfied
        # This is more complex for IBM model, so we use the heuristic above

        return max(1, self._t + 1, int(w_sweet))

    def get_onset_weight(self) -> float:
        """
        Get the onset weight w0.
        """
        return self._params.get("w0", float(self._t + 1))
