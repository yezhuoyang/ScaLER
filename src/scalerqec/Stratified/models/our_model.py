"""
Our S-curve model implementation (Definition 1 from the paper).

Model:
    f_t[μ,α,β](w) = 0.5 / (1 + exp(-(w-μ)/α + β/√(w-t)))

In linear space (after log-logit transform):
    y(w) = log(0.5/P_w - 1) = -(w-μ)/α + β/√(w-t)
         = a*w + b + c/√(w-t)

Parameter relationships:
    a = -1/α
    b = μ/α
    c = β

    α = -1/a
    μ = -b/a = α*b
    β = c
"""

from typing import Any, Callable, Dict, List, Tuple
import numpy as np
from numpy.typing import NDArray

from scalerqec.Stratified.models.base import ScurveModelBase


class OurScurveModel(ScurveModelBase):
    """
    Our S-curve model from Definition 1 in the paper.

    This model uses the formula:
        P_L(w) = 0.5 / (1 + exp(-(w-μ)/α + β/√(w-t)))

    The model has a pole at w=t (the fault-tolerant threshold), ensuring
    P_L = 0 for w ≤ t.

    Parameters:
        alpha (α): Characteristic decay rate (controls slope)
        mu (μ): Threshold parameter (controls position)
        beta (β): Pole strength (controls curvature near fault-tolerant region)

    Internal representation uses (a, b, c) for numerical stability:
        a = -1/α, b = μ/α, c = β
    """

    def __init__(self, t: int = 0, gamma: float = 0.05):
        super().__init__(t=t, gamma=gamma)
        # Initialize default parameters
        self._params = {
            "alpha": 10.0,
            "mu": 0.0,
            "beta": 1.0,
        }
        # Internal (a, b, c) representation
        self._a: float = -0.1
        self._b: float = 0.0
        self._c: float = 1.0

    @property
    def name(self) -> str:
        return "OurModel"

    @property
    def param_names(self) -> Tuple[str, ...]:
        return ("alpha", "mu", "beta")

    def predict(
        self, w: float | NDArray[np.floating[Any]]
    ) -> float | NDArray[np.floating[Any]]:
        """
        Compute P_L(w) = 0.5 / (1 + exp(-(w-μ)/α + β/√(w-t))).

        For w ≤ t, returns 0 (fault-tolerant region).
        """
        w_arr = np.atleast_1d(np.asarray(w, dtype=float))
        result = np.zeros_like(w_arr)

        # Only compute for w > t
        valid = w_arr > self._t
        if np.any(valid):
            w_valid = w_arr[valid]
            alpha = self._params.get("alpha", 10.0)
            mu = self._params.get("mu", 0.0)
            beta = self._params.get("beta", 1.0)

            # Compute exponent: -(w-μ)/α + β/√(w-t)
            sqrt_term = np.sqrt(w_valid - self._t)
            exponent = -(w_valid - mu) / alpha + beta / sqrt_term

            with np.errstate(over="ignore"):
                result[valid] = 0.5 / (1 + np.exp(exponent))

        # Return scalar if input was scalar
        if np.isscalar(w):
            return float(result[0])
        return result

    def transform(
        self, p_w: float | NDArray[np.floating[Any]]
    ) -> float | NDArray[np.floating[Any]]:
        """
        Transform P_L to y(w) = log(0.5/P_w - 1).
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
        Transform y(w) back to P_L = 0.5 / (1 + exp(y)).
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
        Compute y(w) = a*w + b + c/√(w-t) in transformed linear space.
        """
        w_arr = np.atleast_1d(np.asarray(w, dtype=float))
        result = np.full_like(w_arr, np.inf)

        valid = w_arr > self._t
        if np.any(valid):
            w_valid = w_arr[valid]
            delta = np.sqrt(w_valid - self._t)
            # Avoid division by zero
            delta = np.where(delta < 1e-12, 1e-12, delta)
            result[valid] = self._a * w_valid + self._b + self._c / delta

        if np.isscalar(w):
            return float(result[0])
        return result

    def _get_fit_function(self) -> Callable[..., Any]:
        """
        Return the function for curve fitting: y = a*w + b + c/√(w-t).
        """
        t = self._t

        def modified_linear_function(
            x: float | NDArray[np.floating[Any]], a: float, b: float, c: float
        ) -> float | NDArray[np.floating[Any]]:
            eps = 1e-12
            x_arr = np.atleast_1d(np.asarray(x, dtype=float))
            delta = np.sqrt(np.maximum(x_arr - t, eps))
            delta = np.where(np.abs(delta) < eps, eps, delta)
            result = a * x_arr + b + c / delta
            if np.isscalar(x):
                return float(result[0])
            return result

        return modified_linear_function

    def _get_initial_guess(
        self, x_list: List[float], y_list: List[float]
    ) -> Tuple[float, ...]:
        """
        Get initial parameter guess from linear regression.
        """
        if len(x_list) < 2:
            return (-0.1, 0.0, 1.0)

        # Simple linear regression for initial (a, b) guess
        x0, x1 = x_list[0], x_list[-1]
        y0, y1 = y_list[0], y_list[-1]

        if abs(x1 - x0) < 1e-6:
            a = -0.1
        else:
            a = (y1 - y0) / (x1 - x0)

        # Ensure a is negative (slope should decrease)
        if a >= 0:
            a = -0.1

        b = y0 - a * x0

        # Initial guess for c (beta)
        alpha = -1.0 / a if a != 0 else 10.0
        c = abs(alpha)  # beta ~ alpha is a reasonable starting point

        return (a, b, c)

    def _get_bounds(
        self, initial_guess: Tuple[float, ...]
    ) -> Tuple[List[float], List[float]]:
        """
        Get parameter bounds for fitting.
        """
        a, b, c = initial_guess

        # Allow wide range around initial guess
        lower = [
            min(a * 5.0, a * 0.2),
            min(b * 5.0, b * 0.2) if b != 0 else -100.0,
            min(c * 0.2, c * 5.0) if c != 0 else 0.01,
        ]
        upper = [
            max(a * 5.0, a * 0.2),
            max(b * 5.0, b * 0.2) if b != 0 else 100.0,
            max(c * 0.2, c * 5.0) if c != 0 else 100.0,
        ]

        # Ensure lower < upper
        for i in range(3):
            if lower[i] >= upper[i]:
                lower[i], upper[i] = upper[i], lower[i]
            if lower[i] == upper[i]:
                lower[i] = lower[i] - 1.0
                upper[i] = upper[i] + 1.0

        return lower, upper

    def _set_params_from_fit(self, popt: Tuple[float, ...]) -> None:
        """
        Set model parameters from fitted (a, b, c) values.
        """
        self._a, self._b, self._c = popt[0], popt[1], popt[2]

        # Convert to (alpha, mu, beta)
        if abs(self._a) > 1e-12:
            alpha = -1.0 / self._a
            mu = -self._b / self._a  # mu = alpha * b = (-1/a) * b = -b/a
        else:
            alpha = 10.0
            mu = 0.0

        beta = self._c

        self._params = {
            "alpha": alpha,
            "mu": mu,
            "beta": beta,
        }

    def calculate_sweet_spot(self) -> int:
        """
        Calculate the sweet spot using the general formula: d²y/dw² = γ * dy/dw.

        For our model: y = a*w + b + c/√(w-t)
            dy/dw = a - c/(2*(w-t)^{3/2})
            d²y/dw² = (3/4)*c/(w-t)^{5/2}

        Setting d²y/dw² = γ * |dy/dw| and approximating for large w:
            (3/4)*c/(w-t)^{5/2} ≈ γ * |a|
            (w-t)^{5/2} ≈ (3*c)/(4*γ*|a|)
            w-t ≈ [(3*c)/(4*γ*|a|)]^{2/5}

        Returns:
            Sweet spot weight (integer > t)
        """
        beta = self._c
        a = self._a

        if abs(a) < 1e-12 or beta <= 0 or self._gamma <= 0:
            # Fallback: use expected weight
            return self._t + 1

        # General formula: w = t + [(3*β)/(4*γ*|a|)]^(2/5)
        ratio = (3 * beta) / (4 * self._gamma * abs(a))
        w_sweet = self._t + ratio ** (2 / 5)

        return max(self._t + 1, int(w_sweet))

    def set_params(self, **params: float) -> None:
        """
        Set model parameters (alpha, mu, beta) and sync internal (a, b, c).

        Args:
            **params: Parameter values to set (alpha, mu, beta)
        """
        self._params.update(params)
        # Sync internal parameters from external params
        alpha = self._params.get("alpha", 10.0)
        mu = self._params.get("mu", 0.0)
        beta = self._params.get("beta", 1.0)

        # Convert to internal (a, b, c)
        # a = -1/α, b = μ/α, c = β
        if abs(alpha) > 1e-12:
            self._a = -1.0 / alpha
            self._b = mu / alpha
        else:
            self._a = -0.1
            self._b = 0.0
        self._c = beta

    def get_internal_params(self) -> Dict[str, float]:
        """
        Get the internal (a, b, c) parameters.

        Returns:
            Dictionary with 'a', 'b', 'c' keys
        """
        return {"a": self._a, "b": self._b, "c": self._c}

    def set_internal_params(self, a: float, b: float, c: float) -> None:
        """
        Set internal (a, b, c) parameters directly.

        Also updates the external (alpha, mu, beta) parameters.
        """
        self._a = a
        self._b = b
        self._c = c
        self._set_params_from_fit((a, b, c))
