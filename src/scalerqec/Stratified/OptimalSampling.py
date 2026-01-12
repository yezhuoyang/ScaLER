"""
Optimal sampling strategy based on MSE minimization.

This module implements the theory from the paper for optimizing
the next sampling step by minimizing the expected MSE of the logical error rate estimator.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from scipy.optimize import minimize


class MSEOptimizer:
    """
    Optimizer for choosing the next sampling point to minimize expected MSE.

    Based on the theoretical framework where:
    MSE(P_L_hat) = Var(P_L_hat) + [Bias(P_L_hat)]^2

    The goal is to greedily choose the next subspace that minimizes
    the expected MSE given the budget constraint.
    """

    def __init__(self, num_noise: int, error_rate: float, t: int):
        """
        Initialize the MSE optimizer.

        Args:
            num_noise: Total number of noise locations in the circuit
            error_rate: Physical error rate
            t: Fault tolerance parameter ((d-1)/2)
        """
        self.num_noise = num_noise
        self.error_rate = error_rate
        self.t = t

    def binomial_weight(self, w: int) -> float:
        """
        Calculate the binomial probability of weight w.

        Args:
            w: Pauli weight

        Returns:
            P(w errors occur) = C(n,w) * p^w * (1-p)^(n-w)
        """
        from scipy.special import comb
        n = self.num_noise
        p = self.error_rate
        return comb(n, w, exact=False) * (p ** w) * ((1 - p) ** (n - w))

    def estimate_subspace_variance(self, w: int, P_w: float, N_w: int) -> float:
        """
        Estimate the variance of the logical error rate estimate for subspace w.

        For a subspace with N_w samples and M_w errors:
        Var(P_w_hat) = P_w * (1 - P_w) / N_w

        Args:
            w: Pauli weight
            P_w: True logical error rate at weight w (or current estimate)
            N_w: Number of samples at weight w

        Returns:
            Estimated variance
        """
        if N_w == 0:
            return float('inf')
        return P_w * (1 - P_w) / N_w

    def estimate_parameter_covariance(
        self,
        weights: List[int],
        P_w_values: List[float],
        N_w_values: List[int],
        a: float,
        b: float,
        c: float
    ) -> np.ndarray:
        """
        Estimate the covariance matrix of fitted parameters (α, β, μ).

        This uses the Fisher Information Matrix approach.
        The Fisher Information for parameter θ at subspace w is:
        I(θ; w) = N_w * (∂P_w/∂θ)^T * (∂P_w/∂θ) / Var(P_w)

        Args:
            weights: List of sampled weights
            P_w_values: Estimated logical error rates for each weight
            N_w_values: Number of samples for each weight
            a, b, c: Current parameter estimates

        Returns:
            3x3 covariance matrix for (a, b, c)
        """
        # Initialize Fisher Information Matrix
        fisher = np.zeros((3, 3))

        for w, P_w, N_w in zip(weights, P_w_values, N_w_values):
            if N_w == 0 or w <= self.t:
                continue

            # Compute gradients of P_w with respect to (a, b, c)
            # P_w = 0.5 / (1 + exp(-w/α + μ/α + β/sqrt(w-t)))
            # where α = -1/a, μ = -b/a

            alpha = -1.0 / a if a != 0 else 1.0
            mu = -b / a if a != 0 else 0.0
            sqrt_term = np.sqrt(w - self.t) if w > self.t else 1.0

            exponent = -w / alpha + mu / alpha + c / sqrt_term
            exp_val = np.exp(exponent)
            denominator = (1 + exp_val) ** 2

            if denominator == 0:
                continue

            # Gradients (chain rule applied)
            # ∂P_w/∂a
            dP_da = 0.5 * exp_val / denominator * (w / (alpha ** 2) - mu / (alpha ** 2))

            # ∂P_w/∂b
            dP_db = 0.5 * exp_val / denominator * (-1 / alpha)

            # ∂P_w/∂c
            dP_dc = 0.5 * exp_val / denominator * (-1 / sqrt_term)

            gradient = np.array([dP_da, dP_db, dP_dc])

            # Variance of P_w estimate
            var_P_w = self.estimate_subspace_variance(w, P_w, N_w)

            if var_P_w > 0:
                # Add to Fisher Information
                fisher += N_w * np.outer(gradient, gradient) / var_P_w

        # Covariance is inverse of Fisher Information
        try:
            cov = np.linalg.inv(fisher)
        except np.linalg.LinAlgError:
            # Singular matrix, return large covariance
            cov = np.eye(3) * 1e6

        return cov

    def estimate_LER_variance(
        self,
        param_cov: np.ndarray,
        a: float,
        b: float,
        c: float,
        w_min: int,
        w_max: int
    ) -> float:
        """
        Estimate the variance of the logical error rate estimator.

        Uses first-order Taylor expansion:
        Var(P_L_hat) ≈ (∂P_L/∂θ)^T * Cov(θ) * (∂P_L/∂θ)

        Args:
            param_cov: Covariance matrix of parameters
            a, b, c: Parameter values
            w_min, w_max: Range of weights for integration

        Returns:
            Estimated variance of LER
        """
        # Compute gradient of P_L with respect to parameters
        gradient = np.zeros(3)

        alpha = -1.0 / a if a != 0 else 1.0

        for w in range(max(w_min, self.t + 1), w_max + 1):
            binomial_prob = self.binomial_weight(w)

            if binomial_prob < 1e-10:
                continue

            sqrt_term = np.sqrt(w - self.t)
            mu = -b / a if a != 0 else 0.0
            exponent = -w / alpha + mu / alpha + c / sqrt_term
            exp_val = np.exp(exponent)
            denominator = (1 + exp_val) ** 2

            if denominator == 0:
                continue

            # Gradients weighted by binomial probability
            dP_da = 0.5 * exp_val / denominator * (w / (alpha ** 2) - mu / (alpha ** 2))
            dP_db = 0.5 * exp_val / denominator * (-1 / alpha)
            dP_dc = 0.5 * exp_val / denominator * (-1 / sqrt_term)

            gradient[0] += binomial_prob * dP_da
            gradient[1] += binomial_prob * dP_db
            gradient[2] += binomial_prob * dP_dc

        # Var(P_L) = gradient^T * Cov * gradient
        variance = gradient.T @ param_cov @ gradient

        return max(0, variance)  # Ensure non-negative

    def estimate_LER_bias(
        self,
        weights: List[int],
        N_w_values: List[int],
        a: float,
        b: float,
        c: float
    ) -> float:
        """
        Estimate the bias of the logical error rate estimator.

        The bias comes from:
        1. Finite sample bias in estimating P_w for each subspace
        2. Bias in parameter estimation

        For simplicity, we use a conservative estimate based on
        the variance of parameters and sample sizes.

        Args:
            weights: List of sampled weights
            N_w_values: Number of samples for each weight
            a, b, c: Parameter values

        Returns:
            Estimated bias
        """
        # For large sample sizes, bias is typically small
        # We use a heuristic based on minimum sample size
        min_samples = min(N_w_values) if N_w_values else 0

        if min_samples == 0:
            return 0.1  # High bias when no samples

        # Bias decreases with sample size
        # Heuristic: bias ~ 1/sqrt(min_samples)
        bias = 1.0 / np.sqrt(min_samples) * 0.01

        return bias

    def estimate_MSE(
        self,
        weights: List[int],
        P_w_values: List[float],
        N_w_values: List[int],
        a: float,
        b: float,
        c: float,
        w_min: int,
        w_max: int
    ) -> float:
        """
        Estimate the MSE of the logical error rate estimator.

        MSE = Variance + Bias^2

        Args:
            weights: List of sampled weights
            P_w_values: Logical error rates for each weight
            N_w_values: Number of samples for each weight
            a, b, c: Parameter estimates
            w_min, w_max: Integration range

        Returns:
            Estimated MSE
        """
        # Estimate parameter covariance
        param_cov = self.estimate_parameter_covariance(
            weights, P_w_values, N_w_values, a, b, c
        )

        # Estimate variance
        variance = self.estimate_LER_variance(param_cov, a, b, c, w_min, w_max)

        # Estimate bias
        bias = self.estimate_LER_bias(weights, N_w_values, a, b, c)

        mse = variance + bias ** 2

        return mse

    def choose_next_weight(
        self,
        candidate_weights: List[int],
        current_weights: List[int],
        current_P_w: Dict[int, float],
        current_N_w: Dict[int, int],
        a: float,
        b: float,
        c: float,
        w_min: int,
        w_max: int,
        batch_size: int = 1000
    ) -> int:
        """
        Choose the next weight to sample that minimizes expected MSE.

        For each candidate weight w_new, estimate the MSE if we add
        batch_size samples at w_new, and choose the one with minimum MSE.

        Args:
            candidate_weights: List of candidate weights to consider
            current_weights: Weights already sampled
            current_P_w: Current estimated P_w for each weight
            current_N_w: Current number of samples for each weight
            a, b, c: Current parameter estimates
            w_min, w_max: Integration range
            batch_size: Number of samples to add

        Returns:
            The weight that minimizes expected MSE
        """
        best_weight = candidate_weights[0]
        best_mse = float('inf')

        for w_new in candidate_weights:
            # Create hypothetical sampling scenario
            weights = current_weights.copy()
            P_w_values = [current_P_w.get(w, 0.0) for w in weights]
            N_w_values = [current_N_w.get(w, 0) for w in weights]

            # Add new weight with batch_size samples
            if w_new in current_weights:
                idx = current_weights.index(w_new)
                N_w_values[idx] += batch_size
            else:
                weights.append(w_new)
                # Estimate P_w for new weight using current model
                P_w_new = self._predict_P_w(w_new, a, b, c)
                P_w_values.append(P_w_new)
                N_w_values.append(batch_size)

            # Estimate MSE with this sampling strategy
            mse = self.estimate_MSE(
                weights, P_w_values, N_w_values, a, b, c, w_min, w_max
            )

            if mse < best_mse:
                best_mse = mse
                best_weight = w_new

        return best_weight

    def _predict_P_w(self, w: int, a: float, b: float, c: float) -> float:
        """
        Predict P_w using current model parameters.

        Args:
            w: Weight
            a, b, c: Model parameters

        Returns:
            Predicted logical error rate
        """
        if w <= self.t:
            return 0.0

        alpha = -1.0 / a if a != 0 else 1.0
        mu = -b / a if a != 0 else 0.0
        sqrt_term = np.sqrt(w - self.t)

        exponent = -w / alpha + mu / alpha + c / sqrt_term
        P_w = 0.5 / (1 + np.exp(exponent))

        return np.clip(P_w, 0.0, 0.5)
