"""
Optimized Scaler with MSE-based next step selection.

This extends the base Scaler class with an optimized strategy for choosing
the next sampling point by minimizing the expected MSE of the logical error rate.
"""

from typing import List, Tuple
from scalerqec.Stratified.Scaler import Scaler
from scalerqec.Stratified.OptimalSampling import MSEOptimizer


class ScalerOptimized(Scaler):
    """
    Enhanced Scaler that uses MSE minimization to choose the next sampling point.

    Key differences from the base Scaler:
    1. After the first round of sampling, optimizes point selection based on expected MSE
    2. Chooses the subspace that minimizes MSE(P_L_hat) = Var(P_L_hat) + Bias^2
    3. More adaptive to the specific S-curve shape
    """

    def __init__(self, error_rate: float = 0.0, time_budget: int = 30):
        super().__init__(error_rate, time_budget)
        self._mse_optimizer = None
        self._first_round_complete = False

    def _initialize_mse_optimizer(self):
        """Initialize the MSE optimizer once we have circuit information."""
        if self._mse_optimizer is None and self._num_noise > 0:
            self._mse_optimizer = MSEOptimizer(
                num_noise=self._num_noise,
                error_rate=self._error_rate,
                t=self._t
            )

    def next_step(self) -> Tuple[List[int], List[int]]:
        """
        Decide the next step with MSE-based optimization.

        Strategy:
        - First round: Use uniform sampling like base class to get initial estimates
        - Subsequent rounds: Greedily choose subspaces that minimize expected MSE
        """
        # If we cannot take any more shots, stop
        if self._sampling_rate <= 0.0 or self._remaining_time_budget <= 0.0:
            return [], []

        # Total shots we are allowed to take this step
        remaining_shots = int(self._remaining_time_budget * self._sampling_rate)
        if remaining_shots <= 0:
            return [], []

        # Enforce per-step hard cap for memory control
        step_shots = min(remaining_shots, self._MAX_SHOTS_PER_STEP)

        # Initialize MSE optimizer if needed
        self._initialize_mse_optimizer()

        # Choose candidate weights
        wlist = self._choose_candidate_weights()
        if not wlist:
            return [], []

        # Check if we have enough data to do MSE optimization
        num_sampled_weights = sum(1 for w in wlist if self._subspace_sample_used.get(w, 0) > 0)

        # First round or not enough data: use base class strategy
        if not self._first_round_complete or num_sampled_weights < 3:
            self._first_round_complete = True
            return self._allocate_shots_uniform(wlist, step_shots)

        # MSE-based optimization for subsequent rounds
        return self._allocate_shots_mse_optimal(wlist, step_shots)

    def _allocate_shots_uniform(self, wlist: List[int], step_shots: int) -> Tuple[List[int], List[int]]:
        """
        Allocate shots uniformly (used in first round).
        This is the same as the base class strategy.
        """
        # Priority factors for each weight
        factors = []
        center = self._sweet_spot

        for w in wlist:
            f = 1.0

            # 1) Strongly favour the sweet spot and its immediate neighbours
            if center is not None:
                d = abs(w - center)
                if d <= 1:
                    f *= 6.0
                elif d <= 2:
                    f *= 3.0

            # 2) Boost completely unsampled subspaces
            if self._subspace_sample_used.get(w, 0) == 0:
                f *= 4.0

            factors.append(f)

        total_factor = sum(factors)
        if total_factor <= 0.0:
            shots_per_w = max(1, step_shots // len(wlist))
            return wlist, [shots_per_w] * len(wlist)

        min_shots_per_w = max(500, int(step_shots * 0.02))

        slist: List[int] = []
        remaining = step_shots

        for i, (w, f) in enumerate(zip(wlist, factors)):
            if i == len(wlist) - 1:
                s = max(1, remaining)
            else:
                raw = step_shots * (f / total_factor)
                s = max(min_shots_per_w, int(round(raw)))
                remaining -= s
                if remaining < 0:
                    remaining = 0

            slist.append(s)

        print("  [Uniform allocation] weights:", wlist)
        print("  [Uniform allocation] shots:", slist, " (total =", sum(slist), ")")

        return wlist, slist

    def _allocate_shots_mse_optimal(self, wlist: List[int], step_shots: int) -> Tuple[List[int], List[int]]:
        """
        Allocate shots to minimize expected MSE.

        Strategy:
        1. Evaluate MSE for adding samples to each candidate weight
        2. Greedily choose the weight that gives maximum MSE reduction
        3. Allocate most shots to that weight, distribute rest uniformly
        """
        if self._mse_optimizer is None or self._a == 0.0:
            # Fall back to uniform if optimizer not ready
            return self._allocate_shots_uniform(wlist, step_shots)

        # Prepare current state for MSE optimizer
        current_weights = list(self._estimated_subspaceLER.keys())
        current_P_w = self._estimated_subspaceLER
        current_N_w = self._subspace_sample_used

        # Calculate integration range
        import numpy as np
        sigma = int(np.sqrt(self._error_rate * (1.0 - self._error_rate) * self._num_noise))
        if sigma == 0:
            sigma = 1
        ep = int(self._error_rate * self._num_noise)
        w_min = max(self._t + 1, ep - 5 * sigma)
        w_max = min(self._num_noise, ep + 5 * sigma)

        # Choose optimal weight using MSE minimization
        batch_size = min(5000, step_shots // 3)  # Use a reasonable batch size for evaluation

        try:
            optimal_weight = self._mse_optimizer.choose_next_weight(
                candidate_weights=wlist,
                current_weights=current_weights,
                current_P_w=current_P_w,
                current_N_w=current_N_w,
                a=self._a,
                b=self._b,
                c=self._c,
                w_min=w_min,
                w_max=w_max,
                batch_size=batch_size
            )

            print(f"  [MSE-optimal] Chose weight {optimal_weight} to minimize expected MSE")

            # Allocate shots: give more to optimal weight
            slist = []
            for w in wlist:
                if w == optimal_weight:
                    # Allocate 60% of shots to optimal weight
                    s = int(step_shots * 0.6)
                elif abs(w - optimal_weight) <= 1:
                    # Allocate 15% to neighbors
                    s = int(step_shots * 0.15)
                else:
                    # Distribute rest uniformly
                    remaining_weights = len([x for x in wlist if x != optimal_weight and abs(x - optimal_weight) > 1])
                    if remaining_weights > 0:
                        s = max(500, int(step_shots * 0.25 / remaining_weights))
                    else:
                        s = 500

                slist.append(s)

            # Adjust to match total budget
            total_allocated = sum(slist)
            if total_allocated != step_shots and total_allocated > 0:
                scale = step_shots / total_allocated
                slist = [max(1, int(s * scale)) for s in slist]

                # Fix rounding errors
                diff = step_shots - sum(slist)
                if diff > 0:
                    slist[wlist.index(optimal_weight)] += diff
                elif diff < 0:
                    slist[wlist.index(optimal_weight)] -= abs(diff)

            print("  [MSE-optimal allocation] weights:", wlist)
            print("  [MSE-optimal allocation] shots:", slist, " (total =", sum(slist), ")")

            return wlist, slist

        except Exception as e:
            print(f"  [MSE optimization failed: {e}, falling back to uniform]")
            return self._allocate_shots_uniform(wlist, step_shots)


if __name__ == "__main__":
    # Example usage
    filepath = "C:/Users/yezhu/GitRepos/ScaLERQEC/stimprograms/surface/surface7"
    scaler = ScalerOptimized(error_rate=0.001, time_budget=300)
    scaler.calculate_LER_from_file(
        filepath,
        pvalue=0.001,
        codedistance=7,
        figname="test_optimized_surface7_",
        titlename="Optimized Surface 7",
        repeat=1,
    )
