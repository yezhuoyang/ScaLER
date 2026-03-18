"""Legacy stratified sampling LER estimator with S-curve fitting.

This module provides :class:`StratifiedScurveLERcalc`, which extends
the basic stratified approach by fitting the modified-linear S-curve
model ``y(w) = a*w + b + c/sqrt(w-t)`` and using the fitted curve to
interpolate P_L(w) at un-sampled weights.  This reduces the number of
weight subspaces that must be sampled compared to
:class:`~scalerqec.Stratified.stratifiedLER.StratifiedLERcalc`.

For new code, prefer :class:`~scalerqec.Stratified.Scaler.Scaler` which
adds time-budgeted multi-phase sampling and a model-abstraction layer.
"""

from __future__ import annotations
from typing import Optional
from scalerqec.qepg import (
    compile_QEPG,
    return_samples_many_weights_separate_obs_with_QEPG,
    return_samples_with_fixed_QEPG,
    QEPGGraph,
)
from scalerqec.Clifford.clifford import *
import math
import pymatching
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from contextlib import redirect_stdout
import pickle
import time
from ..QEC.noisemodel import NoiseModel
from ..QEC.qeccircuit import StabCode
from ..util import binomial_weight, format_with_uncertainty
from .ScurveModel import *
from .fitting import r_squared

class StratifiedScurveLERcalc:
    """Stratified sampling LER estimator with S-curve fitting (legacy).

    This class combines stratified sampling with S-curve model fitting
    to estimate the logical error rate.  It samples a small number of
    weight subspaces, fits the modified-linear model
    ``y(w) = a*w + b + c/sqrt(w-t)`` in log-logit space, and uses the
    fitted curve to interpolate P_L(w) at un-sampled weights.

    The workflow is:

    1. Binary-search to find the weight range ``[w_err, w_sat]``.
    2. Sample a uniform grid of points in that range.
    3. Fit a linear model, then the full modified-linear model.
    4. Compute the sweet spot and optionally refine near it.
    5. Integrate to obtain the total LER.

    .. note::
        This is the legacy implementation.  For new code, prefer
        :class:`~scalerqec.Stratified.Scaler.Scaler`, which uses a
        time-budget, model-abstraction layer, and multi-phase algorithm.

    Args:
        error_rate: Physical error rate per noise location.
        sampleBudget: Maximum total number of samples (shots).
        k_range: Number of standard deviations for the critical
            region ``[ep - k*sigma, ep + k*sigma]``.
        num_subspace: Number of evenly-spaced weight subspaces to
            sample within the S-curve range.
        beta: Initial guess for the pole-strength parameter.
    """

    def __init__(
        self,
        error_rate: float = 0.0,
        sampleBudget: int = 10000,
        k_range: int = 3,
        num_subspace: int = 5,
        beta: float = 4,
    ):
        self._num_detector: int = 0
        self._num_noise: int = 0
        self._error_rate: float = error_rate
        self._cliffordcircuit: CliffordCircuit = CliffordCircuit(4)

        self._ler: float = 0
        """
        Use a dictionary to store the estimated subspace logical error rate,
        how many samples have been used in each subspace
        """
        self._estimated_subspaceLER: dict[int, float] = {}
        self._subspace_LE_count: dict[int, int] = {}
        self._estimated_subspaceLER_second: dict[int, float] = {}
        self._subspace_sample_used: dict[int, int] = {}

        self._sampleBudget: int = sampleBudget
        self._sample_used: int = 0
        self._circuit_level_code_distance: int = 1
        self._t: int = 1
        self._num_subspace: int = num_subspace
        """
        minw and maxw store the range of subspace we need to fit.
        This is determined by the uncertainty value
        """
        self._minw: int = 1
        self._maxw: int = 10000000000000
        """
        self._saturatew is the weight of the subspace where the 
        logical error get saturated
        """
        self._saturatew: int = 10000000000000
        self._has_logical_errorw: int = 0
        self._estimated_wlist: list[int] = []

        self._stim_str_after_rewrite: str = ""

        self._mu: float = 0
        self._sigma: float = 0

        # In the area we are interested in, the maximum value of the logical error rate
        self._rough_value_for_subspace_LER: float = 0

        self._stratified_succeed: bool = False

        self._k_range: int = k_range
        self._QEPG_graph: Optional[QEPGGraph] = None

        self._R_square_score: float = 0
        self._beta: float = beta

        self._sweet_spot: float = 0.000

        self._min_num_ke_event: int = 100
        self._sample_gap: int = 100
        self._max_sample_gap: int = 1000000
        self._max_subspace_sample: int = 5000000

        self._ratio: float = 0.05
        self._max_PL: float = 0.15

    def set_sample_bound(
        self,
        MIN_NUM_LE_EVENT: int,
        SAMPLE_GAP: int,
        MAX_SAMPLE_GAP: int,
        MAX_SUBSPACE_SAMPLE: int,
    ):
        """Configure the adaptive sampling thresholds.

        Args:
            MIN_NUM_LE_EVENT: Minimum number of logical error events
                required per subspace before it is considered converged.
            SAMPLE_GAP: Minimum batch size when exploring a subspace
                with no errors yet.
            MAX_SAMPLE_GAP: Maximum batch size per sampling round.
            MAX_SUBSPACE_SAMPLE: Hard cap on total samples per subspace.
        """
        self._min_num_ke_event = MIN_NUM_LE_EVENT
        self._sample_gap = SAMPLE_GAP
        self._max_sample_gap = MAX_SAMPLE_GAP
        self._max_subspace_sample = MAX_SUBSPACE_SAMPLE

    def clear_all(self):
        """Reset all internal state for a fresh estimation run."""
        self._estimated_subspaceLER = {}
        self._subspace_LE_count = {}
        self._estimated_subspaceLER_second = {}
        self._subspace_sample_used = {}
        self._sample_used = 0
        self._ler = 0
        self._estimated_wlist = []
        self._saturatew = 10000000000000
        self._minw = self._t + 1
        self._maxw = 10000000000000
        # self._cliffordcircuit=CliffordCircuit(4)
        self._R_square_score = 0

    def calc_logical_error_rate_with_fixed_w(self, shots: int, w: int):
        """Estimate P_L(w) by sampling *shots* error patterns of weight *w*.

        Generates *shots* random fault configurations with exactly *w*
        faults, decodes each one, and returns the fraction that produce
        a logical error.

        Args:
            shots: Number of samples to draw.
            w: Exact Pauli fault weight.

        Returns:
            Estimated conditional logical error probability P_L(w).
        """
        assert self._QEPG_graph is not None, (
            "QEPG graph must be initialized before sampling"
        )
        result = return_samples_with_fixed_QEPG(self._QEPG_graph, w, shots)
        # self._sample_used+=shots
        # if w not in self._subspace_LE_count.keys():
        #     self._subspace_LE_count[w]=0
        #     self._subspace_sample_used[w]=shots
        # else:
        #     self._subspace_sample_used[w]+=shots
        arr = np.asarray(result)
        states = arr[:, :-1]
        observables = arr[:, -1]
        predictions = np.squeeze(self._matcher.decode_batch(states))
        num_errors = np.count_nonzero(observables != predictions)
        # self._subspace_LE_count[w]+=num_errors
        return num_errors / shots

    def binary_search_upper(self, low: int, high: int, shots: int):
        """Find the smallest weight where P_L exceeds the saturation threshold.

        Uses binary search over ``[low, high]`` to find the first
        weight *w* for which ``P_L(w) > _max_PL``.

        Args:
            low: Lower bound of the search range.
            high: Upper bound of the search range.
            shots: Number of samples per probe.

        Returns:
            The smallest weight at which P_L exceeds ``_max_PL``.
        """
        left: int = low
        right: int = high
        epsion: float = self._max_PL
        while left < right:
            mid = (left + right) // 2
            er = self.calc_logical_error_rate_with_fixed_w(shots, mid)
            if er > epsion:
                right = mid
            else:
                left = mid + 1
        return left

    def binary_search_lower(self, low: int, high: int, shots: int = 5000):
        """Find the smallest weight where P_L is noticeably non-zero.

        Uses binary search over ``[low, high]`` to find the first
        weight *w* for which ``P_L(w) > 0.002``.

        Args:
            low: Lower bound of the search range.
            high: Upper bound of the search range.
            shots: Number of samples per probe.

        Returns:
            The smallest weight at which P_L is detectable.
        """
        left: int = low
        right: int = high
        epsion: float = 0.002
        while left < right:
            mid = (left + right) // 2
            er = self.calc_logical_error_rate_with_fixed_w(shots, mid)
            if er > epsion:
                right = mid
            else:
                left = mid + 1
        return left

    def determine_lower_w(self):
        """Determine the first weight with non-zero logical errors.

        Sets ``_has_logical_errorw`` to the smallest weight at which
        P_L is detectable. For very small circuits (N <= 8) defaults
        to 1.
        """
        if self._num_noise <= 8:
            self._has_logical_errorw = 1
        else:
            self._has_logical_errorw = self.binary_search_lower(
                self._t + 1, self._num_noise
            )
        # self._has_logical_errorw=self._t+100

    def determine_saturated_w(self, shots: int = 1000):
        """Determine the weight at which P_L is essentially saturated.

        Sets ``_saturatew`` to the smallest weight for which
        ``P_L(w) > _max_PL``, indicating the S-curve has reached
        its plateau.

        Args:
            shots: Number of samples per binary-search probe.
        """
        # self._saturatew=self._num_detector//30
        if self._num_noise <= 8:
            self._saturatew = self._num_noise
        else:
            self._saturatew = self.binary_search_upper(
                self._minw, self._num_noise, shots
            )
            if self._saturatew < self._has_logical_errorw + 8:
                self._saturatew = self._has_logical_errorw + 8
        # print("Self._saturatew: ",self._saturatew)

    def parse_from_file(self, filepath: str):
        """Load a STIM circuit from *filepath* and prepare for sampling.

        Compiles the circuit into a Clifford representation, builds the
        QEPG graph for weight-stratified sampling, and initialises a
        PyMatching decoder from the detector error model.

        Args:
            filepath: Path to a STIM circuit file.
        """
        stim_str = ""
        with open(filepath, "r", encoding="utf-8") as f:
            stim_str = f.read()

        self._cliffordcircuit.error_rate = self._error_rate
        self._cliffordcircuit.compile_from_stim_circuit_str(stim_str)
        self._num_noise = self._cliffordcircuit.totalnoise
        self._num_detector = len(self._cliffordcircuit.parityMatchGroup)
        self._stim_str_after_rewrite = stim_str

        # Configure a decoder using the circuit.
        self._detector_error_model = (
            self._cliffordcircuit.stimcircuit.detector_error_model(
                decompose_errors=True
            )
        )
        self._matcher = pymatching.Matching.from_detector_error_model(
            self._detector_error_model
        )

        self._QEPG_graph = compile_QEPG(stim_str)

    def determine_range_to_sample(self):
        """Determine the critical weight range [minw, maxw].

        Computes the range of weights whose binomial contribution to
        the LER is significant.  The range is centred on the expected
        number of faults ``N * p`` and extends +/- ``k_range`` standard
        deviations of the binomial distribution.
        """
        sigma = int(
            np.sqrt(self._error_rate * (1 - self._error_rate) * self._num_noise)
        )
        if sigma == 0:
            sigma = 1
        ep = int(self._error_rate * self._num_noise)
        self._minw = max(self._t + 1, ep - self._k_range * sigma)
        self._maxw = max(self._num_subspace, ep + self._k_range * sigma)
        self._maxw = min(self._maxw, self._num_noise)

    def subspace_sampling(self):
        """Adaptively sample weight subspaces near the sweet spot.

        Draws samples at evenly-spaced weights between the sweet spot
        and the first weight with logical errors. Sampling at each
        weight continues until either the minimum logical error event
        count is reached or the total sample budget is exhausted.
        """
        assert self._QEPG_graph is not None, (
            "QEPG graph must be initialized before sampling"
        )
        # wlist_need_to_sample = list(range(self._minw, self._maxw + 1))
        # wlist_need_to_sample=evenly_spaced_ints(self._sweet_spot,self._saturatew,self._num_subspace)

        wlist_need_to_sample = evenly_spaced_ints(
            int(self._sweet_spot), self._has_logical_errorw, self._num_subspace
        )

        # print("wlist_need_to_sample: ",wlist_need_to_sample)
        for weight in wlist_need_to_sample:
            if not weight in self._estimated_wlist:
                self._estimated_wlist.append(weight)
                self._subspace_LE_count[weight] = 0
                self._subspace_sample_used[weight] = 0

        # print(wlist_need_to_sample)
        self._sample_used = 0
        total_LE_count = 0
        while True:
            x_list = [
                x
                for x in self._estimated_subspaceLER.keys()
                if (
                    self._estimated_subspaceLER[x] < 0.5
                    and self._estimated_subspaceLER[x] > 0
                )
            ]

            slist = []
            wlist = []
            """
            Case 1 to end the while loop: We have consumed all of our sample budgets
            """
            if self._sample_used > self._sampleBudget:
                break

            for weight in wlist_need_to_sample:
                """
                When we declare the circuit level code distance, we don't need to sample these subspaces
                """
                if self._subspace_sample_used[weight] > self._max_subspace_sample:
                    continue

                if self._subspace_LE_count[weight] < self._min_num_ke_event:
                    if self._subspace_LE_count[weight] >= 1:
                        """
                        For larger subspaces, when we have already get some logical error, 
                        we can estimate how many we still need to sample
                        """
                        sample_num_required = (
                            int(
                                self._min_num_ke_event / self._subspace_LE_count[weight]
                            )
                            * self._subspace_sample_used[weight]
                        )
                        if sample_num_required > self._max_sample_gap:
                            sample_num_required = self._max_sample_gap
                        slist.append(sample_num_required)
                        self._subspace_sample_used[weight] += sample_num_required
                        self._sample_used += sample_num_required
                    else:
                        """
                        For larger subspaces, if we have not get any logical error, then we double the sample size
                        """
                        sample_num_required = max(
                            self._sample_gap, self._subspace_sample_used[weight] * 10
                        )
                        if sample_num_required > self._max_sample_gap:
                            sample_num_required = self._max_sample_gap
                        slist.append(sample_num_required)
                        self._subspace_sample_used[weight] += sample_num_required
                        self._sample_used += sample_num_required
                    wlist.append(weight)
            """
            Case 2 to end the while loop: We have get 100 logical error events for all these subspaces
            """
            if len(wlist) == 0:
                break

            # print("wlist: ",wlist)
            # print("slist: ",slist)
            # detector_result,obsresult=return_samples_many_weights_separate_obs(self._stim_str_after_rewrite,wlist,slist)
            detector_result, obsresult = (
                return_samples_many_weights_separate_obs_with_QEPG(
                    self._QEPG_graph, wlist, slist
                )
            )
            predictions_result = self._matcher.decode_batch(detector_result)
            # print("Result get!")

            begin_index = 0
            for w_idx, (w, quota) in enumerate(zip(wlist, slist)):
                observables = np.asarray(
                    obsresult[begin_index : begin_index + quota]
                )  # (shots,)
                predictions = np.asarray(
                    predictions_result[begin_index : begin_index + quota]
                ).ravel()

                # 3. count mismatches in vectorised form ---------------------------------
                num_errors = np.count_nonzero(observables != predictions)
                total_LE_count += num_errors
                self._subspace_LE_count[w] += num_errors
                self._estimated_subspaceLER[w] = (
                    self._subspace_LE_count[w] / self._subspace_sample_used[w]
                )

                # print(f"Logical error rate when w={w}: {self._estimated_subspaceLER[w]*binomial_weight(self._num_noise, w,self._error_rate):.6g}")

                begin_index += quota
            # print(self._subspace_LE_count)
            # print("Subspace LE count: ",self._subspace_LE_count)
            # print("self._subspace_sample_used: ",self._subspace_sample_used)

        # print("Samples used:{}".format(self._sample_used))
        # print("circuit level code distance:{}".format(self._circuit_level_code_distance))
        # print(self._subspace_LE_count)

    def subspace_sampling_to_fit_curve(self, sampleBudget):
        """Sample evenly-spaced subspaces for S-curve fitting.

        Distributes *sampleBudget* shots evenly across
        ``num_subspace`` weights between ``w_err`` and ``w_sat``. The
        resulting P_L(w) estimates are used in the subsequent curve
        fitting step.

        Args:
            sampleBudget: Total number of shots to distribute.
        """
        assert self._QEPG_graph is not None, (
            "QEPG graph must be initialized before sampling"
        )

        wlist = evenly_spaced_ints(
            self._has_logical_errorw, self._saturatew, self._num_subspace
        )
        for weight in wlist:
            if not (weight in self._estimated_wlist):
                self._estimated_wlist.append(weight)
        slist = [sampleBudget // self._num_subspace] * len(wlist)

        # detector_result,obsresult=return_samples_many_weights_separate_obs(self._stim_str_after_rewrite,wlist,slist)
        detector_result, obsresult = return_samples_many_weights_separate_obs_with_QEPG(
            self._QEPG_graph, wlist, slist
        )
        predictions_result = self._matcher.decode_batch(detector_result)

        for w, s in zip(wlist, slist):
            if not w in self._subspace_LE_count.keys():
                self._subspace_LE_count[w] = 0
                self._subspace_sample_used[w] = s
                self._estimated_subspaceLER[w] = 0
            else:
                self._subspace_sample_used[w] += s

        begin_index = 0
        for w_idx, (w, quota) in enumerate(zip(wlist, slist)):
            observables = np.asarray(
                obsresult[begin_index : begin_index + quota]
            )  # (shots,)
            predictions = np.asarray(
                predictions_result[begin_index : begin_index + quota]
            ).ravel()

            # 3. count mismatches in vectorised form ---------------------------------
            num_errors = np.count_nonzero(observables != predictions)

            self._subspace_LE_count[w] += num_errors
            self._estimated_subspaceLER[w] = (
                self._subspace_LE_count[w] / self._subspace_sample_used[w]
            )
            # print("Logical error rate when w={} ".format(w)+str(self._estimated_subspaceLER[w]))
            begin_index += quota

    def calculate_R_square_score(self):
        """Compute R^2 for the current S-curve fit against measured data.

        Returns:
            Coefficient of determination (R^2).
        """
        y_observed = [self._estimated_subspaceLER[x] for x in self._estimated_wlist]
        y_predicted = [
            scurve_function(x, self._mu, self._sigma) for x in self._estimated_wlist
        ]
        # y_predicted = [scurve_function_with_distance(x,self._mu,self._sigma) for x in self._estimated_wlist]
        r2 = r_squared(y_observed, y_predicted)
        # print("R^2 score: ", r2)
        return r2

    def calculate_LER(self):
        """Compute the total LER from measured subspace P_L values.

        Sums ``P_L(w) * Binom(N, w) * p^w * (1-p)^(N-w)`` over
        sampled weights only (no interpolation).

        Returns:
            The estimated total logical error rate.
        """
        self._ler = 0
        for weight in range(1, self._num_noise + 1):
            if weight in self._estimated_subspaceLER.keys():
                self._ler += self._estimated_subspaceLER[weight] * binomial_weight(
                    self._num_noise, weight, self._error_rate
                )
        return self._ler

    def get_LER_subspace(self, weight):
        """Return the LER contribution of a single weight subspace.

        Args:
            weight: The fault weight.

        Returns:
            ``P_L(w) * Binom(N, w) * p^w * (1-p)^(N-w)``.
        """
        return self._estimated_subspaceLER[weight] * binomial_weight(
            self._num_noise, weight, self._error_rate
        )

    def fit_linear_area(self):
        """Fit a simple linear model ``y = a*w + b`` in log-logit space.

        Performs a weighted least-squares fit using bias-corrected
        log-logit values and sigma estimates as weights.  Updates
        ``_a``, ``_b``, ``_minw``, ``_maxw``, and saves a diagnostic
        plot to ``total_linear_fit.pdf``.
        """
        x_list = [
            x
            for x in self._estimated_subspaceLER.keys()
            if (
                self._estimated_subspaceLER[x] < 0.5
                and self._estimated_subspaceLER[x] > 0
            )
        ]

        # y_list = [np.log(0.5/self._estimated_subspaceLER[x]-1) for x in x_list]
        y_list = [
            np.log(0.5 / self._estimated_subspaceLER[x] - 1)
            - bias_estimator(self._subspace_sample_used[x], self._subspace_LE_count[x])
            for x in x_list
        ]
        sigma_list = [
            sigma_estimator(self._subspace_sample_used[x], self._subspace_LE_count[x])
            for x in x_list
        ]

        # print(x_list)
        # print(y_list)
        popt, pcov = curve_fit(
            linear_function,
            x_list,
            y_list,
            sigma=sigma_list,  # <-- use sigma_list as weights
        )

        sigma = int(
            np.sqrt(self._error_rate * (1 - self._error_rate) * self._num_noise)
        )
        if sigma == 0:
            sigma = 1
        ep = int(self._error_rate * self._num_noise)
        self._minw = max(self._t + 1, ep - self._k_range * sigma)
        self._maxw = max(2, ep + self._k_range * sigma)

        self._a, self._b = popt[0], popt[1]

        # Plot the fitted line
        x_fit = np.linspace(min(x_list), max(x_list), 1000)

        y_fit = linear_function(x_fit, self._a, self._b)

        # print("Fitted parameters: a={}, b={}".format(self._a,self._b))
        plt.figure()
        plt.plot(x_fit, y_fit, label="Fitted line", color="orange")
        plt.scatter(x_list, y_list, label="Data points", color="blue")

        # ── NEW: highlight linear-fit window ───────────────────────────────────
        plt.axvline(
            self._minw, color="red", linestyle="--", linewidth=1.2, label=r"$w_{\min}$"
        )
        plt.axvline(
            self._maxw,
            color="green",
            linestyle="--",
            linewidth=1.2,
            label=r"$w_{\max}$",
        )
        plt.axvspan(
            self._minw, self._maxw, color="gray", alpha=0.15
        )  # translucent fill
        # ───────────────────────────────────────────────────────────────────────

        alpha = -1 / self._a
        mu = alpha * self._b

        # ── NEW: annotate alpha and mu ─────────────────────────────────────
        textstr = "\n".join((r"$\alpha=%.4f$" % alpha, r"$\mu=%.4f$" % mu))
        plt.text(
            0.05,
            0.95,
            textstr,
            transform=plt.gca().transAxes,
            fontsize=10,
            verticalalignment="top",
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="black"),
        )
        # ───────────────────────────────────────────────────────────────────

        plt.xlabel("Weight")
        plt.ylabel("Linear")
        plt.title("Linear Fit of S-curve")
        plt.legend()
        plt.savefig("total_linear_fit.pdf", dpi=300)
        plt.close()

    def fit_log_S_model(self, filename, savefigure=True, time=None):
        """Fit the modified-linear S-curve model and update the LER.

        Performs a weighted nonlinear least-squares fit of
        ``y(w) = a*w + b + c/sqrt(w-t)`` in log-logit space, computes
        the R^2 score, calculates the sweet spot, and optionally saves
        a diagnostic plot.

        Args:
            filename: Output path for the diagnostic PDF plot.
            savefigure: If ``True``, save the plot to *filename*.
            time: Elapsed wall-clock time (seconds) to annotate on the
                plot. Defaults to ``None``.
        """
        x_list = [
            x
            for x in self._estimated_subspaceLER.keys()
            if (
                self._estimated_subspaceLER[x] < 0.5
                and self._estimated_subspaceLER[x] > 0
                and self._subspace_LE_count[x] >= (self._min_num_ke_event // 5)
            )
        ]

        sigma_list = [
            sigma_estimator(self._subspace_sample_used[x], self._subspace_LE_count[x])
            for x in x_list
        ]
        y_list = [
            np.log(0.5 / self._estimated_subspaceLER[x] - 1)
            - bias_estimator(self._subspace_sample_used[x], self._subspace_LE_count[x])
            for x in x_list
        ]

        # print("Saturated weight: ",self._saturatew)
        # print("LE count: ",self._subspace_LE_count)
        # print("Sample used: ",self._subspace_sample_used)

        non_zero_indices = [x for x in x_list if self._estimated_subspaceLER[x] > 0]
        upper_bound_code_distance = (
            min(non_zero_indices)
            if len(non_zero_indices) > 0
            else self._circuit_level_code_distance * 10
        )

        center = self._saturatew / 2
        sigma = self._saturatew / 7  # center in the middle of that span
        b = self._b
        a = self._a
        c = self._beta
        alpha = -1 / self._a
        initial_guess = (a, b, alpha)

        # print("Initial guess d: ",int((self._circuit_level_code_distance+upper_bound_code_distance)/2))
        sigma = int(
            np.sqrt(self._error_rate * (1 - self._error_rate) * self._num_noise)
        )
        if sigma == 0:
            sigma = 1
        ep = int(self._error_rate * self._num_noise)
        self._minw = max(self._t + 1, ep - self._k_range * sigma)
        self._maxw = max(2, ep + self._k_range * sigma)

        self._num_detector
        self._num_noise

        beta = alpha
        initial_guess = (a, b, beta)
        # ── lower bounds for [param1, param2, param3, param4]

        lower = [
            min(self._a * 5, self._a * 0.2),
            min(self._b * 0.2, self._b * 5),
            min(beta * 0.2, beta * 5),
        ]

        # ── upper bounds for [param1, param2, param3, param4]
        upper = [
            max(self._a * 5, self._a * 0.2),
            max(self._b * 0.2, self._b * 5),
            max(beta * 0.2, beta * 5),
        ]

        popt, pcov = curve_fit(
            modified_linear_function(self._t),
            x_list,
            y_list,
            p0=initial_guess,  # len(initial_guess) must be 4 and within the bounds above
            bounds=(lower, upper),  # <-- tuple with two arrays
            maxfev=50_000,  # or max_nfev in newer SciPy
        )

        self._codedistance = 0
        # Extract the best-fit parameter (alpha)
        self._a, self._b, self._c = popt[0], popt[1], popt[2]

        # print("circuit d:",self._circuit_level_code_distance)
        y_list = [
            np.log(0.5 / self._estimated_subspaceLER[x] - 1)
            - bias_estimator(self._subspace_sample_used[x], self._subspace_LE_count[x])
            for x in x_list
        ]
        y_predicted = [
            modified_linear_function_with_d(x, self._a, self._b, self._c, self._t)
            for x in x_list
        ]
        # y_predicted = [scurve_function_with_distance(x,self._mu,self._sigma) for x in self._estimated_wlist]
        self._R_square_score = r_squared(y_list, y_predicted)
        # print("R^2 score: ", self._R_square_score)

        # Plot the fitted line
        x_fit = np.linspace(self._t + 1, max(x_list), 1000)

        y_fit = modified_linear_function_with_d(
            x_fit, self._a, self._b, self._c, self._t
        )

        self.calc_logical_error_rate_after_curve_fitting()

        alpha = -1 / self._a
        self._sweet_spot = int(
            refined_sweet_spot(alpha, self._c, self._t, ratio=self._ratio)
        )
        if self._sweet_spot < ep:
            self._sweet_spot = ep
        if self._sweet_spot <= self._t:
            self._sweet_spot = self._t + 1
        # self._sweet_spot=self._t+100

        sweet_spot_y = modified_linear_function_with_d(
            self._sweet_spot, self._a, self._b, self._c, self._t
        )

        sample_cost_list = [self._subspace_sample_used[x] for x in x_list]

        # print("Fitted parameters: a={}, b={}, c={}, d={}".format(self._a, self._b, self._c, self._d))

        # Setup the plot
        fig, ax = plt.subplots(figsize=(10, 6))

        # Calculate adaptive ranges for positioning
        x_min_data, x_max_data = min(x_list), max(x_list)
        y_min_data, y_max_data = min(y_list), max(y_list)
        y_range = y_max_data - y_min_data if y_max_data != y_min_data else 1.0
        x_range = x_max_data - x_min_data if x_max_data != x_min_data else 1.0

        # Adaptive saturation region width (proportional to data range)
        saturation_width = max(3, x_range * 0.3)
        x_plot_max = self._saturatew + saturation_width

        # Plot histogram-style bars for the y values
        bar_width = max(0.4, min(0.8, x_range / len(x_list) * 0.6))
        bar_container = ax.bar(
            x_list,
            y_list,
            width=bar_width,
            align="center",
            color="orange",
            edgecolor="orange",
            label="Data histogram",
        )

        # Overlay error bars on top of bars
        ax.errorbar(
            x_list,
            y_list,
            yerr=sigma_list,
            fmt="o",
            color="black",
            capsize=3,
            markersize=1,
            elinewidth=1,
            label="Error bars",
        )

        # Fit curve
        ax.plot(
            x_fit,
            y_fit,
            label=f"Fitted line, R2={self._R_square_score:.4f}",
            color="blue",
            linestyle="--",
        )

        # sweet spot marker
        ax.scatter(
            self._sweet_spot,
            sweet_spot_y,
            color="purple",
            marker="o",
            s=50,
            label="Sweet Spot",
        )

        # Set axis limits before adding regions and text
        y_plot_min = y_min_data - y_range * 0.1
        y_plot_max = y_max_data + y_range * 0.5  # Extra space at top for labels
        ax.set_xlim(-0.5, x_plot_max)
        ax.set_ylim(y_plot_min, y_plot_max)

        # Helper function to position text in data coordinates based on region
        def region_text_y():
            return y_max_data + y_range * 0.35

        # Region: Fault-tolerant (green)
        if self._t > 0:
            ax.axvspan(0, self._t, color="green", alpha=0.15)
            ax.text(
                self._t / 2,
                region_text_y(),
                "Fault\ntolerant",
                ha="center",
                va="bottom",
                color="green",
                fontsize=8,
            )

        # Region: Curve fitting (yellow)
        ax.axvspan(self._t, self._saturatew, color="yellow", alpha=0.10)
        ax.text(
            (self._t + self._saturatew) / 2,
            region_text_y(),
            "Curve fitting",
            ha="center",
            va="bottom",
            fontsize=10,
        )

        # Region: Critical area (gray)
        ax.axvspan(self._minw, self._maxw, color="gray", alpha=0.2)
        ax.axvline(
            self._minw, color="red", linestyle="--", linewidth=1.2, label=r"$w_{\min}$"
        )
        ax.axvline(
            self._maxw,
            color="green",
            linestyle="--",
            linewidth=1.2,
            label=r"$w_{\max}$",
        )
        ax.text(
            (self._minw + self._maxw) / 2,
            region_text_y(),
            r"Critical Region",
            ha="center",
            va="bottom",
            fontsize=9,
        )

        # Region: Saturation (red)
        ax.axvspan(self._saturatew, x_plot_max, color="red", alpha=0.15)
        ax.text(
            (self._saturatew + x_plot_max) / 2,
            region_text_y(),
            "Saturation",
            ha="center",
            va="bottom",
            color="red",
            fontsize=10,
        )

        # Sample cost annotations (scientific notation) - position above bars
        num_points_to_annotate = min(5, len(x_list))
        if num_points_to_annotate > 0:
            indices = np.linspace(
                0, len(x_list) - 1, num=num_points_to_annotate, dtype=int
            )
            for i in indices:
                x, y, s = x_list[i], y_list[i], sample_cost_list[i]
                if s > 0:
                    s_str = "{0:.1e}".format(s)
                    base, exp = s_str.split("e")
                    label = r"${0}\times 10^{{{1}}}$".format(base, int(exp))
                    ax.annotate(
                        label,
                        (x, y),
                        textcoords="offset points",
                        xytext=(0, 8),
                        ha="center",
                        fontsize=7,
                    )

        # Side annotation box
        text_lines = [
            r"$N_{LE}^{Clip}=%d$" % self._min_num_ke_event,
            r"$N_{sub}^{Gap}=%d$" % self._max_sample_gap,
            r"$N_{sub}^{Max}=%d$" % self._max_subspace_sample,
            r"$N_{total}=%d$" % self._sample_used,
            r"$r_{sweet}=%.2f$" % self._ratio,
            r"$\alpha=%.4f$" % alpha,
            r"$\mu =%.4f$" % (alpha * self._b),
            r"$\beta=%.4f$" % self._c,
            r"$w_{\min}=%d$" % self._minw,
            r"$w_{\max}=%d$" % self._maxw,
            r"$w_{sweet}=%d$" % self._sweet_spot,
            r"$\#\mathrm{detector}=%d$" % self._num_detector,
            r"$\#\mathrm{noise}=%d$" % self._num_noise,
            r"$P_L=%.2e$" % self._ler,
        ]
        if time is not None:
            text_lines.append(r"$\mathrm{Time}=%.2f\,\mathrm{s}$" % time)

        # Place annotation box in upper right corner using axes coordinates
        ax.text(
            0.98,
            0.98,
            "\n".join(text_lines),
            transform=ax.transAxes,
            fontsize=7,
            va="top",
            ha="right",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.95),
        )

        # Final formatting
        ax.set_xlabel("Weight")
        ax.set_ylabel(r"$\log\left(\frac{0.5}{\mathrm{LER}} - 1\right)$")
        ax.set_title("Fitted log-S-curve")
        ax.legend(fontsize=8, loc="upper left")
        fig.tight_layout()
        if savefigure:
            fig.savefig(filename, format="pdf", bbox_inches="tight")
        plt.show()
        plt.close()

    def fit_Scurve(self):
        """Fit a basic logistic S-curve ``P_L = 0.5/(1+exp(-(w-mu)/sigma))``.

        Uses ``scipy.optimize.curve_fit`` in probability space (not
        log-logit).  Updates ``_mu`` and ``_sigma``.

        Returns:
            Tuple of ``(code_distance, mu, sigma)``.
        """
        if self._stratified_succeed:
            self._saturatew = self._maxw

        center = self._saturatew / 2
        sigma = self._saturatew / 7  # centre in the middle of that span
        initial_guess = (center, sigma)

        popt, pcov = curve_fit(
            scurve_function,
            self._estimated_wlist,
            [self._estimated_subspaceLER[x] for x in self._estimated_wlist],
            p0=initial_guess,
        )

        self._codedistance = 0
        # Extract the best-fit parameter (alpha)
        self._mu, self._sigma = popt[0], popt[1]
        return self._codedistance, self._mu, self._sigma

    def ground_truth_subspace_sampling(self):
        """Sample all subspaces exhaustively for ground-truth validation.

        Uses the same adaptive sampling logic as :meth:`subspace_sampling`
        but covers the full critical weight range ``[minw, maxw]``.
        Results are stored in ``_ground_*`` attributes for comparison
        against the curve-fitted estimates.
        """
        assert self._QEPG_graph is not None, (
            "QEPG graph must be initialized before sampling"
        )
        sigma = int(
            np.sqrt(self._error_rate * (1 - self._error_rate) * self._num_noise)
        )
        if sigma == 0:
            sigma = 1
        ep = int(self._error_rate * self._num_noise)
        minw = max(self._t + 1, ep - self._k_range * sigma)
        maxw = max(self._num_subspace, ep + self._k_range * sigma)
        maxw = min(maxw, self._num_noise)
        wlist_need_to_sample = list(range(minw, maxw + 1))
        self._ground_sample_used = 0
        self._ground_estimated_subspaceLER = {}
        self._ground_subspace_LE_count = {}
        self._ground_subspace_sample_used = {}
        for weight in wlist_need_to_sample:
            self._ground_subspace_LE_count[weight] = 0
            self._ground_subspace_sample_used[weight] = 0

        while True:
            slist = []
            wlist = []

            if self._ground_sample_used > self._sampleBudget:
                break

            for weight in wlist_need_to_sample:
                if (
                    self._ground_subspace_sample_used[weight]
                    > self._max_subspace_sample
                ):
                    continue
                if self._ground_subspace_LE_count[weight] == 0:
                    wlist.append(weight)
                    sample_num_required = min(
                        self._max_sample_gap,
                        self._ground_subspace_sample_used[weight] * 10,
                    )
                    self._ground_subspace_sample_used[weight] += sample_num_required
                    self._ground_sample_used += sample_num_required
                    slist.append(sample_num_required)
                    continue
                if self._ground_subspace_LE_count[weight] < self._min_num_ke_event:
                    sample_num_required = (
                        int(
                            self._min_num_ke_event
                            / self._ground_subspace_LE_count[weight]
                        )
                        * self._ground_subspace_sample_used[weight]
                    )
                    if sample_num_required > self._max_sample_gap:
                        sample_num_required = self._max_sample_gap
                    self._ground_subspace_sample_used[weight] += sample_num_required
                    self._ground_sample_used += sample_num_required
                    wlist.append(weight)
                    slist.append(sample_num_required)

            if len(wlist) == 0:
                break
            # detector_result,obsresult=return_samples_many_weights_separate_obs(self._stim_str_after_rewrite,wlist,slist)

            # print("Ground truth wlist: ",wlist)
            # print("Ground truth slist: ",slist)

            detector_result, obsresult = (
                return_samples_many_weights_separate_obs_with_QEPG(
                    self._QEPG_graph, wlist, slist
                )
            )
            predictions_result = self._matcher.decode_batch(detector_result)

            begin_index = 0
            for w_idx, (w, quota) in enumerate(zip(wlist, slist)):
                observables = np.asarray(
                    obsresult[begin_index : begin_index + quota]
                )  # (shots,)
                predictions = np.asarray(
                    predictions_result[begin_index : begin_index + quota]
                ).ravel()

                # 3. count mismatches in vectorised form ---------------------------------
                num_errors = np.count_nonzero(observables != predictions)

                self._ground_subspace_LE_count[w] += num_errors
                self._ground_estimated_subspaceLER[w] = (
                    self._ground_subspace_LE_count[w]
                    / self._ground_subspace_sample_used[w]
                )

                # print(f"Logical error rate when w={w}: {self._ground_estimated_subspaceLER[w]*binomial_weight(self._num_noise, w,self._error_rate):.6g}")

                begin_index += quota
            # print(self._ground_subspace_LE_count)
            # print(self._ground_subspace_sample_used)
        # print("Samples used:{}".format(self._ground_sample_used))

    def calc_logical_error_rate_after_curve_fitting(self):
        """Compute the total LER using fitted S-curve for un-sampled weights.

        For weights where direct measurements exist, the measured
        P_L(w) is used.  For un-sampled weights in the critical range,
        the fitted modified-sigmoid function is used to interpolate.

        Returns:
            The estimated total logical error rate.
        """
        self._ler = 0

        sigma = int(
            np.sqrt(self._error_rate * (1 - self._error_rate) * self._num_noise)
        )
        if sigma == 0:
            sigma = 1
        ep = int(self._error_rate * self._num_noise)
        self._minw = max(self._t + 1, ep - self._k_range * sigma)
        self._maxw = max(2, ep + self._k_range * sigma)

        for weight in range(self._minw, self._maxw + 1):
            """
            If the weight is in the estimated list, we use the estimated value
            Else, we use the curve fitting value

            If the weight is less than the minw, we just declare it as 0
            """
            if weight in self._estimated_subspaceLER.keys():
                self._ler += self._estimated_subspaceLER[weight] * binomial_weight(
                    self._num_noise, weight, self._error_rate
                )
                # print("Weight: ",weight," LER: ",self._estimated_subspaceLER[weight]*binomial_weight(self._num_noise, weight,self._error_rate))
            else:
                fitted_subspace_LER = modified_sigmoid_function(
                    weight, self._a, self._b, self._c, self._t
                )
                self._ler += fitted_subspace_LER * binomial_weight(
                    self._num_noise, weight, self._error_rate
                )
                # print("Weight: ",weight," LER: ",fitted_subspace_LER*binomial_weight(self._num_noise, weight,self._error_rate))
            # self._ler+=scurve_function(weight,self._mu,self._sigma)*binomial_weight(self._num_noise,weight,self._error_rate)
        return self._ler

    def plot_scurve(self, filename=None, savefigure=False, title="S-curve"):
        """Plot the fitted S-curve overlaid on measured subspace P_L values.

        Generates a bar chart of measured P_L(w) with error bars, the
        fitted S-curve, and annotated regions (fault-tolerant, curve
        fitting, saturation, critical).

        Args:
            filename: Output path for the PDF. If ``None``, a default
                name is used.
            savefigure: If ``True``, save the figure to *filename*.
            title: Title string for the plot.
        """
        keys = list(self._estimated_subspaceLER.keys())
        values = [self._estimated_subspaceLER[k] for k in keys]
        sigma_list = [
            subspace_sigma_estimator(
                self._subspace_sample_used[k], self._subspace_LE_count[k]
            )
            for k in keys
        ]
        fig, ax = plt.subplots(figsize=(10, 6))

        # Calculate adaptive ranges
        x_min_data, x_max_data = min(keys), max(keys)
        y_max_data = max(values) if values else 0.5
        x_range = x_max_data - x_min_data if x_max_data != x_min_data else 1.0

        # Adaptive saturation region width
        saturation_width = max(3, x_range * 0.3)
        x_plot_max = self._saturatew + saturation_width

        # Adaptive bar width
        bar_width = (
            max(0.4, min(0.8, x_range / len(keys) * 0.6)) if len(keys) > 0 else 0.6
        )

        # bars ── discrete estimate
        ax.bar(
            keys,
            values,
            width=bar_width,
            color="tab:orange",
            alpha=0.8,
            label="Estimated subspace LER by sampling",
        )

        ax.errorbar(
            keys,
            values,
            yerr=sigma_list,
            fmt="none",
            ecolor="black",
            capsize=3,
            elinewidth=1,
            label="LER error bars",
        )

        # smooth S-curve
        x = np.linspace(self._t + 0.1, self._saturatew, 1000)
        y = modified_sigmoid_function(x, self._a, self._b, self._c, self._t)
        ax.plot(
            x,
            y,
            color="tab:blue",
            linewidth=2.0,
            label="Fitted S-curve",
            linestyle="--",
        )

        # Set axis limits
        ax.set_xlim(-0.5, x_plot_max)
        ax.set_ylim(0, y_max_data * 1.3)

        # Region text y position (use axes transform for consistent positioning)
        region_text_y = y_max_data * 1.15

        # Fault-tolerant area
        if self._t > 0:
            ax.axvspan(0, self._t, color="green", alpha=0.15)
            ax.text(
                self._t / 2,
                region_text_y,
                "Fault\ntolerant",
                ha="center",
                va="bottom",
                color="green",
                fontsize=8,
            )

        # Curve fitting region
        ax.axvspan(self._t, self._saturatew, color="yellow", alpha=0.10)
        ax.text(
            (self._t + self._saturatew) / 2,
            region_text_y,
            "Curve fitting",
            ha="center",
            va="bottom",
            fontsize=10,
        )

        # Saturation region
        ax.axvspan(self._saturatew, x_plot_max, color="red", alpha=0.15)
        ax.text(
            (self._saturatew + x_plot_max) / 2,
            region_text_y,
            "Saturation",
            ha="center",
            va="bottom",
            color="red",
            fontsize=10,
        )

        # Region: Critical area (gray)
        ax.axvspan(self._minw, self._maxw, color="gray", alpha=0.2)
        ax.axvline(
            self._minw, color="red", linestyle="--", linewidth=1.2, label=r"$w_{\min}$"
        )
        ax.axvline(
            self._maxw,
            color="green",
            linestyle="--",
            linewidth=1.2,
            label=r"$w_{\max}$",
        )
        ax.text(
            (self._minw + self._maxw) / 2,
            region_text_y,
            r"Critical Region",
            ha="center",
            va="bottom",
            fontsize=9,
        )

        # Labels and legend
        ax.set_xlabel("Weight")
        ax.set_ylabel("Logical Error Rate in subspace")
        ax.set_title(f"S-curve of {title} (PL={self._ler:.2e})")
        ax.legend(loc="upper left", fontsize=8)

        # Integer ticks on x-axis
        ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

        # Layout and save
        plt.tight_layout()
        if savefigure:
            fig.savefig(filename + ".pdf", format="pdf", bbox_inches="tight")
        plt.show()
        plt.close(fig)

    def set_t(self, t):
        """Set the fault-tolerant threshold *t*.

        This is ``(code_distance - 1) / 2`` and determines the pole
        location in the modified-linear model.

        Args:
            t: Fault-tolerant threshold value.
        """
        self._t = t

    def fast_calculate_LER_from_file(
        self, filepath, pvalue, codedistance, figname, titlename, repeat=1
    ):
        """End-to-end LER estimation with S-curve fitting from a STIM file.

        Runs the full pipeline: parse circuit, bracket the S-curve,
        sample subspaces, fit the linear and modified-linear models,
        compute LER, and print summary statistics over *repeat* trials.

        Args:
            filepath: Path to a STIM circuit file.
            pvalue: Physical error rate.
            codedistance: Code distance of the QEC code.
            figname: Base filename for diagnostic plots.
            titlename: Title string for plots.
            repeat: Number of independent trials.
        """
        self._error_rate = pvalue
        self._circuit_level_code_distance = codedistance
        ler_list = []
        sample_used_list = []
        r_squared_list = []
        Nerror_list = []
        time_list = []
        for i in range(repeat):
            self.clear_all()
            self.parse_from_file(filepath)
            start = time.time()
            self.determine_lower_w()
            self.determine_saturated_w()
            self.subspace_sampling()
            self.fit_linear_area()
            tmptime = time.time()

            self.fit_log_S_model(figname + "-R" + str(i) + "Final.pdf", tmptime - start)
            self.calc_logical_error_rate_after_curve_fitting()

            end = time.time()
            time_list.append(end - start)

            self.plot_scurve(figname + ".pdf", titlename)
            r_squared_list.append(self._R_square_score)
            self._sample_used = np.sum(list(self._subspace_sample_used.values()))
            # print("Final LER: ",self._ler)
            # print("Total samples used: ",self._sample_used)
            ler_list.append(self._ler)
            sample_used_list.append(self._sample_used)
            Nerror_list.append(sum(self._subspace_LE_count.values()))

        # Compute means
        self._ler = float(np.mean(ler_list))
        self._sample_used = int(np.mean(sample_used_list))

        # Compute standard deviations
        ler_std = float(np.std(ler_list))
        sample_used_std = float(np.std(sample_used_list))
        r2_mean = float(np.mean(r_squared_list))
        r2_std = float(np.std(r_squared_list))
        Nerror_mean = float(np.mean(Nerror_list))
        Nerror_std = float(np.std(Nerror_list))

        time_mean = float(np.mean(time_list))
        time_std = float(np.std(time_list))

        # Print with scientific ± formatting
        print("k: ", self._k_range)
        print("beta: ", self._beta)
        print("Subspaces: ", self._num_subspace)
        print("R2: ", format_with_uncertainty(r2_mean, r2_std))
        print(
            "Samples(ours): ",
            format_with_uncertainty(self._sample_used, sample_used_std),
        )
        print("Time(our): ", format_with_uncertainty(time_mean, time_std))
        print("PL(ours): ", format_with_uncertainty(self._ler, ler_std))
        print("Nerror(ours): ", format_with_uncertainty(Nerror_mean, Nerror_std))

    def sample_all_subspaces(self, Nclip, Budget, save_path=None):
        """Sample all subspaces from ``t+1`` to ``w_sat`` exhaustively.

        Adaptively samples each weight until either *Nclip* logical
        error events are collected or *Budget* samples have been used
        at that weight.

        Args:
            Nclip: Target number of logical error events per subspace.
            Budget: Maximum number of samples per subspace.
            save_path: If provided, pickle the results to this path.

        Returns:
            Tuple of ``(ler_count, sample_used, subspaceLER)`` dicts
            mapping weight to counts and estimated P_L values.
        """
        assert self._QEPG_graph is not None, (
            "QEPG graph must be initialized before sampling"
        )
        self.determine_saturated_w()
        wlist_to_sample = np.arange(self._t + 1, self._saturatew, step=1)
        sample_used: dict[int, int] = {}
        ler_count: dict[int, int] = {}
        subspaceLER: dict[int, float] = {}
        for w in wlist_to_sample:
            sample_used[int(w)] = 0
            ler_count[int(w)] = 0
            subspaceLER[int(w)] = 0.0
        """
            First round of sampling
            """
        slist = [8000] * len(wlist_to_sample)
        wlist: list[int] | np.ndarray = wlist_to_sample
        detector_result, obsresult = return_samples_many_weights_separate_obs_with_QEPG(
            self._QEPG_graph, list(wlist), slist
        )
        predictions_result = self._matcher.decode_batch(detector_result)
        begin_index = 0
        for w_idx, (w, quota) in enumerate(zip(wlist, slist)):
            observables = np.asarray(
                obsresult[begin_index : begin_index + quota]
            )  # (shots,)
            predictions = np.asarray(
                predictions_result[begin_index : begin_index + quota]
            ).ravel()
            # 3. count mismatches in vectorised form ---------------------------------
            num_errors = np.count_nonzero(observables != predictions)
            ler_count[w] += num_errors
            sample_used[w] += quota
            subspaceLER[w] = ler_count[w] / sample_used[w]
            begin_index += quota

        while True:
            slist = []
            wlist = []

            for weight in wlist_to_sample:
                if sample_used[weight] > Budget:
                    continue
                if ler_count[weight] < Nclip:
                    wlist.append(weight)
                    if ler_count[weight] == 0:
                        slist.append(min(10000, sample_used[weight] * 10))
                    else:
                        slist.append(
                            min(
                                10000,
                                int(
                                    sample_used[weight]
                                    * Nclip
                                    / (Nclip - ler_count[weight])
                                ),
                            )
                        )

            if len(wlist) == 0:
                break

            detector_result, obsresult = (
                return_samples_many_weights_separate_obs_with_QEPG(
                    self._QEPG_graph, wlist, slist
                )
            )
            predictions_result = self._matcher.decode_batch(detector_result)

            begin_index = 0
            for w_idx, (w, quota) in enumerate(zip(wlist, slist)):
                observables = np.asarray(
                    obsresult[begin_index : begin_index + quota]
                )  # (shots,)
                predictions = np.asarray(
                    predictions_result[begin_index : begin_index + quota]
                ).ravel()
                # 3. count mismatches in vectorised form ---------------------------------
                num_errors = np.count_nonzero(observables != predictions)
                ler_count[w] += num_errors
                sample_used[w] += quota
                subspaceLER[w] = ler_count[w] / sample_used[w]
                begin_index += quota

        result = {
            "ler_count": ler_count,
            "sample_used": sample_used,
            "subspaceLER": subspaceLER,
        }

        if save_path is not None:
            with open(save_path, "wb") as f:
                pickle.dump(result, f)

        return ler_count, sample_used, subspaceLER

    def load_all_sample_result(self, filepath):
        """Load previously saved subspace sampling results and plot.

        Reads a pickle file produced by :meth:`sample_all_subspaces`,
        plots the raw S-curve and a fitted log-logit model.

        Args:
            filepath: Path to the pickled results file.
        """
        with open(filepath, "rb") as f:
            result = pickle.load(f)

        ler_count = result["ler_count"]
        sample_used = result["sample_used"]
        subspaceLER = result["subspaceLER"]

        """Plot the S-curve and its discrete estimate."""
        keys = list(subspaceLER.keys())
        values = [subspaceLER[k] for k in keys]

        fig = plt.figure()
        # bars ── discrete estimate
        plt.bar(
            keys,
            values,
            color="tab:blue",  # pick any color you like
            alpha=0.8,
            label="Estimated subspace LER by sampling",
        )

        plt.axvline(
            x=self._error_rate * self._num_noise,
            color="red",
            linestyle="--",
            linewidth=1.2,
            label="Average Error number",
        )  # vertical line at x=0.5

        plt.xlabel("Weight")
        plt.ylabel("Logical Error Rate in subspace")
        plt.title("Scurve plot")
        plt.legend()  # <- shows the two labels

        # ── Force integer ticks on the X axis ──────────────────────────
        plt.gca().xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

        figname = "AllSubspace.pdf"
        plt.tight_layout()  # optional: nicely fit everything
        plt.savefig(figname, dpi=300)
        # plt.show()
        plt.close(fig)

        """Fit the curve and plot the fitted line."""
        # print("circuit d:",self._circuit_level_code_distance)
        x_list = [
            x
            for x in subspaceLER.keys()
            if (subspaceLER[x] < 0.5 and subspaceLER[x] > 0 and ler_count[x] > 100)
        ]

        sigma_list = [sigma_estimator(sample_used[x], ler_count[x]) for x in x_list]

        y_list = [np.log(0.5 / subspaceLER[x] - 1) for x in x_list]

        initial_guess = (0, 0, 0)

        lower = [-np.inf, -np.inf, -np.inf]
        # ── upper bounds for [param1, param2, param3, param4]
        upper = [np.inf, np.inf, np.inf]

        popt, pcov = curve_fit(
            modified_linear_function(self._t),
            x_list,
            y_list,
            p0=initial_guess,  # len(initial_guess) must be 4 and within the bounds above
            bounds=(lower, upper),  # <-- tuple with two arrays
            maxfev=50_000,  # or max_nfev in newer SciPy
        )

        # Extract the best-fit parameter (alpha)
        a, b, c = popt[0], popt[1], popt[2]

        # print("circuit d:",self._circuit_level_code_distance)
        y_list = [np.log(0.5 / subspaceLER[x] - 1) for x in x_list]
        y_predicted = [
            modified_linear_function_with_d(x, a, b, c, self._t) for x in x_list
        ]
        # y_predicted = [scurve_function_with_distance(x,self._mu,self._sigma) for x in self._estimated_wlist]
        R_square_score = r_squared(y_list, y_predicted)
        # print("R^2 score: ", self._R_square_score)

        # Plot the fitted line
        x_fit = np.linspace(self._t + 1, max(x_list), 1000)

        y_fit = modified_linear_function_with_d(x_fit, a, b, c, self._t)

        alpha = -1 / a

        sample_cost_list = [sample_used[x] for x in x_list]

        # print("Fitted parameters: a={}, b={}, c={}, d={}".format(self._a, self._b, self._c, self._d))
        plt.figure()
        plt.errorbar(
            x_list,
            y_list,
            yerr=sigma_list,
            fmt="o",
            color="orange",
            label="Data points with error bars",
            capsize=3,
            markersize=4,
            elinewidth=1,
        )
        plt.plot(x_fit, y_fit, label=f"Fitted line, R2={R_square_score}", color="blue")

        # Select 5 uniformly spaced indices from the available data points
        num_points_to_annotate = 5
        indices = np.linspace(0, len(x_list) - 1, num=num_points_to_annotate, dtype=int)

        # Annotate only the selected points
        for i in indices:
            x, y, s = x_list[i], y_list[i], sample_cost_list[i]
            if s == 0:
                continue
            s_str = "{0:.1e}".format(s)  # Scientific notation like 1.2e+03
            base, exp = s_str.split("e")
            exp = int(exp)
            label = r"${0}\times 10^{{{1}}}$".format(base, exp)
            plt.annotate(
                label,
                (x, y),
                textcoords="offset points",
                xytext=(0, 10),
                ha="center",
                fontsize=7,
            )

        text_lines = [
            r"$\alpha=%.4f$" % alpha,
            r"$\mu =%.4f$" % (alpha * b),
            r"$\beta=%.4f$" % c,
            r"$\#\mathrm{detector}=%d$" % self._num_detector,
            r"$\#\mathrm{noise}=%d$" % self._num_noise,
        ]

        textstr = "\n".join(text_lines)

        fig = plt.gcf()
        fig.subplots_adjust(right=0.75)  # Make room on the right
        fig.text(
            0.78,
            0.5,
            textstr,
            fontsize=7,
            va="center",
            ha="left",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.95),
        )
        # ──────────────────────────────────────────────────────────────────

        plt.xlabel("Weight")
        plt.ylabel(r"$\log\left(\frac{0.5}{\mathrm{LER}} - 1\right)$")
        plt.title("Linear Fit of S-curve")
        plt.legend(fontsize=9)
        plt.tight_layout()
        plt.savefig("yfunction.pdf")
        plt.close()

    def calculate_LER_from_file(
        self, filepath, pvalue, codedistance, figname, titlename, repeat=1
    ):
        self._error_rate = pvalue
        self._circuit_level_code_distance = codedistance
        ler_list = []
        sample_used_list = []
        r_squared_list = []
        Nerror_list = []
        time_list = []
        for i in range(repeat):
            self.clear_all()
            self.parse_from_file(filepath)
            start = time.time()
            # self.determine_range_to_sample()
            # self.subspace_sampling()
            self.determine_lower_w()

            # self.ground_truth_subspace_sampling()
            # self._has_logical_errorw=self._t+4
            self.determine_saturated_w()

            self.subspace_sampling_to_fit_curve(1000 * self._num_subspace)
            """
            Fit the curve first time just to get the estimated sweet spot
            """
            self.fit_linear_area()
            tmptime = time.time()
            self.fit_log_S_model(figname + "-R" + str(i) + "First.pdf", tmptime - start)
            """
            Second round of samples
            """
            self.subspace_sampling()

            self.fit_linear_area()
            tmptime = time.time()
            self.fit_log_S_model(figname + "-R" + str(i) + "Final.pdf", tmptime - start)

            self.calc_logical_error_rate_after_curve_fitting()

            end = time.time()
            time_list.append(end - start)

            self.plot_scurve(figname, titlename)
            r_squared_list.append(self._R_square_score)
            self._sample_used = np.sum(list(self._subspace_sample_used.values()))
            # print("Final LER: ",self._ler)
            # print("Total samples used: ",self._sample_used)
            ler_list.append(self._ler)
            sample_used_list.append(self._sample_used)
            Nerror_list.append(sum(self._subspace_LE_count.values()))

        # Compute means
        self._ler = float(np.mean(ler_list))
        self._sample_used = int(np.mean(sample_used_list))

        # Compute standard deviations
        ler_std = float(np.std(ler_list))
        sample_used_std = float(np.std(sample_used_list))
        r2_mean = float(np.mean(r_squared_list))
        r2_std = float(np.std(r_squared_list))
        Nerror_mean = float(np.mean(Nerror_list))
        Nerror_std = float(np.std(Nerror_list))

        time_mean = float(np.mean(time_list))
        time_std = float(np.std(time_list))

        # Print with scientific ± formatting
        print("k: ", self._k_range)
        print("beta: ", self._beta)
        print("Subspaces: ", self._num_subspace)
        print("R2: ", format_with_uncertainty(r2_mean, r2_std))
        print(
            "Samples(ours): ",
            format_with_uncertainty(self._sample_used, sample_used_std),
        )
        print("Time(our): ", format_with_uncertainty(time_mean, time_std))
        print("PL(ours): ", format_with_uncertainty(self._ler, ler_std))
        print("Nerror(ours): ", format_with_uncertainty(Nerror_mean, Nerror_std))

    def calculate_LER_from_StabCode(
        self,
        qeccirc: StabCode,
        noise_model: NoiseModel,
        figname: str | None = None,
        titlename: str | None = None,
        savefigure: bool = False,
        repeat: int = 1,
    ):
        qeccirc.construct_IR_standard_scheme()
        qeccirc.compile_stim_circuit_from_IR_standard()
        noisy_circuit = noise_model.reconstruct_clifford_circuit(qeccirc.circuit)

        self._error_rate = noise_model.error_rate

        self._circuit_level_code_distance = qeccirc.d
        self._t = max(0, int((self._circuit_level_code_distance - 1) / 2))
        ler_list = []
        sample_used_list = []
        r_squared_list = []
        Nerror_list = []
        time_list = []

        self._cliffordcircuit = noisy_circuit
        self._num_noise = self._cliffordcircuit.totalnoise
        self._num_detector = len(self._cliffordcircuit.parityMatchGroup)
        self._stim_str_after_rewrite = str(noisy_circuit.stimcircuit)
        # Configure a decoder using the circuit.
        self._detector_error_model = (
            self._cliffordcircuit.stimcircuit.detector_error_model(
                decompose_errors=True
            )
        )
        self._matcher = pymatching.Matching.from_detector_error_model(
            self._detector_error_model
        )
        self._QEPG_graph = compile_QEPG(str(noisy_circuit.stimcircuit))

        for i in range(repeat):
            self.clear_all()
            start = time.perf_counter()
            # self.determine_range_to_sample()
            # self.subspace_sampling()
            self.determine_lower_w()

            # self.ground_truth_subspace_sampling()
            # self._has_logical_errorw=self._t+4
            self.determine_saturated_w()

            self.subspace_sampling_to_fit_curve(1000 * self._num_subspace)
            """
            Fit the curve first time just to get the estimated sweet spot
            """
            self.fit_linear_area()
            tmptime = time.perf_counter()
            if savefigure:
                assert figname is not None
                self.fit_log_S_model(
                    figname + "-R" + str(i) + "First.pdf",
                    savefigure=savefigure,
                    time=tmptime - start,
                )
            else:
                self.fit_log_S_model(None, savefigure=savefigure, time=tmptime - start)
            """
            Second round of samples
            """
            self.subspace_sampling()

            self.fit_linear_area()
            tmptime = time.perf_counter()
            if savefigure:
                assert figname is not None
                self.fit_log_S_model(
                    figname + "-R" + str(i) + "Final.pdf",
                    savefigure=savefigure,
                    time=tmptime - start,
                )
            else:
                self.fit_log_S_model(None, savefigure=savefigure, time=tmptime - start)
            self.calc_logical_error_rate_after_curve_fitting()

            end = time.perf_counter()
            time_list.append(end - start)

            self.plot_scurve(figname, savefigure=savefigure, title=titlename)
            r_squared_list.append(self._R_square_score)
            self._sample_used = np.sum(list(self._subspace_sample_used.values()))
            # print("Final LER: ",self._ler)
            # print("Total samples used: ",self._sample_used)
            ler_list.append(self._ler)
            sample_used_list.append(self._sample_used)
            Nerror_list.append(sum(self._subspace_LE_count.values()))

        # Compute means
        self._ler = float(np.mean(ler_list))
        self._sample_used = int(np.mean(sample_used_list))

        # Compute standard deviations
        ler_std = float(np.std(ler_list))
        sample_used_std = float(np.std(sample_used_list))
        r2_mean = float(np.mean(r_squared_list))
        r2_std = float(np.std(r_squared_list))
        Nerror_mean = float(np.mean(Nerror_list))
        Nerror_std = float(np.std(Nerror_list))

        time_mean = float(np.mean(time_list))
        time_std = float(np.std(time_list))

        # Print with scientific ± formatting
        print("k: ", self._k_range)
        print("beta: ", self._beta)
        print("Subspaces: ", self._num_subspace)
        print("R2: ", format_with_uncertainty(r2_mean, r2_std))
        print(
            "Samples(ours): ",
            format_with_uncertainty(self._sample_used, sample_used_std),
        )
        print("Time(our): ", format_with_uncertainty(time_mean, time_std))
        print("PL(ours): ", format_with_uncertainty(self._ler, ler_std))
        print("Nerror(ours): ", format_with_uncertainty(Nerror_mean, Nerror_std))


if __name__ == "__main__":
    p = 0.001
    sample_budget = 100_000_0000

    for d in range(7, 9, 2):
        t = (d - 1) // 2

        stim_path = (
            f"C:/Users/username/Documents/Sampling/stimprograms/surface/surface{d}"
        )
        figname = f"Surface{d}"
        titlename = f"Surface{d}"
        output_filename = f"Surface{d}.txt"

        tmp = StratifiedScurveLERcalc(
            p, sampleBudget=sample_budget, k_range=5, num_subspace=6, beta=4
        )
        tmp.set_t(t)
        tmp.set_sample_bound(
            MIN_NUM_LE_EVENT=100,
            SAMPLE_GAP=100,
            MAX_SAMPLE_GAP=5000,
            MAX_SUBSPACE_SAMPLE=50000,
        )
        with open(output_filename, "w") as f:
            with redirect_stdout(f):
                tmp.calculate_LER_from_file(stim_path, p, 0, figname, titlename, 1)
