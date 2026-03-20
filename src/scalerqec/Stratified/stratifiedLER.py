"""Legacy stratified sampling LER estimator (without S-curve fitting).

This module provides :class:`StratifiedLERcalc`, which estimates the
logical error rate by directly sampling each weight subspace and
aggregating via binomial weights.  No curve fitting is performed, so
every contributing weight must be sampled with enough shots to
accumulate a minimum number of logical error events.

For codes with many noise locations the newer :class:`Scaler` class
(which fits an S-curve model to interpolate unsampled weights) is
strongly preferred.
"""

from __future__ import annotations
from typing import Optional
import numpy as np
from ..qepg import (
    return_samples_many_weights_separate_obs,
    compile_QEPG,
    return_samples_many_weights_separate_obs_with_QEPG,
    QEPGGraph,
    return_samples_with_noise_vector,
)
from ..Clifford.clifford import *
import pymatching
import time
from ..QEC.noisemodel import NoiseModel
from ..QEC.qeccircuit import StabCode
from ..util import binomial_weight, subspace_size, format_with_uncertainty


MIN_NUM_LE_EVENT = 1000  # Increased from 50 for better accuracy
SAMPLE_GAP = 500


class StratifiedLERcalc:
    """Stratified sampling LER estimator (without S-curve fitting).

    This is the legacy stratified sampler that estimates the logical
    error rate by directly sampling each weight subspace until a
    sufficient number of logical error events have been observed, then
    aggregating via the binomial-weight formula::

        LER = sum_w  P_L(w) * Binom(N, w) * p^w * (1-p)^(N-w)

    Unlike :class:`~scalerqec.Stratified.Scaler.Scaler`, this class
    does **not** fit an S-curve model and instead relies on brute-force
    sampling at every weight in a heuristically chosen range. It is
    controlled by a *sample budget* (total number of shots) rather than
    a wall-clock time budget.

    Args:
        error_rate: Physical error rate per noise location.
        sampleBudget: Maximum total number of samples (shots) across
            all weight subspaces.
        num_subspace: Number of weight subspaces to sample.
    """

    def __init__(
        self, error_rate: float = 0, sampleBudget: int = 10000, num_subspace: int = 30,
        decoder=None,
    ):
        self._num_detector: int = 0
        self._num_noise: int = 0
        self._error_rate: float = error_rate
        self._decoder = decoder
        self._cliffordcircuit: CliffordCircuit = CliffordCircuit(4)

        self._ler: float = 0
        """
        Use a dictionary to store the estimated subspace logical error rate
        """
        self._estimated_subspaceLER: dict[int, float] = {}
        self._subspace_LE_count: dict[int, int] = {}
        self._subspace_sample_used: dict[int, int] = {}

        self._sampleBudget: int = sampleBudget
        self._num_subspace: int = num_subspace
        self._minW: int = 0
        self._maxW: int = 0

        self._stim_str_after_rewrite: str = ""

        self._sample_used: int = 0
        self._uncertainty: float = 0

        self._circuit_level_code_distance: int = 1

        self._QEPG_graph: Optional[QEPGGraph] = None

    def parse_from_file(self, filepath: str):
        """Load a STIM circuit from *filepath* and prepare for sampling.

        Compiles the circuit into a Clifford representation, builds the
        QEPG graph for efficient weight-stratified sampling, and
        initialises a PyMatching decoder from the detector error model.

        Args:
            filepath: Path to a STIM circuit file.
        """
        stim_str = ""
        with open(filepath, "r", encoding="utf-8") as f:
            stim_str = f.read()

        self._cliffordcircuit.compile_from_stim_circuit_str(stim_str)
        self._stim_str_after_rewrite = stim_str

        # Compile QEPG graph and get the correct noise count from it
        # IMPORTANT: Use QEPG's noise count, not clifford's totalnoise
        # They differ because QEPG merges noise operations
        self._QEPG_graph = compile_QEPG(stim_str)
        self._num_noise = self._QEPG_graph.get_total_noise()
        self._num_detector = self._QEPG_graph.get_total_detector()

        # Inject noise for decoder (DEM needs noisy circuit)
        from ..QEC.noisemodel import SIDNoiseModel

        noisy_stim = SIDNoiseModel(self._error_rate).inject_noise(
            self._cliffordcircuit.stimcircuit
        )

        # Configure a decoder using the noisy circuit.
        self._detector_error_model = noisy_stim.detector_error_model(
            decompose_errors=True
        )
        if self._decoder is not None:
            self._matcher = self._decoder
        else:
            self._matcher = pymatching.Matching.from_detector_error_model(
                self._detector_error_model
            )

    def sample_all_subspace(self, shots_each_subspace: int = 1000000):
        """Sample every weight subspace with a fixed number of shots.

        This is a brute-force ground-truth method intended for
        validating the accuracy of the stratified algorithm. It samples
        **all** weights from 0 to *N* (the total number of noise
        locations).

        Args:
            shots_each_subspace: Number of samples drawn per weight.
        """
        wlist = list(range(0, self._num_noise + 1))
        slist = [shots_each_subspace] * len(wlist)
        detector_result, obsresult = return_samples_many_weights_separate_obs(
            self._stim_str_after_rewrite, wlist, slist
        )
        predictions_result = self._matcher.decode_batch(detector_result)

        for w in wlist:
            self._subspace_LE_count[w] = 0
            self._estimated_subspaceLER[w] = 0
            self._subspace_sample_used[w] = shots_each_subspace

        begin_index = 0
        for w_idx, (w, quota) in enumerate(zip(wlist, slist)):
            observables = np.asarray(
                obsresult[begin_index : begin_index + quota]
            ).ravel()  # (shots,)
            # 2. batch-decode (decode_batch should accept ndarray) -------------------
            # shape (shots,) or (shots,1)
            predictions = np.asarray(
                predictions_result[begin_index : begin_index + quota]
            ).ravel()

            # 3. count mismatches in vectorised form ---------------------------------
            num_errors = np.count_nonzero(observables != predictions)

            self._subspace_LE_count[w] += num_errors
            self._estimated_subspaceLER[w] = (
                self._subspace_LE_count[w] / self._subspace_sample_used[w]
            )

            # print(f"Subspace logical error {self._estimated_subspaceLER[w]}")
            # print(f"Logical error rate when w={w}: {self._estimated_subspaceLER[w]*binomial_weight(self._num_noise, w,self._error_rate):.6g}")
            begin_index += quota

    def sample_all_subspace_sequential(self, shots_each_subspace: int = 1000000):
        """Sample every weight subspace sequentially for verification.

        Like :meth:`sample_all_subspace`, but uses the sequential C++
        ``return_samples_with_noise_vector`` function instead of the
        batched multi-weight sampler. Primarily used to verify that the
        two sampling backends produce consistent results.

        Args:
            shots_each_subspace: Number of samples drawn per weight.
        """
        wlist = list(range(0, self._num_noise + 1))

        # Initialize counters
        for w in wlist:
            self._subspace_LE_count[w] = 0
            self._estimated_subspaceLER[w] = 0
            self._subspace_sample_used[w] = shots_each_subspace

        # Sample each weight separately using the sequential function
        for w in wlist:
            if w == 0:
                # Weight 0: no errors, no logical errors
                self._subspace_LE_count[w] = 0
                self._estimated_subspaceLER[w] = 0.0
                continue

            # Use return_samples_with_noise_vector for this weight
            # Returns (noise_vectors, results) where results is list of [det0, det1, ..., detN, observable]
            _noise_vectors, results = return_samples_with_noise_vector(
                self._stim_str_after_rewrite, w, shots_each_subspace
            )

            # Convert to numpy array for easier manipulation
            results_array = np.array(
                results, dtype=bool
            )  # shape: (shots, num_detectors+1)

            # Split into detectors and observables
            detectors = results_array[:, :-1]  # All columns except last
            observables = results_array[:, -1]  # Last column

            # Decode detectors
            predictions = self._matcher.decode_batch(detectors)
            predictions = np.asarray(predictions).ravel()
            observables = observables.ravel()

            # Count mismatches
            num_errors = np.count_nonzero(observables != predictions)

            self._subspace_LE_count[w] = num_errors
            self._estimated_subspaceLER[w] = (
                self._subspace_LE_count[w] / self._subspace_sample_used[w]
            )

    def determine_range_to_sample(self):
        """Determine the weight range [minW, maxW] to sample.

        Sets ``_minW`` and ``_maxW`` to cover the weights whose
        binomial contribution ``Binom(N, w) * p^w * (1-p)^(N-w)`` is
        non-negligible.  The range is centred on the expected number
        of faults ``N * p`` and extends +/- 5 standard deviations,
        with a minimum span of 10 weights to ensure adequate coverage.
        """
        sigma = int(
            (self._error_rate * (1 - self._error_rate) * self._num_noise) ** 0.5
        )
        if sigma == 0:
            sigma = 1
        ep = int(self._error_rate * self._num_noise)
        self._minW = max(1, ep - 5 * sigma)
        # Cap at num_noise to avoid exceeding available noise sources
        # The C++ sampler now handles weight == num_noise correctly via removal strategy
        self._maxW = max(2, min(self._num_noise, ep + 5 * sigma))

        # Ensure we sample at least 10 weights for better coverage
        if self._maxW - self._minW + 1 < 10:
            # Expand range symmetrically around the mean
            needed = 10 - (self._maxW - self._minW + 1)
            expand_left = needed // 2
            expand_right = needed - expand_left
            self._minW = max(1, self._minW - expand_left)
            self._maxW = min(self._num_noise, self._maxW + expand_right)

    def subspace_sampling(self):
        """Adaptively sample weight subspaces to estimate P_L(w).

        Iteratively allocates shots to weight subspaces, prioritising
        those that have not yet accumulated enough logical error events
        (``MIN_NUM_LE_EVENT``). Sampling stops when either:

        1. Every subspace has at least ``MIN_NUM_LE_EVENT`` logical
           error events, or
        2. The total sample budget ``_sampleBudget`` is exhausted.

        The method also detects the circuit-level code distance by
        identifying weights where exhaustive sampling finds no errors.
        """
        self.determine_range_to_sample()

        # print("Sampling weights from {} to {}".format(self._minW,self._maxW))
        """
        wlist store the subset of weights we need to sample and get
        correct logical error rate.
        
        In each subspace, we stop sampling until 100 logical error events are detected, or we hit the total budget.
        """
        wlist_need_to_sample = list(range(self._minW, self._maxW + 1))
        self._sample_used = 0
        for weight in wlist_need_to_sample:
            self._subspace_LE_count[weight] = 0
            self._subspace_sample_used[weight] = 0

        # print("Weights need to sample: ")
        # print(wlist_need_to_sample)

        min_non_zero_weight = 1e9
        while True:
            slist = []
            wlist = []
            """
            Case 1 to end the while loop: We have consumed all of our sample budgets
            """
            if self._sample_used > self._sampleBudget:
                break

            for weight in wlist_need_to_sample:
                if weight <= (self._circuit_level_code_distance - 1) / 2:
                    continue
                """
                If the subspace has been sampled enough, but logical error rate is still zero,
                also, the number of samples used in the subspace is comparable with the size of the subspace,
                we declare that the code distance is larger than the current weight.
                """
                if weight + 1 < wlist_need_to_sample[-1]:
                    # if(weight>min_non_zero_weight):
                    #     continue
                    if self._subspace_LE_count[weight] == 0:
                        if (
                            subspace_size(self._num_noise, weight)
                            < 2 * self._subspace_sample_used[weight]
                        ):
                            self._circuit_level_code_distance = weight
                            continue
                        if (
                            self._subspace_LE_count[weight + 1] >= MIN_NUM_LE_EVENT
                            and 2 * self._subspace_sample_used[weight + 1]
                            < self._subspace_sample_used[weight]
                        ):
                            self._circuit_level_code_distance = weight
                            continue
                        else:
                            slist.append(
                                max(
                                    2 * self._subspace_sample_used[weight + 1],
                                    SAMPLE_GAP,
                                )
                            )
                            wlist.append(weight)
                            self._subspace_sample_used[weight] += max(
                                2 * self._subspace_sample_used[weight + 1], SAMPLE_GAP
                            )
                            self._sample_used += max(
                                2 * self._subspace_sample_used[weight + 1], SAMPLE_GAP
                            )
                            continue

                if self._subspace_LE_count[weight] > 0:
                    min_non_zero_weight = min(weight, min_non_zero_weight)

                if self._subspace_LE_count[weight] < MIN_NUM_LE_EVENT:
                    if self._subspace_LE_count[weight] >= 1:
                        sample_num_required = (
                            int(MIN_NUM_LE_EVENT / self._subspace_LE_count[weight])
                            * self._subspace_sample_used[weight]
                        )
                        slist.append(sample_num_required)
                        self._subspace_sample_used[weight] += sample_num_required
                        self._sample_used += sample_num_required
                    else:
                        slist.append(SAMPLE_GAP)
                        self._subspace_sample_used[weight] += SAMPLE_GAP
                        self._sample_used += SAMPLE_GAP
                    wlist.append(weight)
            """
            Case 2 to end the while loop: We have get 100 logical error events for all these subspaces
            """
            if len(wlist) == 0:
                break

            # print("wlist: ",wlist)
            # print("slist: ",slist)
            # detector_result,obsresult=return_samples_many_weights_separate_obs(self._stim_str_after_rewrite,wlist,slist)
            assert self._QEPG_graph is not None, (
                "QEPG graph must be initialized before sampling"
            )

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
                ).ravel()  # (shots,)
                predictions = np.asarray(
                    predictions_result[begin_index : begin_index + quota]
                ).ravel()

                # 3. count mismatches in vectorised form ---------------------------------
                num_errors = np.count_nonzero(observables != predictions)

                self._subspace_LE_count[w] += num_errors
                self._estimated_subspaceLER[w] = (
                    self._subspace_LE_count[w] / self._subspace_sample_used[w]
                )

                # print(f"Logical error rate when w={w}: {self._estimated_subspaceLER[w]*binomial_weight(self._num_noise, w,self._error_rate):.6g}")

                begin_index += quota
            # print(self._subspace_LE_count)
            # print(self._subspace_sample_used)
        # print("Samples used:{}".format(self._sample_used))
        # print("circuit level code distance:{}".format(self._circuit_level_code_distance))

    def calculate_LER(self, debug=False):
        """Compute the total logical error rate from subspace estimates.

        Aggregates the per-weight conditional logical error probabilities
        using the formula::

            LER = sum_w  P_L(w) * Binom(N, w) * p^w * (1-p)^(N-w)

        Args:
            debug: If ``True``, print a detailed table of each weight's
                contribution and the top-5 contributing weights.

        Returns:
            The estimated total logical error rate.
        """
        self._ler: float = 0

        if debug:
            print(f"\n{'=' * 90}")
            print(f"DEBUG: calculate_LER() - Combining subspace results")
            print(f"{'=' * 90}")
            print(f"  Total noise sources: {self._num_noise}")
            print(f"  Error rate: {self._error_rate}")
            print(
                f"\n  {'Weight':<8} {'Samples':<12} {'LE Count':<12} {'P(LE|w)':<15} {'P(w)':<15} {'Contribution':<15}"
            )
            print(f"  {'-' * 8} {'-' * 12} {'-' * 12} {'-' * 15} {'-' * 15} {'-' * 15}")

        contributions = []
        for weight in range(1, self._num_noise + 1):
            if weight in self._estimated_subspaceLER.keys():
                p_le_given_w = self._estimated_subspaceLER[weight]
                p_w = binomial_weight(self._num_noise, weight, self._error_rate)
                contribution = p_le_given_w * p_w
                self._ler += contribution

                if debug:
                    samples = self._subspace_sample_used.get(weight, 0)
                    le_count = self._subspace_LE_count.get(weight, 0)
                    print(
                        f"  {weight:<8} {samples:<12,} {le_count:<12,} {p_le_given_w:<15.6e} {p_w:<15.6e} {contribution:<15.6e}"
                    )
                    contributions.append((weight, contribution))

        if debug:
            print(f"  {'-' * 90}")
            print(f"  Total LER: {self._ler:.6e}")

            # Show top 5 contributors
            contributions.sort(key=lambda x: x[1], reverse=True)
            print(f"\n  Top 5 contributing weights:")
            for i, (w, contrib) in enumerate(contributions[:5]):
                pct = 100 * contrib / self._ler if self._ler > 0 else 0
                print(f"    {i + 1}. Weight {w}: {contrib:.6e} ({pct:.1f}% of total)")
            print(f"{'=' * 90}\n")

        return self._ler

    def get_LER_subspace_no_weight(self, weight: int):
        """Return the raw conditional logical error probability P_L(w).

        Args:
            weight: The fault weight.

        Returns:
            P_L(weight) without the binomial weight factor.
        """
        return self._estimated_subspaceLER[weight]

    def get_LER_subspace(self, weight: int):
        """Return the LER contribution of a single weight subspace.

        Computes ``P_L(w) * Binom(N, w) * p^w * (1-p)^(N-w)``.

        Args:
            weight: The fault weight.

        Returns:
            The weighted contribution of this subspace to the total LER.
        """
        return self._estimated_subspaceLER[weight] * binomial_weight(
            self._num_noise, weight, self._error_rate
        )

    def calculate_LER_from_file(self, filepath: str, pvalue: float, repeat: int):
        """End-to-end LER estimation from a STIM circuit file.

        Parses the circuit, runs stratified sampling, computes the LER,
        and prints summary statistics (mean +/- std) over *repeat*
        independent trials.

        Args:
            filepath: Path to a STIM circuit file.
            pvalue: Physical error rate.
            repeat: Number of independent trials to run.
        """
        self.parse_from_file(filepath)
        self._error_rate = pvalue
        ler_list: list[float] = []
        sample_used_list: list[int] = []
        time_list: list[float] = []

        for i in range(repeat):
            starttime = time.perf_counter()

            self.subspace_sampling()
            self.calculate_LER()
            ler_list.append(self._ler)
            sample_used_list.append(self._sample_used)
            endtime = time.perf_counter()
            time_list.append(endtime - starttime)

        average_LER = sum(ler_list) / len(ler_list)
        average_sample_used = sum(sample_used_list) / len(sample_used_list)
        time_mean = sum(time_list) / len(time_list)
        ler_std: float = float(np.std(ler_list))
        sample_used_std: float = float(np.std(sample_used_list))
        time_std: float = float(np.std(time_list))
        print(
            "Samples(ours): ",
            format_with_uncertainty(average_sample_used, sample_used_std),
        )
        print("Time(our): ", format_with_uncertainty(time_mean, time_std))
        print("PL(ours): ", format_with_uncertainty(average_LER, ler_std))

    def clear_all(self):
        """Reset internal state for a fresh estimation run."""
        pass

    def calculate_LER_from_StabCode(
        self, qeccirc: StabCode, noise_model: NoiseModel, repeat: int = 1
    ):
        """Estimate LER directly from a :class:`StabCode` and noise model.

        Constructs the noisy STIM circuit internally, then runs the
        stratified sampling pipeline.

        Args:
            qeccirc: A stabiliser code with stabilisers and logical
                operators already configured.
            noise_model: The noise model to apply to the circuit.
            repeat: Number of independent trials.
        """
        qeccirc.construct_IR_standard_scheme()
        qeccirc.compile_stim_circuit_from_IR_standard()
        noisy_circuit = noise_model.reconstruct_clifford_circuit(qeccirc.circuit)
        self._error_rate = noise_model.error_rate
        self._circuit_level_code_distance = qeccirc.d
        self._cliffordcircuit = noisy_circuit
        self._num_noise = self._cliffordcircuit.totalnoise
        self._num_detector = len(self._cliffordcircuit.parityMatchGroup)
        self._stim_str_after_rewrite = str(self._cliffordcircuit.stimcircuit)
        # Configure a decoder using the circuit.
        self._detector_error_model = (
            self._cliffordcircuit.stimcircuit.detector_error_model(
                decompose_errors=True
            )
        )
        if self._decoder is not None:
            self._matcher = self._decoder
        else:
            self._matcher = pymatching.Matching.from_detector_error_model(
                self._detector_error_model
            )
        self._QEPG_graph = compile_QEPG(self._stim_str_after_rewrite)

        ler_list: list[float] = []
        sample_used_list: list[int] = []
        time_list: list[float] = []

        for i in range(repeat):
            starttime = time.perf_counter()

            self.subspace_sampling()
            self.calculate_LER()
            ler_list.append(self._ler)
            sample_used_list.append(self._sample_used)
            endtime = time.perf_counter()
            time_list.append(endtime - starttime)

        average_LER = sum(ler_list) / len(ler_list)
        average_sample_used = sum(sample_used_list) / len(sample_used_list)
        time_mean = sum(time_list) / len(time_list)
        ler_std: float = float(np.std(ler_list))
        sample_used_std: float = float(np.std(sample_used_list))
        time_std: float = float(np.std(time_list))
        print(
            "Samples(ours): ",
            format_with_uncertainty(average_sample_used, sample_used_std),
        )
        print("Time(our): ", format_with_uncertainty(time_mean, time_std))
        print("PL(ours): ", format_with_uncertainty(average_LER, ler_std))


if __name__ == "__main__":
    # tmp=StratifiedLERcalc(0.001,sampleBudget=15000000,num_subspace=5)
    # filepath="C:/Users/username/Documents/Sampling/stimprograms/small/simple"
    # tmp.parse_from_file(filepath)
    # tmp.sample_all_subspace(11*1000000)

    # LER=tmp.calculate_LER()

    # print(LER)

    # for weight in range(1,12):
    #     #print("LER in the subspace {} is {}".format(weight,tmp.get_LER_subspace_no_weight(weight)))
    #     print("LER in the subspace {} is {}".format(weight,tmp.get_LER_subspace(weight)))

    qeccirc = StabCode(n=3, k=1, d=3)
    # Stabilizer generators
    qeccirc.add_stab("ZZI")
    qeccirc.add_stab("IZZ")
    qeccirc.set_logical_Z(0, "ZZZ")
    noise_model = NoiseModel(0.001)  # Set the noise model
    # Set stabilizer parity measurement scheme, round of repetition
    qeccirc.scheme = "Standard"
    qeccirc.rounds = 2

    stratifiedcalculator = StratifiedLERcalc()
    stratifiedcalculator.calculate_LER_from_StabCode(qeccirc, noise_model, repeat=1)
