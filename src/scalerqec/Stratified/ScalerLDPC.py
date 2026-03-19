from __future__ import annotations

# ScaLERLDPC: ScaLER for LDPC codes using BP+OSD decoder
# This class extends Scaler to support LDPC codes by replacing
# PyMatching with the BPOSD decoder from stimbposd.

import numpy as np
import time
from typing import List, Optional, Dict
from scalerqec.Stratified.Scaler import Scaler
from scalerqec.Stratified.models import ModelType
from scalerqec.qepg import (
    compile_QEPG,
    return_samples_many_weights_separate_obs_with_QEPG,
    return_samples_with_fixed_QEPG,
)

try:
    from stimbposd import BPOSD
except ImportError:
    raise ImportError(
        "stimbposd is required for ScalerLDPC. Install it with: pip install stimbposd"
    )


class ScalerLDPC(Scaler):
    """
    ScaLER for LDPC codes using BP+OSD decoder.

    This class extends Scaler to support LDPC codes by replacing
    PyMatching with the BPOSD decoder from stimbposd.

    The BPOSD decoder uses Belief Propagation followed by Ordered
    Statistics Decoding, which works well for LDPC codes where
    minimum-weight perfect matching is not applicable.

    User-facing hyperparameters (in addition to Scaler's):
      - max_bp_iters: Maximum BP iterations (default: 20)
      - bp_method: BP method - 'product_sum' or 'min_sum' (default: 'product_sum')
      - osd_method: OSD method - 'osd_e', 'osd_cs', 'osd_0' (default: 'osd_cs')
      - osd_order: OSD order (default: 0, use higher for better accuracy but slower)
    """

    # Maximum samples per batch for BPOSD (slower decoder)
    _MAX_BATCH_SIZE_LDPC: int = 500

    def __init__(
        self,
        error_rate: float = 0.0,
        time_budget: int = 30,
        model_type: ModelType = ModelType.OUR_MODEL,
        gamma: float = 1,
        num_subspaces_phase2: int = 12,
        binary_search_shots: int = 100,
        # BPOSD-specific parameters
        max_bp_iters: int = 20,
        bp_method: str = "product_sum",
        osd_method: str = "osd0",  # Use 'osd0' for order-0 OSD (fast)
        osd_order: int = 0,
    ):
        """
        Initialize the ScalerLDPC.

        Args:
            error_rate: Physical error rate (probability of error per gate)
            time_budget: Time budget in seconds
            model_type: Which S-curve model to use
            gamma: Sweet spot tuning parameter for d²y/dw² = γ * dy/dw
            num_subspaces_phase2: Number of subspaces to sample between w_sweet and w_err (default: 12)
            binary_search_shots: Number of shots per binary search iteration (default: 100)
            max_bp_iters: Maximum number of BP iterations
            bp_method: BP method ('product_sum' or 'min_sum')
            osd_method: OSD method ('osd0', 'osd_e', 'osd_cs')
            osd_order: OSD order (only used with osd_e/osd_cs, higher for better accuracy)
        """
        super().__init__(
            error_rate, time_budget, model_type, gamma, num_subspaces_phase2
        )

        # Binary search configuration
        self._binary_search_shots: int = binary_search_shots

        # BPOSD-specific parameters
        self._max_bp_iters: int = max_bp_iters
        self._bp_method: str = bp_method
        self._osd_method: str = osd_method
        self._osd_order: int = osd_order

        # Replace _matcher with _decoder
        self._decoder: Optional[BPOSD] = None

    # ------------------------------------------------------------------
    #  Circuit / QEPG setup - Override to use BPOSD
    # ------------------------------------------------------------------

    def parse_from_file(self, filepath: str) -> None:
        """
        Read the circuit, parse from the file, compile stim circuit and QEPG graph.
        Uses BPOSD decoder instead of PyMatching.
        """
        print(f"\n[ScalerLDPC] Parsing circuit from: {filepath}")
        with open(filepath, "r", encoding="utf-8") as f:
            stim_str = f.read()

        print("[ScalerLDPC] Compiling circuit...")
        self._cliffordcircuit.compile_from_stim_circuit_str(stim_str)
        self._num_noise = self._cliffordcircuit.totalnoise
        self._num_detector = len(self._cliffordcircuit.parityMatchGroup)
        self._stim_str_after_rewrite = stim_str

        print(f"[ScalerLDPC] Circuit compiled:")
        print(f"  - Number of noise locations: {self._num_noise}")
        print(f"  - Number of detectors: {self._num_detector}")

        # Inject noise for decoder (DEM needs noisy circuit)
        from scalerqec.QEC.noisemodel import SIDNoiseModel

        noisy_stim = SIDNoiseModel(self._error_rate).inject_noise(
            self._cliffordcircuit.stimcircuit
        )

        # Configure BPOSD decoder using the noisy detector error model
        print("[ScalerLDPC] Building detector error model...")
        self._detector_error_model = noisy_stim.detector_error_model(
            decompose_errors=False
        )

        # Initialize BPOSD decoder
        print(f"[ScalerLDPC] Initializing BPOSD decoder...")
        print(f"  - max_bp_iters: {self._max_bp_iters}")
        print(f"  - bp_method: {self._bp_method}")
        print(f"  - osd_method: {self._osd_method}")
        print(f"  - osd_order: {self._osd_order}")
        self._decoder = BPOSD(
            self._detector_error_model,
            max_bp_iters=self._max_bp_iters,
            bp_method=self._bp_method,
            osd_method=self._osd_method,
            osd_order=self._osd_order,
        )
        print("[ScalerLDPC] BPOSD decoder initialized successfully!")

        # Compile QEPG graph once
        print("[ScalerLDPC] Compiling QEPG graph...")
        self._QEPG_graph = compile_QEPG(stim_str)
        print("[ScalerLDPC] QEPG graph compiled!")

        # Set _matcher to None to ensure we don't accidentally use it
        self._matcher = None

    # ------------------------------------------------------------------
    #  Binary search helpers - Override with verbose output
    # ------------------------------------------------------------------

    def binary_search_upper(self, low: int, high: int, shots: int) -> int:
        """Find the smallest w in [low, high] such that PL(w) > _max_PL. With verbose output."""
        print(f"\n[ScalerLDPC] Binary search UPPER (saturation point)")
        print(f"  Range: [{low}, {high}], shots={shots}, threshold={self._max_PL}")
        left = low
        right = high
        epsilon = self._max_PL
        iteration = 0
        while left < right:
            mid = (left + right) // 2
            iteration += 1
            print(f"  Iter {iteration}: testing w={mid}...", end=" ", flush=True)
            er = self.calc_logical_error_rate_with_fixed_w(shots, mid)
            print(f"PL={er:.4f}", end=" ")
            if er > epsilon:
                print(f"-> search left [{left}, {mid}]")
                right = mid
            else:
                print(f"-> search right [{mid + 1}, {right}]")
                left = mid + 1
        print(f"  Result: w_sat = {left}")
        return left

    def binary_search_lower(
        self, low: int, high: int, shots: int | None = None, epsilon: float = 0.002
    ) -> int:
        """Find the smallest w in [low, high] such that PL(w) > epsilon. With verbose output."""
        if shots is None:
            shots = self._binary_search_shots
        print(f"\n[ScalerLDPC] Binary search LOWER (first logical error)")
        print(f"  Range: [{low}, {high}], shots={shots}, threshold={epsilon}")
        left = low
        right = high
        iteration = 0
        while left < right:
            mid = (left + right) // 2
            iteration += 1
            print(f"  Iter {iteration}: testing w={mid}...", end=" ", flush=True)
            er = self.calc_logical_error_rate_with_fixed_w(shots, mid)
            print(f"PL={er:.4f}", end=" ")
            if er > epsilon:
                print(f"-> search left [{left}, {mid}]")
                right = mid
            else:
                print(f"-> search right [{mid + 1}, {right}]")
                left = mid + 1
        print(f"  Result: w_has_error = {left}")
        return left

    def determine_lower_w(self) -> None:
        """Determine the first weight where PL is noticeably non-zero. With verbose output."""
        print("\n[ScalerLDPC] Determining lower weight (w_has_error)...")
        if self._num_noise <= 8:
            self._has_logical_errorw = self._t + 1
            print(
                f"  Small circuit, using w_has_error = t+1 = {self._has_logical_errorw}"
            )
        else:
            self._has_logical_errorw = self.binary_search_lower(
                self._t + 1, self._num_noise
            )

    def determine_saturated_w(self) -> None:
        """Determine the weight where PL is essentially saturated. With verbose output."""
        print("\n[ScalerLDPC] Determining saturated weight (w_sat)...")
        if self._num_noise <= 8:
            self._saturatew = self._num_noise
            print(f"  Small circuit, using w_sat = num_noise = {self._saturatew}")
        else:
            self._saturatew = self.binary_search_upper(
                self._has_logical_errorw,
                self._num_noise,
                shots=self._binary_search_shots,
            )
            if self._saturatew < self._has_logical_errorw + 8:
                self._saturatew = min(self._num_noise, self._has_logical_errorw + 8)
                print(f"  Adjusted w_sat to ensure minimum range: {self._saturatew}")

    # ------------------------------------------------------------------
    #  Logical error rate calculation - Override to use BPOSD
    # ------------------------------------------------------------------

    def calc_logical_error_rate_with_fixed_w(self, shots: int, w: int) -> float:
        """
        Calculate the logical error rate with fixed Pauli weight w.
        Uses BPOSD decoder.
        """
        assert self._QEPG_graph is not None, (
            "QEPG graph must be initialized before sampling"
        )
        assert self._decoder is not None, (
            "BPOSD decoder must be initialized before decoding"
        )

        result = return_samples_with_fixed_QEPG(self._QEPG_graph, w, shots)
        arr = np.asarray(result)
        states = arr[:, :-1]
        observables = arr[:, -1]

        # Use BPOSD decoder
        predictions = self._decoder.decode_batch(states)
        # Handle the case where predictions might be 2D
        if predictions.ndim > 1:
            predictions = predictions.ravel()

        num_errors = np.count_nonzero(observables != predictions)
        return num_errors / shots

    # ------------------------------------------------------------------
    #  Sampling rate measurement - Override to use BPOSD
    # ------------------------------------------------------------------

    def measure_sample_rates(self) -> None:
        """
        Measure the sampling rate of the given circuit using BPOSD decoder.
        """
        assert self._QEPG_graph is not None, (
            "QEPG graph must be initialized before sampling"
        )
        assert self._decoder is not None, (
            "BPOSD decoder must be initialized before decoding"
        )

        wlist = [max(1, self._num_noise // 2)]
        slist = [100]  # Reduced from 1000 for faster BPOSD decoding

        start_time = time.perf_counter()
        print("Start time for sampling rate measurement:", start_time)
        detector_result, _obs = return_samples_many_weights_separate_obs_with_QEPG(
            self._QEPG_graph, wlist, slist
        )
        # Also measure decoding time
        _predictions = self._decoder.decode_batch(detector_result)
        end_time = time.perf_counter()
        print("End time for sampling rate measurement:", end_time)
        elapsed = end_time - start_time
        if elapsed <= 0:
            elapsed = 1e-6
        self._sampling_rate = 100.0 / elapsed  # Updated to match slist[0]=100
        self._remaining_time_budget -= elapsed

        print(
            "Elapsed time for sampling rate measurement: {:.6f} seconds".format(elapsed)
        )
        print(f"Measured sampling rate: {self._sampling_rate:.2f} shots/second")

    def profile_optimal_batch_size(
        self,
        test_weight: int | None = None,
        batch_sizes: List[int] | None = None,
        save_plot: str | None = None,
    ) -> tuple[int, dict]:
        """
        Profile sampling+decoding throughput at different batch sizes to find optimal.
        Uses BPOSD decoder.
        """
        assert self._QEPG_graph is not None, "QEPG graph must be initialized"
        assert self._decoder is not None, "BPOSD decoder must be initialized"

        if test_weight is None:
            test_weight = max(1, self._num_noise // 2)

        if batch_sizes is None:
            batch_sizes = [100, 500, 1000, 2000, 5000, 10000, 20000, 50000]

        print("=" * 70)
        print(f"Profiling batch sizes at weight w={test_weight}")
        print("=" * 70)
        print(
            f"{'Batch Size':>12} | {'Sample Time':>11} | {'Decode Time':>11} | {'Total Time':>11} | {'Throughput':>12}"
        )
        print("-" * 70)

        results: Dict[int, dict] = {}
        best_throughput = 0.0
        optimal_batch = batch_sizes[0]

        for batch in batch_sizes:
            wlist = [test_weight]
            slist = [batch]

            # Measure sampling time
            sample_start = time.perf_counter()
            detector_result, obsresult = (
                return_samples_many_weights_separate_obs_with_QEPG(
                    self._QEPG_graph, wlist, slist
                )
            )
            sample_end = time.perf_counter()
            sample_time = sample_end - sample_start

            # Measure decoding time with BPOSD
            decode_start = time.perf_counter()
            _predictions = self._decoder.decode_batch(detector_result)
            decode_end = time.perf_counter()
            decode_time = decode_end - decode_start

            total_time = sample_time + decode_time
            throughput = batch / total_time if total_time > 0 else 0

            results[batch] = {
                "throughput": throughput,
                "sample_time": sample_time,
                "decode_time": decode_time,
                "total_time": total_time,
            }

            print(
                f"{batch:>12,} | {sample_time:>11.4f}s | {decode_time:>11.4f}s | {total_time:>11.4f}s | {throughput:>12,.0f}/s"
            )

            if throughput > best_throughput:
                best_throughput = throughput
                optimal_batch = batch

            # Early stopping: if throughput drops significantly, we've passed the peak
            if throughput < best_throughput * 0.7 and batch > optimal_batch:
                print(f"  (Stopping early: throughput dropped below 70% of peak)")
                break

        print("=" * 70)
        print(
            f"Optimal batch size: {optimal_batch:,} (throughput: {best_throughput:,.0f} shots/s)"
        )

        # Save plot if requested
        if save_plot:
            self._plot_batch_profile(results, optimal_batch, save_plot)

        return optimal_batch, results

    # ------------------------------------------------------------------
    #  One multi-weight sampling step - Override to use BPOSD
    # ------------------------------------------------------------------

    def _sampling_step(self, wlist: List[int], slist: List[int]) -> float:
        """
        Perform one multi-weight sampling call and update subspace statistics.
        Returns the elapsed time for this step (seconds).
        Uses BPOSD decoder with capped batch size for performance.
        """
        assert self._QEPG_graph is not None, (
            "QEPG graph must be initialized before sampling"
        )
        assert self._decoder is not None, (
            "BPOSD decoder must be initialized before decoding"
        )

        if not wlist:
            return 0.0

        # Cap shots per weight at _MAX_BATCH_SIZE_LDPC for BPOSD performance
        slist_capped = [min(s, self._MAX_BATCH_SIZE_LDPC) for s in slist]
        total_shots = int(sum(slist_capped))
        print(
            f"Sampling weights and shots (capped at {self._MAX_BATCH_SIZE_LDPC}):",
            list(zip(wlist, slist_capped)),
        )
        print("  Total shots this step:", total_shots)

        start_time = time.perf_counter()
        detector_result, obsresult = return_samples_many_weights_separate_obs_with_QEPG(
            self._QEPG_graph, wlist, slist_capped
        )

        # Use BPOSD decoder
        predictions_result = self._decoder.decode_batch(detector_result)
        end_time = time.perf_counter()
        elapsed = end_time - start_time

        # Handle predictions shape
        if predictions_result.ndim > 1:
            predictions_result = predictions_result.ravel()

        begin_index = 0
        for w, shots in zip(wlist, slist_capped):
            end_index = begin_index + shots
            observables = np.asarray(obsresult[begin_index:end_index])
            predictions = np.asarray(predictions_result[begin_index:end_index]).ravel()

            num_errors = int(np.count_nonzero(observables != predictions))

            self._subspace_LE_count[w] = self._subspace_LE_count.get(w, 0) + num_errors
            self._subspace_sample_used[w] = self._subspace_sample_used.get(w, 0) + shots
            self._estimated_subspaceLER[w] = (
                self._subspace_LE_count[w] / self._subspace_sample_used[w]
            )

            begin_index = end_index

        print("  PL values:")
        for w in wlist:
            pl = self._estimated_subspaceLER.get(w, 0)
            le = self._subspace_LE_count.get(w, 0)
            samples = self._subspace_sample_used.get(w, 0)
            print(f"    w={w}: PL={pl:.4f}, LE_count={le}, samples={samples}")

        if elapsed > 0 and total_shots > 0:
            inst_rate = total_shots / elapsed
            if self._sampling_rate <= 0:
                self._sampling_rate = inst_rate
            else:
                # Exponential moving average
                alpha = 0.5
                self._sampling_rate = (
                    alpha * inst_rate + (1.0 - alpha) * self._sampling_rate
                )

        print(f"  Step elapsed: {elapsed:.6f} s")
        print(f"  Updated sampling rate: {self._sampling_rate:.2f} shots/second")
        return elapsed

    # ------------------------------------------------------------------
    #  BPOSD-specific methods
    # ------------------------------------------------------------------

    def get_decoder_params(self) -> Dict:
        """Get the current BPOSD decoder parameters."""
        return {
            "max_bp_iters": self._max_bp_iters,
            "bp_method": self._bp_method,
            "osd_method": self._osd_method,
            "osd_order": self._osd_order,
        }

    def set_decoder_params(
        self,
        max_bp_iters: int | None = None,
        bp_method: str | None = None,
        osd_method: str | None = None,
        osd_order: int | None = None,
    ) -> None:
        """
        Update BPOSD decoder parameters.
        Note: This will reinitialize the decoder if it has already been created.
        """
        if max_bp_iters is not None:
            self._max_bp_iters = max_bp_iters
        if bp_method is not None:
            self._bp_method = bp_method
        if osd_method is not None:
            self._osd_method = osd_method
        if osd_order is not None:
            self._osd_order = osd_order

        # Reinitialize decoder if it exists
        if self._decoder is not None and self._detector_error_model is not None:
            self._decoder = BPOSD(
                self._detector_error_model,
                max_bp_iters=self._max_bp_iters,
                bp_method=self._bp_method,
                osd_method=self._osd_method,
                osd_order=self._osd_order,
            )
