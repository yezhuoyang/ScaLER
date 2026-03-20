"""
MonteLDPC: Monte Carlo LER estimation for LDPC codes using BP+OSD decoder

This module provides Monte Carlo simulation for LDPC codes where PyMatching
cannot be used. It uses the BPOSD decoder from stimbposd.
"""

from __future__ import annotations
from typing import Optional
import numpy as np
import time

from ..Clifford.clifford import CliffordCircuit
from ..Clifford.stimparser import rewrite_stim_code
from ..QEC.noisemodel import SIDNoiseModel

try:
    from stimbposd import BPOSD
except ImportError:
    BPOSD = None


def format_with_uncertainty(value: float, std: float) -> str:
    """Format a value and its standard deviation."""
    if value == 0:
        return f"0(+{std:.2e})"
    exponent = int(np.floor(np.log10(abs(value))))
    coeff = value / (10**exponent)
    std_coeff = std / (10**exponent)
    return f"{coeff:.2f}(+{std_coeff:.2f})*10^{exponent}"


# Sampling parameters - smaller batches for slower BPOSD decoder
SAMPLE_GAP_INITIAL = 100
MAX_SAMPLE_GAP = 10000  # Smaller than regular Monte Carlo due to BPOSD speed


class MonteLDPC:
    """
    Monte Carlo LER calculator for LDPC codes using BPOSD decoder.

    This class provides Monte Carlo simulation for LDPC codes where
    minimum-weight perfect matching (PyMatching) cannot be used.
    """

    def __init__(
        self,
        time_budget: int = 10,
        samplebudget: int = 100000,
        MIN_NUM_LE_EVENT: int = 100,
        # BPOSD parameters
        max_bp_iters: int = 20,
        bp_method: str = "product_sum",
        osd_method: str = "osd0",
        osd_order: int = 0,
        decoder=None,
    ) -> None:
        """
        Initialize the MonteLDPC calculator.

        Args:
            time_budget: Time budget in seconds
            samplebudget: Maximum number of samples
            MIN_NUM_LE_EVENT: Minimum number of logical errors to collect
            max_bp_iters: Maximum BP iterations for BPOSD
            bp_method: BP method ('product_sum' or 'min_sum')
            osd_method: OSD method ('osd0', 'osd_e', 'osd_cs')
            osd_order: OSD order (higher for better accuracy)
            decoder: Any object with a ``decode_batch`` method. If ``None``,
                a BPOSD decoder is built automatically.
        """
        self._num_LER: int = 0
        self._sample_used: float = 0.0
        self._sample_needed: int = 0
        self._uncertainty: float = 0.0
        self._estimated_LER: float = 0.0
        self._samplebudget: int = samplebudget
        self._min_num_ke_event: int = MIN_NUM_LE_EVENT
        self._time_budget: int = time_budget

        # BPOSD parameters
        self._max_bp_iters: int = max_bp_iters
        self._bp_method: str = bp_method
        self._osd_method: str = osd_method
        self._osd_order: int = osd_order

        if decoder is None and BPOSD is None:
            raise ImportError(
                "stimbposd is required for MonteLDPC. "
                "Install it with: pip install stimbposd"
            )
        self._decoder = decoder

    def calculate_LER_from_file(
        self, samplebudget: int, filepath: str, pvalue: float, repeat: int = 1
    ) -> float:
        """
        Calculate LER from a stim circuit file using BPOSD decoder.

        Args:
            samplebudget: Maximum number of samples
            filepath: Path to the stim circuit file
            pvalue: Physical error rate
            repeat: Number of repetitions

        Returns:
            Estimated logical error rate
        """
        circuit = CliffordCircuit(2)
        self._samplebudget = samplebudget

        # Read and compile circuit
        with open(filepath, "r", encoding="utf-8") as f:
            stim_str = f.read()

        stim_circuit = rewrite_stim_code(stim_str)
        circuit.compile_from_stim_circuit_str(stim_circuit)
        noisy_stim = SIDNoiseModel(pvalue).inject_noise(circuit.stimcircuit)

        # Create sampler and decoder
        sampler = noisy_stim.compile_detector_sampler()
        detector_error_model = noisy_stim.detector_error_model(decompose_errors=False)

        if self._decoder is None:
            self._decoder = BPOSD(
                detector_error_model,
                max_bp_iters=self._max_bp_iters,
                bp_method=self._bp_method,
                osd_method=self._osd_method,
                osd_order=self._osd_order,
            )

        Ler_list: list[float] = []
        samples_list: list[float] = []
        time_list: list[float] = []
        ler_count_list: list[int] = []

        for i in range(repeat):
            start = time.time()
            self._num_LER = 0
            self._sample_used = 0
            current_sample_gap = SAMPLE_GAP_INITIAL

            while self._num_LER < self._min_num_ke_event:
                if self._num_LER == 0 and self._sample_used > 0:
                    current_sample_gap *= 2
                    current_sample_gap = min(current_sample_gap, MAX_SAMPLE_GAP)
                elif self._num_LER > 0:
                    current_sample_gap = min(
                        int(self._min_num_ke_event / self._num_LER)
                        * int(self._sample_used),
                        MAX_SAMPLE_GAP,
                    )

                self._sample_used += current_sample_gap

                # Sample and decode with BPOSD
                detection_events, observable_flips = sampler.sample(
                    current_sample_gap, separate_observables=True
                )
                predictions = self._decoder.decode_batch(detection_events)

                # Handle array dimensions
                if predictions.ndim > 1:
                    predictions = predictions.ravel()
                if observable_flips.ndim > 1:
                    observable_flips = observable_flips.ravel()

                num_errors = np.count_nonzero(observable_flips != predictions)
                self._num_LER += num_errors

                self._estimated_LER = self._num_LER / self._sample_used

                if self._sample_used > self._samplebudget:
                    if self._num_LER > 0:
                        self._sample_needed = int(
                            self._sample_used * (100 / self._num_LER)
                        )
                    else:
                        self._sample_needed = -1
                    break

                self._sample_needed = int(self._sample_used)

            ler_count_list.append(self._num_LER)
            Ler_list.append(self._estimated_LER)
            samples_list.append(self._sample_used)
            elapsed = time.time() - start
            time_list.append(elapsed)

        ler_count_average = float(np.mean(ler_count_list))
        std_ler_count = float(np.std(ler_count_list))

        self._estimated_LER = float(np.mean(Ler_list))
        self._sample_used = float(np.mean(samples_list))

        std_ler = float(np.std(Ler_list))
        std_sample = float(np.std(samples_list))
        time_mean = float(np.mean(time_list))
        time_std = float(np.std(time_list))

        print("Time(BPOSD): ", format_with_uncertainty(time_mean, time_std))
        print("PL(BPOSD): ", format_with_uncertainty(self._estimated_LER, std_ler))
        print(
            "Nerror(BPOSD): ", format_with_uncertainty(ler_count_average, std_ler_count)
        )
        print("Sample(BPOSD): ", format_with_uncertainty(self._sample_used, std_sample))

        return self._estimated_LER

    def calculate_LER_with_time_budget(
        self, time_budget: float, filepath: str, pvalue: float
    ) -> dict:
        """
        Calculate LER with a strict time budget.

        Args:
            time_budget: Time budget in seconds
            filepath: Path to the stim circuit file
            pvalue: Physical error rate

        Returns:
            Dictionary with ler, le_count, samples_used, time_elapsed, budget_exhausted
        """
        circuit = CliffordCircuit(2)

        with open(filepath, "r", encoding="utf-8") as f:
            stim_str = f.read()

        stim_circuit = rewrite_stim_code(stim_str)
        circuit.compile_from_stim_circuit_str(stim_circuit)
        noisy_stim = SIDNoiseModel(pvalue).inject_noise(circuit.stimcircuit)

        sampler = noisy_stim.compile_detector_sampler()
        detector_error_model = noisy_stim.detector_error_model(decompose_errors=False)

        if self._decoder is None:
            self._decoder = BPOSD(
                detector_error_model,
                max_bp_iters=self._max_bp_iters,
                bp_method=self._bp_method,
                osd_method=self._osd_method,
                osd_order=self._osd_order,
            )

        start_time = time.perf_counter()
        ler_count = 0
        samples_used = 0
        current_sample_gap = SAMPLE_GAP_INITIAL
        budget_exhausted = False

        while True:
            elapsed = time.perf_counter() - start_time
            remaining = time_budget - elapsed

            if remaining <= 0:
                budget_exhausted = True
                break

            if ler_count >= self._min_num_ke_event:
                break

            # Adaptive sample gap
            if ler_count == 0 and samples_used > 0:
                current_sample_gap *= 2
                current_sample_gap = min(current_sample_gap, MAX_SAMPLE_GAP)
            elif ler_count > 0:
                current_sample_gap = min(
                    int(self._min_num_ke_event / ler_count) * samples_used,
                    MAX_SAMPLE_GAP,
                )

            # Sample and decode with BPOSD
            detection_events, observable_flips = sampler.sample(
                current_sample_gap, separate_observables=True
            )
            predictions = self._decoder.decode_batch(detection_events)

            if predictions.ndim > 1:
                predictions = predictions.ravel()
            if observable_flips.ndim > 1:
                observable_flips = observable_flips.ravel()

            num_errors = np.count_nonzero(observable_flips != predictions)

            ler_count += num_errors
            samples_used += current_sample_gap

        elapsed = time.perf_counter() - start_time

        if samples_used > 0:
            estimated_ler = ler_count / samples_used
        else:
            estimated_ler = 0.0

        self._num_LER = ler_count
        self._sample_used = samples_used
        self._estimated_LER = estimated_ler

        return {
            "ler": estimated_ler,
            "le_count": ler_count,
            "samples_used": samples_used,
            "time_elapsed": elapsed,
            "budget_exhausted": budget_exhausted,
        }

    def get_sample_used(self) -> float:
        return self._sample_used

    def get_estimated_LER(self) -> float:
        return self._estimated_LER


if __name__ == "__main__":
    # Test with a BB code
    filepath = "stimprograms/ldpc/bbcode_72_12_6_rounds18"
    p = 0.0005

    print(f"Testing MonteLDPC with {filepath} at p={p}")
    calculator = MonteLDPC(
        MIN_NUM_LE_EVENT=50,
        max_bp_iters=20,
        osd_method="osd0",
    )
    ler = calculator.calculate_LER_from_file(1000000, filepath, p, repeat=1)
    print(f"Estimated LER: {ler:.6e}")
