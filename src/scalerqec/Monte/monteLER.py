"""Monte Carlo logical error rate estimation using stim circuits.

This module implements adaptive Monte Carlo sampling strategies for
estimating the logical error rate (LER) of quantum error-correcting codes
described as stim circuits. It supports three sampling backends:

* **Stim detector sampling** with pymatching decoding (default).
* **QEPG-accelerated sampling** via a compiled C++ graph that injects
  errors and produces detector/observable outcomes without re-parsing
  the circuit each time.
* **Sinter distributed sampling** that fans out across CPU cores.

Adaptive batching is used throughout: the batch size starts at
:data:`SAMPLE_GAP_INITIAL` and is scaled up when few logical errors have
been observed, capped at :data:`MAX_SAMPLE_GAP`. Sampling terminates once
a configurable minimum number of logical error events is reached or the
sample/time budget is exhausted.
"""

from __future__ import annotations
from typing import Optional
from ..Clifford.clifford import CliffordCircuit
from ..Clifford.stimparser import rewrite_stim_code
import pymatching
import stim
import time
import os
from contextlib import redirect_stdout
import numpy as np

from ..qepg import compile_QEPG, return_samples_Monte_separate_obs_with_QEPG, QEPGGraph
from ..QEC.noisemodel import NoiseModel
from ..QEC.qeccircuit import StabCode
from .noise_model_parser import extract_noise_model, NonuniformNoiseModel

import sinter


def count_logical_errors(circuit: stim.Circuit, num_shots: int) -> int:
    """Sample a stim circuit and count logical errors using pymatching.

    The function compiles a detector sampler from the circuit, draws
    ``num_shots`` samples, builds a pymatching decoder from the detector
    error model, and counts shots where the decoded observable outcome
    differs from the actual observable flip.

    Args:
        circuit: A stim circuit that includes noise operations and
            observable/detector annotations.
        num_shots: Number of Monte Carlo shots to sample.

    Returns:
        The number of shots in which at least one observable was
        incorrectly decoded (i.e., logical error count).
    """
    # Sample the circuit.
    sampler = circuit.compile_detector_sampler()
    detection_events, observable_flips = sampler.sample(
        num_shots, separate_observables=True
    )

    # Configure a decoder using the circuit.
    detector_error_model = circuit.detector_error_model(decompose_errors=False)
    matcher = pymatching.Matching.from_detector_error_model(detector_error_model)

    # Run the decoder.
    predictions = matcher.decode_batch(detection_events)

    # Count the mistakes.
    num_errors = 0
    for shot in range(num_shots):
        actual_for_shot = observable_flips[shot]
        predicted_for_shot = predictions[shot]
        if not np.array_equal(actual_for_shot, predicted_for_shot):
            num_errors += 1
    return num_errors


def format_with_uncertainty(value: float, std: float) -> str:
    """Format a value with its uncertainty in scientific notation.

    Produces a human-readable string of the form
    ``1.23(+0.45)*10^k`` where *k* is chosen so that the coefficient
    is in [1, 10). When ``value`` is zero the output is ``0(+<std>)``.

    Args:
        value: The central value to format.
        std: The standard deviation (uncertainty) associated with
            ``value``.

    Returns:
        A string encoding both the value and its uncertainty in
        scientific notation, e.g. ``"1.23(+0.45)*10^-3"``.
    """
    if value == 0:
        return f"0(+{std:.2e})"
    exponent = int(np.floor(np.log10(abs(value))))
    coeff = value / (10**exponent)
    std_coeff = std / (10**exponent)
    return f"{coeff:.2f}(+{std_coeff:.2f})*10^{exponent}"


SAMPLE_GAP_INITIAL: int = 100
"""Initial batch size (number of shots) for the first sampling round."""

MAX_SAMPLE_GAP: int = 500000
"""Maximum batch size that adaptive batching is allowed to grow to."""


class MonteLERcalc:
    """Adaptive Monte Carlo estimator for logical error rates.

    This class wraps multiple sampling strategies (stim detector sampling,
    QEPG-accelerated sampling, and sinter distributed sampling) behind a
    common adaptive-batching loop. The loop starts with a small batch
    (:data:`SAMPLE_GAP_INITIAL`), decodes each batch with pymatching, and
    increases the batch size when few logical errors have been observed.
    Sampling stops when:

    * At least ``MIN_NUM_LE_EVENT`` logical errors have been recorded, **or**
    * The cumulative number of shots exceeds ``samplebudget``, **or**
    * (for time-budgeted runs) the wall-clock time exceeds ``time_budget``.

    When multiple averaging runs are requested via the ``repeat`` parameter,
    the final LER and sample count are averaged across runs.

    Attributes:
        _num_LER: Running count of logical errors in the current run.
        _sample_used: Total number of shots consumed (averaged across
            repeats after completion).
        _sample_needed: Estimated total samples needed to reach
            ``MIN_NUM_LE_EVENT`` errors, or -1 if no errors were seen.
        _uncertainty: Standard error of the LER estimate (set by
            :meth:`calculate_standard_error`).
        _estimated_LER: Current point estimate of the logical error rate.
        _samplebudget: Maximum number of shots per run.
        _min_num_ke_event: Minimum number of logical error events required
            before stopping.
        _QEPG: Compiled QEPG graph for accelerated sampling, or ``None``.
        _time_budget: Time budget in seconds for time-limited runs.
    """

    def __init__(
        self,
        time_budget: int = 10,
        samplebudget: int = 100000,
        MIN_NUM_LE_EVENT: int = 20,
    ) -> None:
        """Initialise the Monte Carlo LER calculator.

        Args:
            time_budget: Wall-clock time budget in seconds for
                time-limited sampling via
                :meth:`calculate_LER_with_time_budget`.
            samplebudget: Maximum number of Monte Carlo shots per
                averaging run before sampling is terminated.
            MIN_NUM_LE_EVENT: Minimum number of logical error events to
                collect before declaring convergence and stopping.
        """
        self._num_LER: int = 0
        self._sample_used: float = 0.0
        self._sample_needed: int = 0
        self._uncertainty: float = 0.0
        self._estimated_LER: float = 0.0
        self._samplebudget: int = samplebudget
        self._min_num_ke_event: int = MIN_NUM_LE_EVENT
        self._QEPG: Optional[QEPGGraph] = None
        self._time_budget: int = time_budget

    def calculate_LER_from_StabCode(
        self, qeccirc: StabCode, noise_model: NoiseModel, repeat: int = 1
    ) -> None:
        """Estimate the LER from a StabCode object using QEPG acceleration.

        Constructs an intermediate representation (IR) from the stabiliser
        code, compiles it to a stim circuit, injects noise according to
        ``noise_model``, and then compiles a QEPG graph for fast Monte
        Carlo sampling. The QEPG backend is compiled once and reused
        across all shots and repeats, avoiding repeated circuit parsing.

        Results (mean and standard deviation of LER, error count, sample
        count, and elapsed time) are printed to stdout.

        Args:
            qeccirc: A stabiliser code object. Its
                ``construct_IR_standard_scheme`` and
                ``compile_stim_circuit_from_IR_standard`` methods are
                called to produce the stim circuit.
            noise_model: Noise model whose ``error_rate`` is used for
                QEPG sampling and whose
                ``reconstruct_clifford_circuit`` method injects noise
                into the circuit.
            repeat: Number of independent averaging runs. The reported
                LER is the mean across all runs.

        Side Effects:
            Updates ``_estimated_LER``, ``_sample_used``, and ``_QEPG``
            on the instance. Prints timing and LER statistics to stdout.
        """
        qeccirc.construct_IR_standard_scheme()
        qeccirc.compile_stim_circuit_from_IR_standard()

        noisy_circuit = noise_model.reconstruct_clifford_circuit(qeccirc.circuit)
        stim_circuit = noisy_circuit.stimcircuit
        self._QEPG = compile_QEPG(str(stim_circuit))

        detector_error_model = stim_circuit.detector_error_model(decompose_errors=False)
        matcher = pymatching.Matching.from_detector_error_model(detector_error_model)

        error_rate = noise_model.error_rate
        Ler_list: list[float] = []
        samples_list: list[float] = []
        time_list: list[float] = []
        ler_count_list: list[int] = []
        for i in range(repeat):
            start = time.perf_counter()
            ler_count = 0
            sampleused = 0

            detector_result, obsresult = return_samples_Monte_separate_obs_with_QEPG(
                self._QEPG, error_rate, SAMPLE_GAP_INITIAL
            )
            predictions_result = matcher.decode_batch(detector_result)
            observables = np.asarray(obsresult).ravel()  # (shots,)
            predictions = np.asarray(predictions_result).ravel()
            num_errors = np.count_nonzero(observables != predictions)
            # 3. count mismatches in vectorised form ---------------------------------

            ler_count += num_errors
            sampleused += SAMPLE_GAP_INITIAL
            while (
                ler_count < self._min_num_ke_event and sampleused < self._samplebudget
            ):
                if ler_count == 0:
                    current_sample_gap = sampleused * 10
                    current_sample_gap = min(current_sample_gap, MAX_SAMPLE_GAP)
                else:
                    current_sample_gap = min(
                        int(self._min_num_ke_event / ler_count) * sampleused,
                        MAX_SAMPLE_GAP,
                    )

                detector_result, obsresult = (
                    return_samples_Monte_separate_obs_with_QEPG(
                        self._QEPG, error_rate, current_sample_gap
                    )
                )
                predictions_result = matcher.decode_batch(detector_result)
                observables = np.asarray(obsresult).ravel()  # (shots,)
                predictions = np.asarray(predictions_result).ravel()
                num_errors = np.count_nonzero(observables != predictions)
                ler_count += num_errors
                sampleused += current_sample_gap

            ler_count_list.append(ler_count)
            Ler_list.append(ler_count / sampleused)
            samples_list.append(sampleused)
            elapsed = time.perf_counter() - start
            time_list.append(elapsed)

        ler_count_average = float(np.mean(ler_count_list))
        # print("Average number of logical errors: ", ler_count_average)
        std_ler_count = float(np.std(ler_count_list))
        self._estimated_LER = float(np.mean(Ler_list))
        self._sample_used = float(np.mean(samples_list))
        """
        Standard deviation
        """
        std_ler = float(np.std(Ler_list))
        std_sample = float(np.std(samples_list))
        # self.calculate_standard_error()
        time_mean = float(np.mean(time_list))
        time_std = float(np.std(time_list))
        print("Time(STIM): ", format_with_uncertainty(time_mean, time_std))
        print("PL(STIM): ", format_with_uncertainty(self._estimated_LER, std_ler))
        print(
            "Nerror(STIM): ", format_with_uncertainty(ler_count_average, std_ler_count)
        )
        print("Sample(STIM): ", format_with_uncertainty(self._sample_used, std_sample))

    def calculate_LER_from_stim_circuit(
        self,
        stim_circuit_str: str,
        samplebudget: int = 100000,
        repeat: int = 1,
    ) -> float:
        """Estimate the LER of any Stim circuit with non-uniform noise.

        This method accepts a raw Stim circuit string that may contain
        arbitrary noise channels (DEPOLARIZE1 with varying rates, X_ERROR,
        Y_ERROR, Z_ERROR, DEPOLARIZE2) and estimates the logical error
        rate using QEPG-accelerated non-uniform Monte Carlo sampling.

        The method:
        1. Normalizes the circuit (strips noise) for QEPG compilation.
        2. Extracts the per-noise-source probability model from the
           original circuit.
        3. Compiles a Python QEPG from the normalized circuit.
        4. Builds a pymatching decoder from the original circuit's
           detector error model.
        5. Runs adaptive batching with non-uniform sampling.

        Args:
            stim_circuit_str: Raw Stim circuit string with noise directives.
            samplebudget: Maximum number of Monte Carlo shots per run.
            repeat: Number of independent averaging runs.

        Returns:
            The estimated logical error rate (averaged over ``repeat`` runs).
        """
        from ..Clifford.QEPGpython import QEPGpython

        self._samplebudget = samplebudget

        # 1. Normalize the circuit (strip noise for QEPG compilation)
        normalized = rewrite_stim_code(stim_circuit_str)

        # 2. Extract per-source noise model from the original circuit
        noise_model = extract_noise_model(stim_circuit_str)

        # 3. Compile QEPG from normalized circuit (Python + C++ backend)
        # The error_rate here is a placeholder — it populates DEPOLARIZE1 in the
        # compiled circuit structure, but actual sampling uses noise_model probabilities.
        circuit = CliffordCircuit(2)
        circuit.error_rate = 0.001
        circuit.compile_from_stim_circuit_str(normalized)

        graph = QEPGpython(circuit)
        graph.backword_graph_construction()

        # Attach C++ QEPG graph for accelerated sampling
        try:
            from scalerqec import qepg as qepg_cpp
            graph._cpp_graph = qepg_cpp.compile_QEPG(normalized)
        except Exception:
            graph._cpp_graph = None

        qepg_noise = circuit.totalnoise
        if noise_model.num_noise != qepg_noise:
            import warnings
            warnings.warn(
                f"Noise model has {noise_model.num_noise} sources but QEPG has "
                f"{qepg_noise}. This may indicate a bug in _GATE_NOISE_COUNT. "
                f"Resizing noise_probs to match QEPG.",
                stacklevel=2,
            )
            # Resize noise_probs to match QEPG if needed
            if noise_model.num_noise < qepg_noise:
                padded = np.zeros((qepg_noise, 3), dtype=np.float64)
                padded[:noise_model.num_noise] = noise_model.noise_probs
                noise_model.noise_probs = padded
                noise_model.num_noise = qepg_noise
            else:
                noise_model.noise_probs = noise_model.noise_probs[:qepg_noise]
                noise_model.num_noise = qepg_noise

        # 4. Build decoder from the ORIGINAL noisy circuit's DEM
        stim_circuit = stim.Circuit(stim_circuit_str)
        detector_error_model = stim_circuit.detector_error_model(
            decompose_errors=True
        )
        matcher = pymatching.Matching.from_detector_error_model(detector_error_model)

        # 5. Adaptive batching with non-uniform sampling
        Ler_list: list[float] = []
        samples_list: list[float] = []
        time_list: list[float] = []
        ler_count_list: list[int] = []

        for _ in range(repeat):
            start = time.perf_counter()
            ler_count = 0
            sampleused = 0

            # Initial batch
            det, obs = graph.sample_nonuniform_batch(
                noise_model.noise_probs,
                SAMPLE_GAP_INITIAL,
                correlated_pairs=noise_model.correlated_pairs or None,
            )
            predictions = matcher.decode_batch(det)
            observables = obs.ravel()
            predictions_flat = np.asarray(predictions).ravel()
            num_errors = np.count_nonzero(observables != predictions_flat)
            ler_count += num_errors
            sampleused += SAMPLE_GAP_INITIAL

            while ler_count < self._min_num_ke_event and sampleused < self._samplebudget:
                if ler_count == 0:
                    current_sample_gap = sampleused * 10
                    current_sample_gap = min(current_sample_gap, MAX_SAMPLE_GAP)
                else:
                    current_sample_gap = min(
                        int(self._min_num_ke_event / ler_count) * sampleused,
                        MAX_SAMPLE_GAP,
                    )

                det, obs = graph.sample_nonuniform_batch(
                    noise_model.noise_probs,
                    current_sample_gap,
                    correlated_pairs=noise_model.correlated_pairs or None,
                )
                predictions = matcher.decode_batch(det)
                observables = obs.ravel()
                predictions_flat = np.asarray(predictions).ravel()
                num_errors = np.count_nonzero(observables != predictions_flat)
                ler_count += num_errors
                sampleused += current_sample_gap

            ler_count_list.append(ler_count)
            Ler_list.append(ler_count / sampleused)
            samples_list.append(sampleused)
            elapsed = time.perf_counter() - start
            time_list.append(elapsed)

        ler_count_average = float(np.mean(ler_count_list))
        std_ler_count = float(np.std(ler_count_list))
        self._estimated_LER = float(np.mean(Ler_list))
        self._sample_used = float(np.mean(samples_list))
        std_ler = float(np.std(Ler_list))
        std_sample = float(np.std(samples_list))
        time_mean = float(np.mean(time_list))
        time_std = float(np.std(time_list))

        print("Time(NonUniform): ", format_with_uncertainty(time_mean, time_std))
        print("PL(NonUniform): ", format_with_uncertainty(self._estimated_LER, std_ler))
        print(
            "Nerror(NonUniform): ",
            format_with_uncertainty(ler_count_average, std_ler_count),
        )
        print(
            "Sample(NonUniform): ",
            format_with_uncertainty(self._sample_used, std_sample),
        )
        return self._estimated_LER

    def calculate_LER_from_my_random_sampler(
        self, samplebudget: int, filepath: str, pvalue: float, repeat: int = 1
    ) -> float:
        """Estimate the LER using the custom QEPG-based error injector.

        Reads a stim circuit from ``filepath``, compiles a QEPG graph for
        fast C++-backed error injection, and runs adaptive Monte Carlo
        sampling with pymatching decoding. The QEPG graph is compiled
        once and reused for all shots, providing significant speedup over
        repeated stim sampling.

        Results (mean and standard deviation of LER, error count, sample
        count, and elapsed time) are printed to stdout.

        Args:
            samplebudget: Maximum number of Monte Carlo shots per
                averaging run.
            filepath: Path to a stim circuit file (text format).
            pvalue: Physical error rate (depolarising probability)
                injected into every noise location.
            repeat: Number of independent averaging runs. The reported
                LER is the mean across all runs.

        Returns:
            The estimated logical error rate (averaged over ``repeat``
            runs).

        Side Effects:
            Updates ``_estimated_LER``, ``_sample_used``, ``_samplebudget``,
            and ``_QEPG`` on the instance. Prints statistics to stdout.
        """
        circuit = CliffordCircuit(2)
        circuit.error_rate = pvalue
        self._samplebudget = samplebudget

        stim_str = ""
        with open(filepath, "r", encoding="utf-8") as f:
            stim_str = f.read()
        self._QEPG = compile_QEPG(stim_str)

        stim_circuit = rewrite_stim_code(stim_str)
        circuit.compile_from_stim_circuit_str(stim_circuit)
        new_stim_circuit = circuit.stimcircuit

        detector_error_model = new_stim_circuit.detector_error_model(
            decompose_errors=False
        )
        matcher = pymatching.Matching.from_detector_error_model(detector_error_model)

        Ler_list: list[float] = []
        samples_list: list[float] = []
        time_list: list[float] = []
        ler_count_list: list[int] = []
        for i in range(repeat):
            ler_count = 0
            sampleused = 0
            start = time.time()

            detector_result, obsresult = return_samples_Monte_separate_obs_with_QEPG(
                self._QEPG, pvalue, SAMPLE_GAP_INITIAL
            )
            predictions_result = matcher.decode_batch(detector_result)
            observables = np.asarray(obsresult).ravel()  # (shots,)
            predictions = np.asarray(predictions_result).ravel()
            num_errors = np.count_nonzero(observables != predictions)
            # 3. count mismatches in vectorised form ---------------------------------

            ler_count += num_errors
            sampleused += SAMPLE_GAP_INITIAL
            while (
                ler_count < self._min_num_ke_event and sampleused < self._samplebudget
            ):
                if ler_count == 0:
                    current_sample_gap = sampleused * 10
                    current_sample_gap = min(current_sample_gap, MAX_SAMPLE_GAP)
                else:
                    current_sample_gap = min(
                        int(self._min_num_ke_event / ler_count) * sampleused,
                        MAX_SAMPLE_GAP,
                    )

                detector_result, obsresult = (
                    return_samples_Monte_separate_obs_with_QEPG(
                        self._QEPG, pvalue, current_sample_gap
                    )
                )
                predictions_result = matcher.decode_batch(detector_result)
                observables = np.asarray(obsresult).ravel()  # (shots,)
                predictions = np.asarray(predictions_result).ravel()
                num_errors = np.count_nonzero(observables != predictions)
                ler_count += num_errors
                sampleused += current_sample_gap

            ler_count_list.append(ler_count)
            Ler_list.append(ler_count / sampleused)
            samples_list.append(sampleused)
            elapsed = time.time() - start
            time_list.append(elapsed)

        ler_count_average = float(np.mean(ler_count_list))
        # print("Average number of logical errors: ", ler_count_average)
        std_ler_count = float(np.std(ler_count_list))

        self._estimated_LER = float(np.mean(Ler_list))
        self._sample_used = float(np.mean(samples_list))
        """
        Standard deviation
        """
        std_ler = float(np.std(Ler_list))
        std_sample = float(np.std(samples_list))
        # self.calculate_standard_error()
        time_mean = float(np.mean(time_list))
        time_std = float(np.std(time_list))

        print("Time(STIM): ", format_with_uncertainty(time_mean, time_std))
        print("PL(STIM): ", format_with_uncertainty(self._estimated_LER, std_ler))
        print(
            "Nerror(STIM): ", format_with_uncertainty(ler_count_average, std_ler_count)
        )
        print("Sample(STIM): ", format_with_uncertainty(self._sample_used, std_sample))
        return self._estimated_LER

    def calculate_LER_from_file_sinter(
        self, samplebudget: int, filepath: str, pvalue: float, repeat: int = 1
    ) -> float:
        """Estimate the LER using sinter for distributed sampling.

        Reads a stim circuit from ``filepath``, injects depolarising
        noise at rate ``pvalue``, and delegates sampling and decoding to
        ``sinter.collect`` which fans work across all available CPU
        cores. Sinter handles its own adaptive stopping (via
        ``max_errors``).

        Results (mean and standard deviation of LER, error count, sample
        count, and elapsed time) are printed to stdout.

        Args:
            samplebudget: Maximum number of Monte Carlo shots passed to
                sinter as ``max_shots``.
            filepath: Path to a stim circuit file (text format).
            pvalue: Physical error rate (depolarising probability)
                injected into every noise location.
            repeat: Number of independent averaging runs. The reported
                LER is the mean across all runs.

        Returns:
            The estimated logical error rate (averaged over ``repeat``
            runs).

        Side Effects:
            Updates ``_estimated_LER``, ``_sample_used``, ``_num_LER``,
            and ``_samplebudget`` on the instance. Prints statistics to
            stdout.
        """
        circuit = CliffordCircuit(2)
        circuit.error_rate = pvalue
        self._samplebudget = samplebudget

        stim_str = ""
        with open(filepath, "r", encoding="utf-8") as f:
            stim_str = f.read()

        stim_circuit = rewrite_stim_code(stim_str)
        circuit.compile_from_stim_circuit_str(stim_circuit)
        new_stim_circuit = circuit.stimcircuit

        Ler_list: list[float] = []
        samples_list: list[float] = []
        time_list: list[float] = []
        ler_count_list: list[int] = []
        for i in range(repeat):
            start = time.time()
            self._num_LER = 0
            self._sample_used = 0

            mytask = sinter.Task(
                circuit=new_stim_circuit,
                json_metadata={
                    "p": pvalue,
                    "d": 0,
                },
            )
            samples = sinter.collect(
                num_workers=os.cpu_count(),
                max_shots=samplebudget,
                max_errors=self._min_num_ke_event,
                tasks=[mytask],
                decoders=["pymatching"],
                count_observable_error_combos=False,
            )

            self._num_LER = int(samples[0].errors)
            ler_count_list.append(self._num_LER)
            self._sample_used = int(samples[0].shots)

            Ler_list.append(self._num_LER / self._sample_used)
            samples_list.append(self._sample_used)
            elapsed = time.time() - start
            time_list.append(elapsed)

        ler_count_average = float(np.mean(ler_count_list))
        # print("Average number of logical errors: ", ler_count_average)
        std_ler_count = float(np.std(ler_count_list))

        self._estimated_LER = float(np.mean(Ler_list))
        self._sample_used = float(np.mean(samples_list))
        """
        Standard deviation
        """
        std_ler = float(np.std(Ler_list))
        std_sample = float(np.std(samples_list))
        # self.calculate_standard_error()
        time_mean = float(np.mean(time_list))
        time_std = float(np.std(time_list))

        print("Time(STIM): ", format_with_uncertainty(time_mean, time_std))
        print("PL(STIM): ", format_with_uncertainty(self._estimated_LER, std_ler))
        print(
            "Nerror(STIM): ", format_with_uncertainty(ler_count_average, std_ler_count)
        )
        print("Sample(STIM): ", format_with_uncertainty(self._sample_used, std_sample))
        return self._estimated_LER

    def calculate_LER_from_file(
        self, samplebudget: int, filepath: str, pvalue: float, repeat: int = 1
    ) -> float:
        """Estimate the LER by sampling a stim circuit from a file.

        This is the main entry point for standard stim-based Monte Carlo
        estimation. The method reads a stim circuit from ``filepath``,
        injects depolarising noise at rate ``pvalue``, compiles a
        detector sampler, and runs adaptive batching with pymatching
        decoding until the minimum number of logical error events is
        reached or the sample budget is exhausted.

        Adaptive batching behaviour:

        * Starts with :data:`SAMPLE_GAP_INITIAL` shots.
        * If no errors have been seen, doubles the batch size each
          iteration (capped at :data:`MAX_SAMPLE_GAP`).
        * If some errors have been seen, estimates the remaining
          samples needed and uses that as the next batch size (capped
          at :data:`MAX_SAMPLE_GAP`).

        Results (mean and standard deviation of LER, error count, sample
        count, and elapsed time) are printed to stdout.

        Args:
            samplebudget: Maximum number of Monte Carlo shots per
                averaging run.
            filepath: Path to a stim circuit file (text format).
            pvalue: Physical error rate (depolarising probability)
                injected into every noise location.
            repeat: Number of independent averaging runs. The reported
                LER is the mean across all runs.

        Returns:
            The estimated logical error rate (averaged over ``repeat``
            runs).

        Side Effects:
            Updates ``_estimated_LER``, ``_sample_used``,
            ``_sample_needed``, ``_num_LER``, and ``_samplebudget`` on
            the instance. Prints statistics to stdout.
        """
        circuit = CliffordCircuit(2)
        circuit.error_rate = pvalue
        self._samplebudget = samplebudget

        stim_str = ""
        with open(filepath, "r", encoding="utf-8") as f:
            stim_str = f.read()

        stim_circuit = rewrite_stim_code(stim_str)
        circuit.compile_from_stim_circuit_str(stim_circuit)
        new_stim_circuit = circuit.stimcircuit

        sampler = new_stim_circuit.compile_detector_sampler()
        detector_error_model = new_stim_circuit.detector_error_model(
            decompose_errors=False
        )
        matcher = pymatching.Matching.from_detector_error_model(detector_error_model)

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
                        int(self._min_num_ke_event / self._num_LER) * self._sample_used,
                        MAX_SAMPLE_GAP,
                    )
                self._sample_used += current_sample_gap
                # self._num_LER+=count_logical_errors(new_stim_circuit,SAMPLE_GAP)

                detection_events, observable_flips = sampler.sample(
                    current_sample_gap, separate_observables=True
                )
                predictions = matcher.decode_batch(detection_events)
                # 3. count mismatches in vectorised form ---------------------------------
                num_errors = np.count_nonzero(observable_flips != predictions)
                self._num_LER += num_errors

                self._estimated_LER = self._num_LER / self._sample_used
                # self.calculate_standard_error()
                # print("Current LER: ", self._num_LER)
                # print("Current logical error rate: ", self._num_LER/self._sample_used)
                # print("Current stdandard error: ", self._uncertainty)
                # print("Current sample used: ", self._sample_used)

                if self._sample_used > self._samplebudget:
                    # print("Sample budget reached, stop sampling")
                    if self._num_LER > 0:
                        self._sample_needed = int(
                            self._sample_used * (100 / self._num_LER)
                        )
                    else:
                        self._sample_needed = -1
                    break
                self._sample_needed = self._sample_used
            ler_count_list.append(self._num_LER)
            Ler_list.append(self._estimated_LER)
            samples_list.append(self._sample_used)
            elapsed = time.time() - start
            time_list.append(elapsed)

        ler_count_average = float(np.mean(ler_count_list))
        # print("Average number of logical errors: ", ler_count_average)
        std_ler_count = float(np.std(ler_count_list))

        self._estimated_LER = float(np.mean(Ler_list))
        self._sample_used = float(np.mean(samples_list))
        """
        Standard deviation
        """
        std_ler = float(np.std(Ler_list))
        std_sample = float(np.std(samples_list))
        # self.calculate_standard_error()
        time_mean = float(np.mean(time_list))
        time_std = float(np.std(time_list))

        print("Time(STIM): ", format_with_uncertainty(time_mean, time_std))
        print("PL(STIM): ", format_with_uncertainty(self._estimated_LER, std_ler))
        print(
            "Nerror(STIM): ", format_with_uncertainty(ler_count_average, std_ler_count)
        )
        print("Sample(STIM): ", format_with_uncertainty(self._sample_used, std_sample))
        return self._estimated_LER

    def calculate_LER_with_time_budget(
        self, time_budget: float, filepath: str, pvalue: float
    ) -> dict:
        """
        Calculate the logical error rate with a strict time budget.
        Returns whatever results are available when time budget is exhausted.

        Args:
            time_budget: Time budget in seconds
            filepath: Path to the stim circuit file
            pvalue: Error rate (probability)

        Returns:
            Dictionary with:
                - ler: Estimated logical error rate
                - le_count: Number of logical errors observed
                - samples_used: Total number of samples used
                - time_elapsed: Actual time elapsed
                - budget_exhausted: True if stopped due to time budget
        """
        circuit = CliffordCircuit(2)
        circuit.error_rate = pvalue

        stim_str = ""
        with open(filepath, "r", encoding="utf-8") as f:
            stim_str = f.read()

        stim_circuit = rewrite_stim_code(stim_str)
        circuit.compile_from_stim_circuit_str(stim_circuit)
        new_stim_circuit = circuit.stimcircuit

        sampler = new_stim_circuit.compile_detector_sampler()
        detector_error_model = new_stim_circuit.detector_error_model(
            decompose_errors=False
        )
        matcher = pymatching.Matching.from_detector_error_model(detector_error_model)

        start_time = time.perf_counter()
        ler_count = 0
        samples_used = 0
        current_sample_gap = SAMPLE_GAP_INITIAL
        budget_exhausted = False

        while True:
            elapsed = time.perf_counter() - start_time
            remaining = time_budget - elapsed

            # Check if time budget is exhausted
            if remaining <= 0:
                budget_exhausted = True
                break

            # Check if we have enough LE events
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

            # Sample and decode
            detection_events, observable_flips = sampler.sample(
                current_sample_gap, separate_observables=True
            )
            predictions = matcher.decode_batch(detection_events)
            num_errors = np.count_nonzero(observable_flips != predictions)

            ler_count += num_errors
            samples_used += current_sample_gap

        elapsed = time.perf_counter() - start_time

        # Calculate LER
        if samples_used > 0:
            estimated_ler = ler_count / samples_used
        else:
            estimated_ler = 0.0

        # Store in instance variables
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

    def calculate_standard_error(self) -> float:
        """Compute the standard error of the current LER estimate.

        Uses the binomial standard-error formula:

            SE = sqrt(p * (1 - p)) / N

        where *p* is the estimated LER and *N* is the total number of
        shots consumed. The result is stored in ``_uncertainty``.

        Returns:
            The standard error of the LER estimate.

        Side Effects:
            Updates ``_estimated_LER`` (recomputed from ``_num_LER``
            and ``_sample_used``) and ``_uncertainty`` on the instance.
        """
        self._estimated_LER = self._num_LER / self._sample_used
        self._uncertainty = float(
            np.sqrt(self._estimated_LER * (1 - self._estimated_LER)) / self._sample_used
        )
        return self._uncertainty

    def get_sample_used(self) -> float:
        """Return the total number of Monte Carlo shots consumed.

        Returns:
            The cumulative sample count (averaged across repeats if
            multiple runs were performed).
        """
        return self._sample_used


if __name__ == "__main__":
    base_dir = "C:/Users/username/Documents/Sampling/stimprograms/"
    result_dir = "C:/Users/username/Documents/Sampling/"

    p = 0.0005
    code_type = "hexagon"
    rel = "hexagon/hexagon"
    dlist = [11]
    for d in dlist:
        stim_path = base_dir + rel + str(d)
        # 3) build your output filename:
        out_fname = (
            result_dir + str(p) + "-" + str(code_type) + str(d) + "-resultMonte.txt"
        )  # e.g. "surface3-result.txt"
        # 4) redirect prints for just this file:
        with open(out_fname, "w") as outf, redirect_stdout(outf):
            print(f"---- Processing {stim_path} ----")

            calculator = MonteLERcalc(20)
            # pass the string path into your function:
            ler = calculator.calculate_LER_from_my_random_sampler(
                500000000, str(stim_path), p, 5
            )

    p = 0.001
    code_type = "hexagon"
    rel = "hexagon/hexagon"
    dlist = [15]
    for d in dlist:
        stim_path = base_dir + rel + str(d)
        # 3) build your output filename:
        out_fname = (
            result_dir + str(p) + "-" + str(code_type) + str(d) + "-resultMonte.txt"
        )  # e.g. "surface3-result.txt"
        # 4) redirect prints for just this file:
        with open(out_fname, "w") as outf, redirect_stdout(outf):
            print(f"---- Processing {stim_path} ----")

            calculator = MonteLERcalc(10)
            # pass the string path into your function:
            ler = calculator.calculate_LER_from_my_random_sampler(
                500000000, str(stim_path), p, 5
            )
