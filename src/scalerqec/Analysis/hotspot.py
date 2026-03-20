"""Hotspot analysis: conditional error attribution by noise category.

Given a noise label map and a **user-chosen decoder**, this module computes
conditional probabilities like::

    P(category fired | logical error)

which reveals which noise categories contribute most to logical failures.

The decoder is **pluggable** — any object with a ``decode_batch(detector_outcomes)``
method works (pymatching, BPOSD, custom decoders).

The analysis pipeline has three phases:

1. **Sample** (C++): generate detector/observable outcomes + per-shot category bitmask
2. **Decode** (user's decoder): ``predictions = decoder.decode_batch(detectors)``
3. **Aggregate** (C++): compute conditional probabilities from bitmask + logical errors

Typical usage::

    from scalerqec.Analysis.hotspot import HotspotAnalyzer

    analyzer = HotspotAnalyzer(qepg_graph, label_map, decoder=matcher)
    result = analyzer.analyze(noise_probs, num_shots=1_000_000)
    analyzer.print_report(result)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Protocol

import numpy as np


class Decoder(Protocol):
    """Protocol for decoders usable with hotspot analysis.

    Any object with a ``decode_batch`` method accepting a 2D uint8/bool
    array of detector outcomes and returning predicted observable outcomes.
    Both pymatching.Matching and ldpc.bposd_decoder satisfy this.
    """

    def decode_batch(self, detector_outcomes: np.ndarray) -> np.ndarray:
        ...


@dataclass
class CategoryReport:
    """Conditional probability report for a single noise category.

    Attributes:
        category: The noise label category name.
        count: Number of noise sources with this label.
        p_fired: P(any source in category fired) -- marginal.
        p_fired_given_error: P(any source in category fired | logical error).
        p_error_given_fired: P(logical error | any source in category fired).
        lift: ``p_fired_given_error / p_fired``.
        num_logical_errors: Number of shots with logical errors.
    """

    category: str
    count: int
    p_fired: float
    p_fired_given_error: float
    p_error_given_fired: float
    lift: float
    num_logical_errors: int


@dataclass
class HotspotResult:
    """Complete hotspot analysis result.

    Attributes:
        num_shots: Total Monte Carlo shots.
        num_logical_errors: Decoded logical errors (using the user's decoder).
        ler: Logical error rate = num_logical_errors / num_shots.
        category_reports: Per-category conditional probability reports.
        config_reports: Multi-error configuration breakdown.
    """

    num_shots: int
    num_logical_errors: int
    ler: float
    category_reports: list[CategoryReport]
    config_reports: list[dict[str, Any]]


class HotspotAnalyzer:
    """Decoder-agnostic hotspot analysis with C++ acceleration.

    Three-phase pipeline:

    1. **C++ labeled_sample**: non-uniform noise sampling via QEPG with
       per-shot category bitmask tracking (near-zero overhead).
    2. **User's decoder**: ``decoder.decode_batch(detector_outcomes)``
       to get predicted observables.
    3. **C++ hotspot_aggregate**: compute P(category fired | logical error),
       lift, and multi-error configuration breakdown.

    The decoder is pluggable — pass any object with ``decode_batch()``.
    This guarantees the LER reported by hotspot analysis matches the LER
    you would get from standard sampling + decoding.

    Args:
        qepg_graph: Pre-compiled QEPG graph (from ``qepg.compile_QEPG()``).
        label_map: C++ NoiseLabelMap (from ``qepg.auto_label()``).
        decoder: Any decoder with ``decode_batch(det_outcomes) -> predictions``.
            Examples: ``pymatching.Matching``, ``ldpc.bposd_decoder``.
    """

    def __init__(
        self,
        qepg_graph: Any,
        label_map: Any,
        decoder: Decoder,
    ) -> None:
        self._graph = qepg_graph
        self._label_map = label_map
        self._decoder = decoder

    def analyze(
        self,
        noise_probs: np.ndarray,
        num_shots: int,
        corr_sources_a: np.ndarray | None = None,
        corr_sources_b: np.ndarray | None = None,
        corr_probs: np.ndarray | None = None,
        min_config_count: int = 5,
    ) -> HotspotResult:
        """Run the full three-phase hotspot analysis pipeline.

        Phase 1 (C++): Sample noise, compute detector/observable outcomes
        and per-shot category bitmasks.

        Phase 2 (decoder): Decode detector outcomes to get predicted
        observables, then compute logical errors.

        Phase 3 (C++): Aggregate hotspot statistics from bitmasks +
        decoded logical error mask.

        Args:
            noise_probs: Shape ``(num_noise, 3)`` with ``[px, py, pz]``.
            num_shots: Number of Monte Carlo shots.
            corr_sources_a: Source A indices for DEPOLARIZE2 pairs.
            corr_sources_b: Source B indices for DEPOLARIZE2 pairs.
            corr_probs: Probabilities for DEPOLARIZE2 pairs.
            min_config_count: Minimum count for config report entries.

        Returns:
            A :class:`HotspotResult` with LER, category reports, and
            configuration reports. The LER is computed using the decoder
            you provided — it will match standard sampling + decoding.
        """
        from scalerqec import qepg as qepg_cpp

        # Default empty correlated pair arrays
        if corr_sources_a is None:
            corr_sources_a = np.array([], dtype=np.uint64)
        if corr_sources_b is None:
            corr_sources_b = np.array([], dtype=np.uint64)
        if corr_probs is None:
            corr_probs = np.array([], dtype=np.float64)

        noise_probs = np.ascontiguousarray(noise_probs, dtype=np.float64)

        # Phase 1: C++ labeled sampling
        det_outcomes, obs_outcomes, bitmasks = qepg_cpp.labeled_sample(
            self._graph,
            self._label_map,
            noise_probs,
            corr_sources_a,
            corr_sources_b,
            corr_probs,
            num_shots,
        )

        # Phase 2: Decode with user's decoder
        predictions = self._decoder.decode_batch(det_outcomes)
        predictions = np.asarray(predictions).ravel().astype(np.uint8)
        obs_flat = np.asarray(obs_outcomes).ravel().astype(np.uint8)
        logical_errors = (obs_flat != predictions).astype(np.uint8)

        num_logical_errors = int(np.count_nonzero(logical_errors))
        ler = num_logical_errors / num_shots if num_shots > 0 else 0.0

        # Phase 3: C++ aggregation
        cpp_result = qepg_cpp.hotspot_aggregate(
            bitmasks,
            logical_errors,
            self._label_map,
            min_config_count,
        )

        # Convert C++ result to Python dataclasses
        category_reports = []
        for cr in cpp_result.category_reports:
            category_reports.append(
                CategoryReport(
                    category=cr.name,
                    count=cr.count,
                    p_fired=cr.p_fired,
                    p_fired_given_error=cr.p_fired_given_error,
                    p_error_given_fired=cr.p_error_given_fired,
                    lift=cr.lift,
                    num_logical_errors=num_logical_errors,
                )
            )

        config_reports = []
        for ce in cpp_result.config_reports:
            cats = (
                ce.configuration.split(" + ")
                if ce.configuration != "(no noise)"
                else []
            )
            config_reports.append({
                "configuration": ce.configuration,
                "count": ce.count,
                "fraction": ce.fraction,
                "categories": cats,
            })

        return HotspotResult(
            num_shots=num_shots,
            num_logical_errors=num_logical_errors,
            ler=ler,
            category_reports=category_reports,
            config_reports=config_reports,
        )

    @staticmethod
    def print_report(result: HotspotResult) -> None:
        """Print a formatted hotspot report to stdout.

        Args:
            result: A :class:`HotspotResult` from :meth:`analyze`.
        """
        reports = result.category_reports
        configs = result.config_reports

        if not reports:
            print("No logical errors observed -- nothing to report.")
            return

        num_errors = result.num_logical_errors
        print(f"\n{'='*70}")
        print(
            f"  Hotspot Analysis  "
            f"({num_errors} logical errors, "
            f"LER = {result.ler:.6f}, "
            f"{result.num_shots} shots)"
        )
        print(f"{'='*70}")
        print(
            f"  {'Category':<25s} {'Count':>5s}  "
            f"{'P(fire)':>8s}  {'P(fire|err)':>11s}  "
            f"{'P(err|fire)':>11s}  {'Lift':>6s}"
        )
        print(f"  {'-'*25} {'-'*5}  {'-'*8}  {'-'*11}  {'-'*11}  {'-'*6}")
        for r in reports:
            print(
                f"  {r.category:<25s} {r.count:>5d}  "
                f"{r.p_fired:>8.4f}  {r.p_fired_given_error:>11.4f}  "
                f"{r.p_error_given_fired:>11.4f}  {r.lift:>6.2f}"
            )
        print(f"{'='*70}")

        if configs:
            print(
                f"\n  Multi-error configurations "
                f"(among {num_errors} logical errors):"
            )
            print(f"  {'-'*55}")
            print(
                f"  {'Configuration':<40s} {'Count':>6s}  {'Fraction':>8s}"
            )
            print(f"  {'-'*40} {'-'*6}  {'-'*8}")
            for c in configs:
                print(
                    f"  {c['configuration']:<40s} {c['count']:>6d}  "
                    f"{c['fraction']:>8.1%}"
                )
            print(f"  {'-'*55}\n")
