"""Tests for hotspot analysis with C++ backend."""

import numpy as np
import pytest

from scalerqec.Clifford.clifford import CliffordCircuit
from scalerqec.Clifford.stimparser import rewrite_stim_code
from scalerqec.Analysis.hotspot import HotspotAnalyzer, HotspotResult, CategoryReport
from scalerqec.qepg import compile_QEPG
import pymatching
import stim


def _build_repetition_code():
    """Build a noisy repetition code d=3 and return (qepg_graph, label_map, decoder, noise_probs)."""
    from scalerqec.QEC import StabCode
    from scalerqec.QEC.noisemodel import NoiseModel
    from scalerqec import qepg as qepg_cpp

    qeccirc = StabCode(n=3, k=1, d=3)
    qeccirc.add_stab("ZZI")
    qeccirc.add_stab("IZZ")
    qeccirc.set_logical_Z(0, "ZZZ")
    qeccirc.scheme = "Standard"
    qeccirc.rounds = 3

    noise_model = NoiseModel(0.01)  # higher rate for test visibility
    qeccirc.construct_circuit()
    stim_circuit = noise_model.inject_noise(qeccirc.stimcirc)

    # Build QEPG graph
    stim_str = str(stim_circuit)
    qepg_graph = compile_QEPG(rewrite_stim_code(stim_str, keep_noise=True))

    # Build CliffordCircuit for auto_label
    cc = CliffordCircuit(0)
    cc.compile_from_noisy_stim_circuit_str(stim_str)
    cpp_cc = qepg_cpp.CliffordCircuit()
    cpp_cc.compile_from_rewrited_stim_string(rewrite_stim_code(stim_str, keep_noise=True))

    # Build label map
    label_map = qepg_cpp.auto_label(cpp_cc, 3)  # 3 data qubits

    # Decoder
    dem = stim_circuit.detector_error_model(decompose_errors=False)
    decoder = pymatching.Matching.from_detector_error_model(dem)

    # Noise probs
    N = cc.totalnoise
    p = 0.01
    noise_probs = np.full((N, 3), [p / 3, p / 3, p / 3])

    return qepg_graph, label_map, decoder, noise_probs


class TestHotspotAnalyzer:
    def test_analyze_returns_result(self):
        qepg_graph, label_map, decoder, noise_probs = _build_repetition_code()

        analyzer = HotspotAnalyzer(qepg_graph, label_map, decoder=decoder)
        result = analyzer.analyze(noise_probs, num_shots=50000)

        assert isinstance(result, HotspotResult)
        assert result.num_shots == 50000
        assert result.num_logical_errors >= 0
        assert 0.0 <= result.ler <= 1.0
        assert isinstance(result.category_reports, list)

    def test_category_reports_structure(self):
        qepg_graph, label_map, decoder, noise_probs = _build_repetition_code()

        analyzer = HotspotAnalyzer(qepg_graph, label_map, decoder=decoder)
        result = analyzer.analyze(noise_probs, num_shots=100000)

        if result.num_logical_errors > 0:
            assert len(result.category_reports) > 0
            for r in result.category_reports:
                assert isinstance(r, CategoryReport)
                assert 0.0 <= r.p_fired <= 1.0
                assert 0.0 <= r.p_fired_given_error <= 1.0
                assert r.lift >= 0.0

    def test_print_report_no_crash(self, capsys):
        qepg_graph, label_map, decoder, noise_probs = _build_repetition_code()

        analyzer = HotspotAnalyzer(qepg_graph, label_map, decoder=decoder)
        result = analyzer.analyze(noise_probs, num_shots=50000)
        HotspotAnalyzer.print_report(result)

        captured = capsys.readouterr()
        assert "Hotspot Analysis" in captured.out or "No logical errors" in captured.out

    def test_zero_noise_no_errors(self):
        qepg_graph, label_map, decoder, noise_probs = _build_repetition_code()

        # Zero noise
        noise_probs_zero = np.zeros_like(noise_probs)

        analyzer = HotspotAnalyzer(qepg_graph, label_map, decoder=decoder)
        result = analyzer.analyze(noise_probs_zero, num_shots=1000)

        assert result.num_logical_errors == 0
        assert result.ler == 0.0
