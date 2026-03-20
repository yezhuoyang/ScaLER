"""Verify that all three compilation approaches produce matching LER.

For each (distance, noise_model) combination, this test:
1. Approach 3 (IR): SurfaceCode + NoiseModel → construct_circuit()
2. Approach 2 (noiseless + inject): compile_from_stim_circuit_str → inject_noise
3. Approach 1 (noisy parse): compile_from_noisy_stim_circuit_str

Then samples LER via our Monte Carlo interface and verifies agreement.
"""

import os
import tempfile

import numpy as np
import pytest
import stim
import pymatching

from scalerqec.QEC.surface import SurfaceCode
from scalerqec.QEC.noisemodel import (
    NoiseModel,
    SD6NoiseModel,
    SI1000NoiseModel,
    SIDNoiseModel,
)
from scalerqec.Clifford.clifford import CliffordCircuit
from scalerqec.Clifford.stimparser import rewrite_stim_code
from scalerqec.Monte.monteLER import count_logical_errors


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _build_noiseless_surface_code(d: int, rounds: int | None = None) -> SurfaceCode:
    """Build a noiseless SurfaceCode and compile to stim circuit."""
    code = SurfaceCode(distance=d, rounds=rounds or d)
    code.scheme = "Standard"
    code.construct_circuit()
    return code


def _ler_from_stim_circuit(circuit: stim.Circuit, shots: int) -> float:
    """Estimate LER using stim sampler + pymatching."""
    num_errors = count_logical_errors(circuit, shots)
    return num_errors / shots


def _build_via_approach3(d: int, noise_model: NoiseModel, rounds: int | None = None) -> stim.Circuit:
    """Approach 3: StabCode IR compilation with NoiseModel."""
    code = SurfaceCode(distance=d, rounds=rounds or d)
    code.scheme = "Standard"
    code.noisemodel = noise_model
    code.construct_circuit()
    return code.stimcirc


def _build_via_approach2(d: int, noise_model: NoiseModel, rounds: int | None = None) -> stim.Circuit:
    """Approach 2: Noiseless compile + inject noise."""
    code = _build_noiseless_surface_code(d, rounds)
    noiseless_stim = code.stimcirc
    noisy_stim = noise_model.inject_noise(noiseless_stim)
    return noisy_stim


def _build_via_approach1(d: int, noise_model: NoiseModel, rounds: int | None = None) -> stim.Circuit:
    """Approach 1: Build noisy circuit, then parse with noisy parser."""
    # First get the noisy stim string (same as approach 2)
    noisy_stim = _build_via_approach2(d, noise_model, rounds)
    noisy_str = str(noisy_stim)

    # Now parse via the noisy interface (builds QEPG gate list)
    circ = CliffordCircuit(2)
    circ.compile_from_noisy_stim_circuit_str(noisy_str)

    # The stim circuit stored in CliffordCircuit is the parsed noisy one
    return circ.stimcircuit


# ---------------------------------------------------------------------------
# DEM equivalence check — strongest verification that circuits match
# ---------------------------------------------------------------------------

def _dem_instructions(circuit: stim.Circuit) -> list[str]:
    """Extract sorted DEM instruction strings for comparison."""
    dem = circuit.detector_error_model(decompose_errors=True)
    return sorted(str(dem).strip().splitlines())


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

NOISE_MODELS = {
    "SID": lambda p: SIDNoiseModel(p),
    "SD6": lambda p: SD6NoiseModel(p),
    "SI1000": lambda p: SI1000NoiseModel(p),
}

ERROR_RATE = 0.005
SHOTS = 100_000


class TestThreeApproachesDEMEquivalence:
    """Verify that all three approaches produce identical DEMs.

    If the DEMs match, the LER must match (same error model).
    """

    @pytest.mark.parametrize("d", [3, 5])
    @pytest.mark.parametrize("noise_name", ["SID", "SD6", "SI1000"])
    def test_dem_equivalence(self, d, noise_name):
        """All three approaches produce the same DEM."""
        noise_model = NOISE_MODELS[noise_name](ERROR_RATE)

        circ1 = _build_via_approach1(d, noise_model)
        circ2 = _build_via_approach2(d, noise_model)
        circ3 = _build_via_approach3(d, noise_model)

        dem1 = _dem_instructions(circ1)
        dem2 = _dem_instructions(circ2)
        dem3 = _dem_instructions(circ3)

        assert dem1 == dem2, f"Approach 1 vs 2 DEM mismatch for d={d} {noise_name}"
        assert dem2 == dem3, f"Approach 2 vs 3 DEM mismatch for d={d} {noise_name}"

    @pytest.mark.parametrize("noise_name", ["SID", "SD6", "SI1000"])
    def test_dem_equivalence_d7(self, noise_name):
        """DEM equivalence for d=7 surface code."""
        noise_model = NOISE_MODELS[noise_name](ERROR_RATE)

        circ1 = _build_via_approach1(7, noise_model)
        circ2 = _build_via_approach2(7, noise_model)
        circ3 = _build_via_approach3(7, noise_model)

        dem1 = _dem_instructions(circ1)
        dem2 = _dem_instructions(circ2)
        dem3 = _dem_instructions(circ3)

        assert dem1 == dem2, f"Approach 1 vs 2 DEM mismatch for d=7 {noise_name}"
        assert dem2 == dem3, f"Approach 2 vs 3 DEM mismatch for d=7 {noise_name}"


class TestThreeApproachesLER:
    """Verify that LER sampled from each approach matches.

    Uses our Monte Carlo interface (count_logical_errors) to sample
    LER and checks they agree within statistical tolerance.
    """

    @pytest.mark.parametrize("d", [3, 5])
    @pytest.mark.parametrize("noise_name", ["SID", "SD6", "SI1000"])
    def test_ler_match(self, d, noise_name):
        """LER from all three approaches must agree within tolerance."""
        noise_model = NOISE_MODELS[noise_name](ERROR_RATE)

        circ1 = _build_via_approach1(d, noise_model)
        circ2 = _build_via_approach2(d, noise_model)
        circ3 = _build_via_approach3(d, noise_model)

        ler1 = _ler_from_stim_circuit(circ1, SHOTS)
        ler2 = _ler_from_stim_circuit(circ2, SHOTS)
        ler3 = _ler_from_stim_circuit(circ3, SHOTS)

        print(f"\n  d={d} {noise_name}: "
              f"A1={ler1:.4f}, A2={ler2:.4f}, A3={ler3:.4f}")

        # All should be > 0 (circuit actually has errors)
        assert ler1 > 0, f"Approach 1 LER is 0 for d={d} {noise_name}"
        assert ler2 > 0, f"Approach 2 LER is 0 for d={d} {noise_name}"
        assert ler3 > 0, f"Approach 3 LER is 0 for d={d} {noise_name}"

        # With 100K shots, statistical error ~ sqrt(ler/shots)
        # Allow 50% relative tolerance to account for MC variance
        mean_ler = (ler1 + ler2 + ler3) / 3
        for label, ler in [("A1", ler1), ("A2", ler2), ("A3", ler3)]:
            ratio = ler / mean_ler
            assert 0.5 < ratio < 2.0, (
                f"{label} LER={ler:.5f} vs mean={mean_ler:.5f} "
                f"(ratio={ratio:.2f}) for d={d} {noise_name}"
            )


class TestMonteLERCalcInterface:
    """Verify using MonteLERcalc.calculate_LER_from_stim_circuit."""

    @pytest.mark.parametrize("noise_name", ["SID", "SD6", "SI1000"])
    def test_monte_ler_calc_d3(self, noise_name):
        """MonteLERcalc gives consistent LER across approaches for d=3."""
        from scalerqec.Monte.monteLER import MonteLERcalc

        noise_model = NOISE_MODELS[noise_name](ERROR_RATE)

        circ2 = _build_via_approach2(3, noise_model)
        circ3 = _build_via_approach3(3, noise_model)

        calc = MonteLERcalc(time_budget=60, samplebudget=SHOTS, MIN_NUM_LE_EVENT=20)

        ler2 = calc.calculate_LER_from_stim_circuit(str(circ2), samplebudget=SHOTS)
        ler3 = calc.calculate_LER_from_stim_circuit(str(circ3), samplebudget=SHOTS)

        print(f"\n  d=3 {noise_name} MonteLERcalc: A2={ler2:.4f}, A3={ler3:.4f}")

        assert ler2 > 0
        assert ler3 > 0
        mean_ler = (ler2 + ler3) / 2
        for label, ler in [("A2", ler2), ("A3", ler3)]:
            ratio = ler / mean_ler
            assert 0.5 < ratio < 2.0, (
                f"{label} LER={ler:.5f} vs mean={mean_ler:.5f}"
            )


class TestMonteLERFromFile:
    """Verify MonteLERcalc.calculate_LER_from_file (approach 2 with SID)."""

    @pytest.mark.parametrize("d", [3, 5])
    def test_from_file_matches_direct(self, d):
        """calculate_LER_from_file matches direct stim sampling."""
        from scalerqec.Monte.monteLER import MonteLERcalc

        # Build noiseless circuit and save to temp file
        code = _build_noiseless_surface_code(d)
        noiseless_str = rewrite_stim_code(str(code.stimcirc))

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".stim", delete=False
        ) as f:
            f.write(noiseless_str)
            tmppath = f.name

        try:
            # MonteLERcalc from file (uses SID internally)
            calc = MonteLERcalc(
                time_budget=60, samplebudget=SHOTS, MIN_NUM_LE_EVENT=20
            )
            ler_file = calc.calculate_LER_from_file(
                samplebudget=SHOTS, filepath=tmppath, pvalue=ERROR_RATE
            )

            # Direct approach: inject SID noise and sample
            noisy_circ = SIDNoiseModel(ERROR_RATE).inject_noise(code.stimcirc)
            ler_direct = _ler_from_stim_circuit(noisy_circ, SHOTS)

            print(f"\n  d={d} SID from_file={ler_file:.4f}, direct={ler_direct:.4f}")

            assert ler_file > 0
            assert ler_direct > 0
            mean_ler = (ler_file + ler_direct) / 2
            for label, ler in [("file", ler_file), ("direct", ler_direct)]:
                ratio = ler / mean_ler
                assert 0.5 < ratio < 2.0, (
                    f"{label} LER={ler:.5f} vs mean={mean_ler:.5f}"
                )
        finally:
            os.unlink(tmppath)


# ---------------------------------------------------------------------------
# Standalone runner — prints a comparison table
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 80)
    print("Three-Approach LER Verification")
    print(f"Error rate: {ERROR_RATE}, Shots: {SHOTS}")
    print("=" * 80)

    for d in [3, 5, 7]:
        print(f"\n{'-' * 70}")
        print(f"  Surface Code d={d}  (n={d*d}, rounds={d})")
        print(f"{'-' * 70}")
        print(f"  {'Noise':<10} {'Approach 1':>12} {'Approach 2':>12} {'Approach 3':>12}  {'Match?':>8}")

        for noise_name, noise_fn in NOISE_MODELS.items():
            noise_model = noise_fn(ERROR_RATE)

            circ1 = _build_via_approach1(d, noise_model)
            circ2 = _build_via_approach2(d, noise_model)
            circ3 = _build_via_approach3(d, noise_model)

            # DEM check
            dem1 = _dem_instructions(circ1)
            dem2 = _dem_instructions(circ2)
            dem3 = _dem_instructions(circ3)
            dem_ok = (dem1 == dem2 == dem3)

            # LER sampling
            shots = SHOTS if d <= 5 else SHOTS * 2
            ler1 = _ler_from_stim_circuit(circ1, shots)
            ler2 = _ler_from_stim_circuit(circ2, shots)
            ler3 = _ler_from_stim_circuit(circ3, shots)

            mean = (ler1 + ler2 + ler3) / 3
            if mean > 0:
                ratios = [abs(l / mean - 1) for l in [ler1, ler2, ler3]]
                ler_ok = all(r < 0.5 for r in ratios)
            else:
                ler_ok = (ler1 == ler2 == ler3 == 0)

            status = "OK" if (dem_ok and ler_ok) else "FAIL"
            dem_tag = "DEM=OK" if dem_ok else "DEM=FAIL"

            print(f"  {noise_name:<10} {ler1:>12.6f} {ler2:>12.6f} {ler3:>12.6f}  {status:>5} ({dem_tag})")

    print(f"\n{'=' * 80}")
    print("Done.")
