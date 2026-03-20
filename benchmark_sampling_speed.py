"""Benchmark: three sampling backends for surface codes.

Compares pure sampling throughput (samples/second):
1. Stim detector sampler
2. QEPG Monte Carlo (Poisson error injection)
3. QEPG fixed-weight (stratified sampling at weight = code distance)
"""

import time
import json
import numpy as np
import stim

# Force unbuffered output
def print(*args, **kwargs):
    kwargs.setdefault("flush", True)
    __builtins__["print"](*args, **kwargs) if isinstance(__builtins__, dict) else __builtins__.print(*args, **kwargs)

from scalerqec.Clifford.stimparser import rewrite_stim_code
from scalerqec.qepg import (
    compile_QEPG,
    return_samples_Monte_separate_obs_with_QEPG,
    return_samples_with_fixed_QEPG_numpy,
)

DISTANCES = [3, 5, 7, 9, 11]
ERROR_RATE = 0.001
WARMUP_SHOTS = 100
BENCH_SHOTS = 100_000

results = {}

for d in DISTANCES:
    print(f"\n{'='*60}")
    print(f"  Surface code d={d}  (n={d*d} data qubits, rounds={d})")
    print(f"{'='*60}")

    filepath = f"stimprograms/surface/surface{d}"
    with open(filepath, "r") as f:
        stim_str = f.read()
    stim_circuit = stim.Circuit(stim_str)

    num_detectors = stim_circuit.num_detectors
    num_qubits = stim_circuit.num_qubits
    print(f"  Qubits: {num_qubits}, Detectors: {num_detectors}")

    # ---- 1. Stim detector sampler ----
    sampler = stim_circuit.compile_detector_sampler()
    sampler.sample(WARMUP_SHOTS, separate_observables=True)

    t0 = time.perf_counter()
    sampler.sample(BENCH_SHOTS, separate_observables=True)
    stim_time = time.perf_counter() - t0
    stim_rate = BENCH_SHOTS / stim_time

    print(f"  Stim:            {stim_rate/1e6:.2f}M samples/s  ({stim_time:.3f} s)")

    # ---- 2. QEPG Monte Carlo ----
    qepg_stim_str = rewrite_stim_code(stim_str, keep_noise=True)
    qepg_graph = compile_QEPG(qepg_stim_str)
    return_samples_Monte_separate_obs_with_QEPG(qepg_graph, ERROR_RATE, WARMUP_SHOTS)

    t0 = time.perf_counter()
    return_samples_Monte_separate_obs_with_QEPG(qepg_graph, ERROR_RATE, BENCH_SHOTS)
    qepg_mc_time = time.perf_counter() - t0
    qepg_mc_rate = BENCH_SHOTS / qepg_mc_time

    print(f"  QEPG Monte Carlo:{qepg_mc_rate/1e6:.2f}M samples/s  ({qepg_mc_time:.3f} s)  [{qepg_mc_rate/stim_rate:.1f}x]")

    # ---- 3. QEPG fixed-weight (stratified) ----
    weight = d  # sample at weight = code distance
    return_samples_with_fixed_QEPG_numpy(qepg_graph, weight, WARMUP_SHOTS)

    t0 = time.perf_counter()
    return_samples_with_fixed_QEPG_numpy(qepg_graph, weight, BENCH_SHOTS)
    qepg_fw_time = time.perf_counter() - t0
    qepg_fw_rate = BENCH_SHOTS / qepg_fw_time

    print(f"  QEPG Fixed-w={weight:>2d}:  {qepg_fw_rate/1e6:.2f}M samples/s  ({qepg_fw_time:.3f} s)  [{qepg_fw_rate/stim_rate:.1f}x]")

    results[str(d)] = {
        "num_qubits": num_qubits,
        "num_detectors": num_detectors,
        "stim_rate": stim_rate,
        "stim_time_s": stim_time,
        "qepg_mc_rate": qepg_mc_rate,
        "qepg_mc_time_s": qepg_mc_time,
        "qepg_fw_rate": qepg_fw_rate,
        "qepg_fw_time_s": qepg_fw_time,
        "speedup_mc": qepg_mc_rate / stim_rate,
        "speedup_fw": qepg_fw_rate / stim_rate,
    }

with open("benchmark_results.json", "w") as f:
    json.dump(results, f, indent=2)

print(f"\n{'='*60}")
print("  Summary")
print(f"{'='*60}")
print(f"  {'d':>3s}  {'Stim (M/s)':>12s}  {'QEPG-MC (M/s)':>15s}  {'QEPG-FW (M/s)':>15s}  {'MC':>5s}  {'FW':>5s}")
print(f"  {'---':>3s}  {'----------':>12s}  {'-------------':>15s}  {'-------------':>15s}  {'--':>5s}  {'--':>5s}")
for d in DISTANCES:
    r = results[str(d)]
    print(
        f"  {d:>3d}  {r['stim_rate']/1e6:>12.2f}  {r['qepg_mc_rate']/1e6:>15.2f}  {r['qepg_fw_rate']/1e6:>15.2f}"
        f"  {r['speedup_mc']:>4.1f}x  {r['speedup_fw']:>4.1f}x"
    )

print(f"\nResults saved to benchmark_results.json")
