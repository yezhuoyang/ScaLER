"""
Benchmark: QEPG (SIMD-optimized) vs Stim sampling speed.

Compares:
  1. Fixed-weight sampling rate (QEPG) vs Stim detector sampler
  2. Monte Carlo sampling rate (QEPG) vs Stim detector sampler
  3. Multi-thread scaling (1, 2, 4, 8 threads)

Usage:
    python benchmark/benchmark_simd_vs_stim.py
"""
import os
import sys
import time
import numpy as np
import stim

# Force non-interactive backend before any matplotlib import
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from scalerqec.qepg import (
    compile_QEPG,
    return_samples_many_weights_separate_obs_with_QEPG,
    return_samples_Monte_separate_obs_with_QEPG,
    return_samples_numpy,
)
from scalerqec.Clifford.clifford import CliffordCircuit
from scalerqec.Clifford.stimparser import rewrite_stim_code

STIM_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "stimprograms", "surface")
NUM_SHOTS = 1_000_000
ERROR_RATE = 0.001


def load_circuit(distance):
    """Load surface code circuit from stimprograms directory."""
    filepath = os.path.join(STIM_DIR, f"surface{distance}")
    if not os.path.exists(filepath):
        return None, None, None, None

    with open(filepath, "r") as f:
        stim_str = f.read()

    circuit = CliffordCircuit(2)
    rewritten = rewrite_stim_code(stim_str)
    circuit.compile_from_stim_circuit_str(rewritten)
    stim_circuit = circuit.stimcircuit
    total_noise = circuit.totalnoise
    prog_str = str(stim_circuit)
    return prog_str, stim_circuit, total_noise, circuit


def benchmark_stim_detector_sampler(stim_circuit, shots):
    """Benchmark Stim's detector sampler (the fair comparison)."""
    sampler = stim_circuit.compile_detector_sampler()
    t0 = time.perf_counter()
    sampler.sample(shots=shots, separate_observables=True)
    elapsed = time.perf_counter() - t0
    return elapsed


def benchmark_stim_sampler(stim_circuit, shots):
    """Benchmark Stim's circuit sampler (raw measurement sampling)."""
    sampler = stim_circuit.compile_sampler()
    t0 = time.perf_counter()
    sampler.sample(shots=shots)
    elapsed = time.perf_counter() - t0
    return elapsed


def benchmark_qepg_fixed_weight(graph, weight, shots):
    """Benchmark QEPG fixed-weight sampling."""
    t0 = time.perf_counter()
    return_samples_many_weights_separate_obs_with_QEPG(graph, [weight], [shots])
    elapsed = time.perf_counter() - t0
    return elapsed


def benchmark_qepg_monte(graph, error_rate, shots):
    """Benchmark QEPG Monte Carlo sampling."""
    t0 = time.perf_counter()
    return_samples_Monte_separate_obs_with_QEPG(graph, error_rate, shots)
    elapsed = time.perf_counter() - t0
    return elapsed


def run_speed_comparison():
    """Compare sampling speed across surface code distances."""
    distances = [3, 5, 7, 9, 11, 13, 15]
    shots = NUM_SHOTS

    print("=" * 80)
    print(f"SAMPLING SPEED COMPARISON: QEPG vs Stim  ({shots:,} shots)")
    print("=" * 80)

    results = {
        'distance': [], 'noise_locs': [],
        'qepg_fixed_time': [], 'qepg_fixed_rate': [],
        'qepg_monte_time': [], 'qepg_monte_rate': [],
        'stim_det_time': [], 'stim_det_rate': [],
        'stim_raw_time': [], 'stim_raw_rate': [],
    }

    for d in distances:
        prog_str, stim_circuit, total_noise, circ = load_circuit(d)
        if prog_str is None:
            print(f"  [SKIP] surface{d} not found")
            continue

        weight = max(1, int(total_noise * ERROR_RATE))

        # Compile QEPG
        t0 = time.perf_counter()
        graph = compile_QEPG(prog_str)
        compile_time = time.perf_counter() - t0

        print(f"\n--- Surface d={d}  (noise_locs={total_noise}, weight={weight}, compile={compile_time:.3f}s) ---")

        # Warmup run
        return_samples_many_weights_separate_obs_with_QEPG(graph, [weight], [1000])

        # QEPG fixed-weight
        t_fw = benchmark_qepg_fixed_weight(graph, weight, shots)
        rate_fw = shots / t_fw

        # QEPG Monte Carlo
        t_mc = benchmark_qepg_monte(graph, ERROR_RATE, shots)
        rate_mc = shots / t_mc

        # Stim detector sampler (fair comparison)
        t_det = benchmark_stim_detector_sampler(stim_circuit, shots)
        rate_det = shots / t_det

        # Stim raw sampler
        t_raw = benchmark_stim_sampler(stim_circuit, shots)
        rate_raw = shots / t_raw

        results['distance'].append(d)
        results['noise_locs'].append(total_noise)
        results['qepg_fixed_time'].append(t_fw)
        results['qepg_fixed_rate'].append(rate_fw)
        results['qepg_monte_time'].append(t_mc)
        results['qepg_monte_rate'].append(rate_mc)
        results['stim_det_time'].append(t_det)
        results['stim_det_rate'].append(rate_det)
        results['stim_raw_time'].append(t_raw)
        results['stim_raw_rate'].append(rate_raw)

        print(f"  QEPG fixed-weight:  {t_fw:.4f}s  ({rate_fw/1e6:.2f}M samples/s)")
        print(f"  QEPG Monte Carlo:   {t_mc:.4f}s  ({rate_mc/1e6:.2f}M samples/s)")
        print(f"  Stim det-sampler:   {t_det:.4f}s  ({rate_det/1e6:.2f}M samples/s)")
        print(f"  Stim raw-sampler:   {t_raw:.4f}s  ({rate_raw/1e6:.2f}M samples/s)")
        print(f"  Speedup (QEPG-FW vs Stim-det): {t_det/t_fw:.1f}x")
        print(f"  Speedup (QEPG-MC vs Stim-det): {t_det/t_mc:.1f}x")

    return results


def run_thread_scaling():
    """Test multi-thread scaling using subprocesses.

    OMP_NUM_THREADS must be set before the OpenMP runtime initializes,
    so each thread count is tested in a separate subprocess.
    """
    import subprocess, json

    distances_to_test = [7, 11]
    thread_counts = [1, 2, 4, 8, 16, 32]
    shots = NUM_SHOTS

    print("\n" + "=" * 80)
    print(f"MULTI-THREAD SCALING TEST ({shots:,} shots)")
    print("=" * 80)

    all_scaling = {}

    for d in distances_to_test:
        print(f"\n  --- Surface d={d} ---")

        scaling_results = {
            'threads': [], 'fixed_time': [], 'fixed_rate': [],
            'monte_time': [], 'monte_rate': [],
        }

        # Reduce MC shots for large codes to keep benchmark reasonable
        mc_shots = shots if d <= 7 else shots // 10

        for n_threads in thread_counts:
            env = os.environ.copy()
            env['OMP_NUM_THREADS'] = str(n_threads)

            script = f"""
import os, time, json
from scalerqec.qepg import compile_QEPG, return_samples_many_weights_separate_obs_with_QEPG, return_samples_Monte_separate_obs_with_QEPG
from scalerqec.Clifford.clifford import CliffordCircuit
from scalerqec.Clifford.stimparser import rewrite_stim_code

with open('stimprograms/surface/surface{d}') as f: s=f.read()
c=CliffordCircuit(2); r=rewrite_stim_code(s); c.compile_from_stim_circuit_str(r)
prog=str(c.stimcircuit); noise=c.totalnoise; w=max(1,int(noise*{ERROR_RATE}))

g=compile_QEPG(prog)
return_samples_many_weights_separate_obs_with_QEPG(g,[w],[1000])

t0=time.perf_counter()
return_samples_many_weights_separate_obs_with_QEPG(g,[w],[{shots}])
t1=time.perf_counter()

t2=time.perf_counter()
return_samples_Monte_separate_obs_with_QEPG(g,{ERROR_RATE},{mc_shots})
t3=time.perf_counter()

print(json.dumps({{'fw': t1-t0, 'mc': t3-t2}}))
"""
            result = subprocess.run(
                [sys.executable, '-c', script],
                capture_output=True, text=True, env=env,
                cwd=os.path.dirname(os.path.dirname(__file__)),
                timeout=600,
            )
            if result.returncode != 0:
                print(f"    {n_threads} threads: FAILED - {result.stderr[:200]}")
                continue

            data = json.loads(result.stdout.strip().split('\n')[-1])
            t_fw = data['fw']
            t_mc = data['mc']
            rate_fw = shots / t_fw
            rate_mc = mc_shots / t_mc

            scaling_results['threads'].append(n_threads)
            scaling_results['fixed_time'].append(t_fw)
            scaling_results['fixed_rate'].append(rate_fw)
            scaling_results['monte_time'].append(t_mc)
            scaling_results['monte_rate'].append(rate_mc)

            print(f"    {n_threads:2d} thread(s):  FW={t_fw:.4f}s ({rate_fw/1e6:.2f}M/s)  MC={t_mc:.4f}s ({rate_mc/1e6:.2f}M/s)")

        # Print scaling efficiency
        if len(scaling_results['threads']) > 1:
            base_fw = scaling_results['fixed_rate'][0]
            base_mc = scaling_results['monte_rate'][0]
            print(f"\n    Thread scaling efficiency (d={d}):")
            for i, nt in enumerate(scaling_results['threads']):
                speedup_fw = scaling_results['fixed_rate'][i] / base_fw
                speedup_mc = scaling_results['monte_rate'][i] / base_mc
                eff_fw = speedup_fw / nt * 100
                eff_mc = speedup_mc / nt * 100
                print(f"      {nt:2d} threads: FW {speedup_fw:.2f}x (eff={eff_fw:.0f}%)  MC {speedup_mc:.2f}x (eff={eff_mc:.0f}%)")

        all_scaling[d] = scaling_results

    return all_scaling


def plot_results(speed_results, all_scaling):
    """Generate comparison plots."""
    n_scaling = len(all_scaling) if all_scaling else 0
    n_cols = 2 + n_scaling
    fig, axes = plt.subplots(1, n_cols, figsize=(6 * n_cols, 5))

    # --- Plot 1: Sampling rate vs distance ---
    ax = axes[0]
    ds = speed_results['distance']
    ax.semilogy(ds, [r/1e6 for r in speed_results['qepg_fixed_rate']], 'o-', label='QEPG Fixed-Weight', linewidth=2)
    ax.semilogy(ds, [r/1e6 for r in speed_results['qepg_monte_rate']], 's-', label='QEPG Monte Carlo', linewidth=2)
    ax.semilogy(ds, [r/1e6 for r in speed_results['stim_det_rate']], '^--', label='Stim Det-Sampler', linewidth=2)
    ax.semilogy(ds, [r/1e6 for r in speed_results['stim_raw_rate']], 'v--', label='Stim Raw-Sampler', linewidth=2)
    ax.set_xlabel('Code Distance')
    ax.set_ylabel('Sampling Rate (M samples/s)')
    ax.set_title('Sampling Speed vs Code Distance')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- Plot 2: Speedup ratio vs distance ---
    ax = axes[1]
    speedup_fw = [speed_results['stim_det_time'][i] / speed_results['qepg_fixed_time'][i] for i in range(len(ds))]
    speedup_mc = [speed_results['stim_det_time'][i] / speed_results['qepg_monte_time'][i] for i in range(len(ds))]
    ax.plot(ds, speedup_fw, 'o-', label='QEPG-FW vs Stim-Det', linewidth=2)
    ax.plot(ds, speedup_mc, 's-', label='QEPG-MC vs Stim-Det', linewidth=2)
    ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Code Distance')
    ax.set_ylabel('Speedup (x)')
    ax.set_title('QEPG Speedup over Stim')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- Plot 3+: Thread scaling for each distance ---
    if all_scaling:
        for idx, (d, scaling_results) in enumerate(all_scaling.items()):
            ax = axes[2 + idx]
            threads = scaling_results['threads']
            base_fw = scaling_results['fixed_rate'][0]
            base_mc = scaling_results['monte_rate'][0]
            ax.plot(threads, [r / base_fw for r in scaling_results['fixed_rate']], 'o-', label='Fixed-Weight', linewidth=2)
            ax.plot(threads, [r / base_mc for r in scaling_results['monte_rate']], 's-', label='Monte Carlo', linewidth=2)
            max_t = max(threads)
            ax.plot([1, max_t], [1, max_t], 'k--', alpha=0.3, label='Ideal linear')
            ax.set_xlabel('Number of Threads')
            ax.set_ylabel('Speedup (x)')
            ax.set_title(f'Thread Scaling (d={d})')
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)

    plt.tight_layout()
    outpath = os.path.join(os.path.dirname(__file__), 'benchmark_simd_vs_stim.png')
    plt.savefig(outpath, dpi=150)
    print(f"\nPlot saved to: {outpath}")


if __name__ == "__main__":
    speed_results = run_speed_comparison()
    scaling_results = run_thread_scaling()
    plot_results(speed_results, scaling_results)
