# SIMD & Memory Layout Optimization for QEPG Sampling

**Date:** 2026-03-18
**Branch:** `Type`
**Author:** Claude (AI pair-programming assistant)

## Objective

Accelerate QEPG sampling speed through SIMD optimization, cache-friendly memory layout, allocation elimination, and faster RNG. Goal: beat Stim's single-thread sampling speed without multi-core acceleration.

## Changes Summary

### Phase 1: Eliminate Per-Sample Allocations

**Files modified:** `QEPG/src/sampler.hpp`, `QEPG/src/sampler.cpp`

- **Bitmap collision detection**: Replaced `std::unordered_set<size_t>` with `std::vector<uint8_t>` bitmap for O(1) collision checking in Floyd sampling. Clears only used positions (O(weight) instead of O(N)).
- **Pre-allocated scratch buffers**: Added `std::vector<singlePauli> scratch_sample_` member reused across samples, eliminating per-sample heap allocation.
- **Branchless Pauli dispatch**: Replaced 3-way if/else chain with direct index arithmetic:
  ```cpp
  size_t row_idx = sample[s].qindex + (sample[s].type - 1) * n_noise;
  result ^= dm[row_idx];
  ```

### Phase 2: Contiguous Matrix Memory Layout

**Files created:** `QEPG/src/flat_bit_table.hpp`
**Files modified:** `QEPG/src/QEPG.hpp`, `QEPG/src/QEPG.cpp`

- **FlatBitTable**: Contiguous 2D bit matrix with 64-byte cache-line aligned rows in a single allocation. Eliminates pointer chasing from `std::vector<Row>` (each Row was a separate heap allocation).
- After `backward_graph_construction()`, the parity propagation matrix transpose is copied into a FlatBitTable for sampling use.
- Added move semantics to `QEPG` class (FlatBitTable is non-copyable).

### Phase 3: SIMD-Accelerated XOR

**Files created:** `QEPG/src/simd_util.hpp`
**Files modified:** `QEPG/src/dynamic_bitset.hpp`

- **Compile-time SIMD dispatch** for XOR operations:
  - AVX2: 256-bit (4 x uint64_t per instruction)
  - SSE2: 128-bit (2 x uint64_t per instruction) — MSVC x64 always has this
  - ARM NEON: 128-bit for Apple Silicon
  - Scalar fallback
- **Cross-platform aligned memory allocation**: `_aligned_malloc` (MSVC), `posix_memalign` (POSIX), `std::aligned_alloc` (C++17)
- **DynamicBitset raw constructor**: `DynamicBitset(num_bits, const block_type* src, n_words)` for zero-copy Row creation from SIMD result buffers.

### Phase 5: Faster RNG (xoshiro256++)

**Files modified:** `QEPG/src/sampler.hpp`, `QEPG/src/sampler.cpp`

- Replaced `std::mt19937` (2.5KB state, cache pollution) with `Xoshiro256pp` (32 bytes, fits in registers).
- ~3x faster per call, passes BigCrush.
- Uses SplitMix64 for seeding from a single uint64_t.
- Monte Carlo path uses threshold comparison on upper 32 bits instead of `std::bernoulli_distribution`.

### OpenMP Multi-Threading

- Each OpenMP thread creates its own `sampler` object (with bitmap + scratch) inside `#pragma omp parallel` — avoids `thread_local` issues on MSVC.
- Each thread gets its own `Xoshiro256pp` RNG seeded with `global_seed ^ thread_id_hash`.

## Benchmark Results

### QEPG vs Stim Speed Comparison (1M shots, all cores)

| Distance | Noise Locs | QEPG Fixed-Weight | Stim Det-Sampler | Speedup |
|----------|-----------|-------------------|-----------------|---------|
| d=3      | 585       | 28.28M/s          | 2.91M/s         | **9.7x**  |
| d=5      | 3,145     | 13.21M/s          | 0.38M/s         | **34.5x** |
| d=7      | 9,121     | 7.34M/s           | 0.16M/s         | **46.9x** |
| d=9      | 19,953    | 3.45M/s           | 0.06M/s         | **54.6x** |
| d=11     | 37,081    | 1.92M/s           | 0.03M/s         | **57.1x** |
| d=13     | 61,945    | 0.76M/s           | 0.02M/s         | **38.0x** |
| d=15     | 95,985    | 0.26M/s           | 0.01M/s         | **20.6x** |

QEPG fixed-weight sampling is **20-57x faster** than Stim's detector sampler.

### Multi-Thread Scaling (Surface d=7, 1M shots)

| Threads | FW Rate   | FW Speedup | FW Efficiency | MC Speedup | MC Efficiency |
|---------|----------|------------|---------------|------------|---------------|
| 1       | 2.14M/s  | 1.00x      | 100%          | 1.00x      | 100%          |
| 2       | 3.51M/s  | 1.64x      | 82%           | 1.95x      | 98%           |
| 4       | 5.18M/s  | 2.43x      | 61%           | 3.58x      | 90%           |
| 8       | 6.17M/s  | 2.89x      | 36%           | 6.71x      | 84%           |
| 16      | 8.14M/s  | 3.81x      | 24%           | 12.45x     | 78%           |
| 32      | 7.76M/s  | 3.63x      | 11%           | 18.60x     | 58%           |

### Multi-Thread Scaling (Surface d=11, 1M shots)

| Threads | FW Rate   | FW Speedup | FW Efficiency | MC Speedup | MC Efficiency |
|---------|----------|------------|---------------|------------|---------------|
| 1       | 0.26M/s  | 1.00x      | 100%          | 1.00x      | 100%          |
| 2       | 0.39M/s  | 1.48x      | 74%           | 2.04x      | 102%          |
| 4       | 0.79M/s  | 3.01x      | 75%           | 3.73x      | 93%           |
| 8       | 1.05M/s  | 4.00x      | 50%           | 6.73x      | 84%           |
| 16      | 1.46M/s  | 5.57x      | 35%           | 13.27x     | 83%           |
| 32      | 1.53M/s  | 5.84x      | 18%           | 18.06x     | 56%           |

**Key observations:**
- Monte Carlo scales nearly linearly up to 16 threads (~83% efficiency)
- Fixed-weight scales well up to 8 threads, then diminishes (work per sample is very light for small weights)
- Larger codes (d=11) show better FW scaling than smaller codes (d=7) due to more work per sample

### Benchmark Plot

![Benchmark Results](benchmark_simd_vs_stim.png)

## Test Results

All 392 tests pass (2 skipped, 1 xfailed):
```
=========== 392 passed, 2 skipped, 1 xfailed, 2 warnings in 26.69s ============
```

## Files Created/Modified

| File | Action | Description |
|------|--------|-------------|
| `QEPG/src/simd_util.hpp` | **NEW** | SIMD XOR dispatch + aligned memory allocation |
| `QEPG/src/flat_bit_table.hpp` | **NEW** | Contiguous 2D bit matrix with cache-line aligned rows |
| `QEPG/src/sampler.hpp` | Modified | Xoshiro256pp RNG, bitmap collision, scratch buffers, branchless dispatch |
| `QEPG/src/sampler.cpp` | Modified | All sampling methods updated for new optimizations |
| `QEPG/src/QEPG.hpp` | Modified | FlatBitTable member, move semantics |
| `QEPG/src/QEPG.cpp` | Modified | Build FlatBitTable after graph construction |
| `QEPG/src/dynamic_bitset.hpp` | Modified | Raw-words constructor for zero-copy Row creation |
| `tests/conftest.py` | Modified | Non-interactive matplotlib backend (Agg) |
| `benchmark/benchmark_simd_vs_stim.py` | **NEW** | Comprehensive benchmark script |

## Known Issues / Future Work

1. **Direct-to-numpy path**: `generate_many_output_samples_to_numpy` and `_Monte_to_numpy` methods exist in sampler.cpp but crash with access violation. Currently using intermediate Row approach (works correctly).
2. **Phase 4 (Batch Monte Carlo)**: Not implemented. Would improve Monte Carlo single-thread performance by batching Bernoulli sampling and eliminating per-location RNG iteration.
3. **Monte Carlo bottleneck**: MC iterates all noise locations per sample O(N), making it slower than fixed-weight for large codes. Phase 4 batch approach would fix this.
