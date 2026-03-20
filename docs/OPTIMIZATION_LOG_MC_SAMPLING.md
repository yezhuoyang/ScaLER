# Optimization Log: QEPG Monte Carlo Sampling

**Date:** 2026-03-19
**Author:** Claude (with Zhuoyang Ye)
**Result:** QEPG-MC went from **14-23x slower** than Stim to **3-9x faster** (single-threaded)

---

## 1. Problem Statement

We want to generate Monte Carlo detector samples for a QEC circuit. Given:
- N noise locations in the circuit
- Per-location error probability p (uniform) or per-location (px, py, pz) (non-uniform)
- For each shot: determine which noise locations fire, compute detector/observable outcomes via GF(2) XOR

**Goal:** Beat Stim's detector sampler throughput.

---

## 2. The Old Algorithm (Dense Bernoulli Scan)

```
For each shot:
    For each noise location i = 0, 1, ..., N-1:     // O(N) loop
        Draw uniform random r
        If r < p:
            Choose Pauli type t in {X, Y, Z}
            XOR row[i + (t-1)*N] into result_buf    // O(words_per_row) SIMD
    Unpack result_buf into NumPy arrays
```

**Cost per shot:** O(N) random number generations + O(k * W) XOR operations

where:
- N = total noise locations (e.g., 585 for d=3, 2745 for d=11)
- k = number of locations that actually fire (expected value E[k] = Np)
- W = words per detector row (e.g., 2 for 72 detectors)

**The bottleneck is the O(N) Bernoulli scan.** At p = 0.001, N = 585 (d=3):
- Expected fires: k = Np = 0.585 (less than 1 error per shot on average!)
- But we still iterate over all 585 locations, generating 585 random numbers per shot
- The XOR propagation is negligible (0-1 XOR operations on average)

This means **99.9% of the per-shot work produces no output** — we're paying O(N) to find ~O(1) active locations.

---

## 3. The New Algorithm (Sparse Poisson + Floyd)

### Key Insight

Instead of asking "does each location fire?" (N questions), ask:
1. "How many locations fire?" (1 question) → sample k
2. "Which k locations fire?" (k questions) → pick k distinct positions

### Step 1: Sample the error count k

The number of errors in a Bernoulli(N, p) process follows a Binomial(N, p) distribution.

For small p, the Poisson approximation is excellent:

$$\text{Binomial}(N, p) \approx \text{Poisson}(\lambda), \quad \lambda = Np$$

**Error bound:** The total variation distance between Binomial(N,p) and Poisson(Np) satisfies:

$$d_{TV}(\text{Bin}(N,p), \text{Poi}(Np)) \leq \frac{Np^2}{1-p} \approx Np^2$$

For p = 0.001, N = 585: error ≤ 585 × 10^{-6} ≈ 0.0006. Negligible.

**Poisson sampling (Knuth's algorithm for small λ):**

```
k = 0
product = 1.0
L = exp(-λ)
do:
    k += 1
    product *= uniform_random()
while product > L
return k - 1
```

Expected iterations: λ + 1. For λ = 0.585, this is ~1.6 iterations. For λ = 2.7 (d=11), ~3.7 iterations.

Compare to the old O(N) = O(585) to O(2745) iterations!

### Step 2: Pick k distinct locations (Floyd's Insertion Algorithm)

Given k, we need k distinct indices from {0, 1, ..., N-1}. Floyd's algorithm does this in O(k) expected time using a collision bitmap:

```
For each of k positions:
    Repeat:
        pos = uniform_random(0, N-1)
    Until pos is not already chosen (bitmap check)
    Mark pos as chosen
    Assign random Pauli type {X, Y, Z}
```

Expected random numbers per position: N/(N-j) for the j-th position ≈ 1 when k << N.

Total expected random numbers: k + O(k²/N) ≈ k for k << N.

### Step 3: XOR propagation (unchanged)

For each of the k active locations, XOR their detector row into the result buffer. Cost: O(k * W) with SIMD acceleration.

### Total cost comparison

| | Old (Dense) | New (Sparse) |
|--|------------|-------------|
| Random numbers per shot | N | λ + 1 + 2k ≈ 3λ + 1 |
| XOR operations per shot | k * W | k * W |
| **Total per shot** | **O(N)** | **O(k) = O(Np)** |

**Speedup factor:** N / (Np) = 1/p

At p = 0.001: theoretical speedup ≈ 1000x on the sampling step alone!

In practice the speedup is smaller because:
1. XOR propagation and NumPy unpacking have fixed overhead
2. Poisson sampling has higher per-iteration cost than a simple threshold compare
3. Floyd's bitmap operations have cache effects
4. The per-shot overhead (zero result buffer, unpack to NumPy) is O(W + n_det), not O(1)

---

## 4. The Second Optimization: Fused NumPy Output

### Old data path

```
Python → C++:
    1. Allocate vector<QEPG::Row> samplecontainer (N_shots × Row objects)
       Each Row = DynamicBitset = vector<uint64_t> (heap allocation per Row)
    2. For each shot: generate sample → XOR into result_buf → copy to Row
    3. Allocate NumPy arrays
    4. For each shot: unpack Row bits → write to NumPy buffers
Python ← C++: return NumPy arrays
```

**Problems:**
- N_shots heap allocations for Row objects (each Row allocates a vector<uint64_t>)
- Two-pass: first write to Rows, then unpack to NumPy
- GIL held during Row allocation

### New data path (fused)

```
Python → C++:
    1. Allocate NumPy arrays (single allocation, contiguous memory)
    2. Release GIL
    3. For each shot: generate sample → XOR into result_buf → unpack directly to NumPy
Python ← C++: return NumPy arrays (move semantics, zero-copy)
```

**Improvements:**
- Zero intermediate heap allocations (no Row objects)
- Single-pass: sample → XOR → unpack in one loop iteration
- GIL released for the entire compute phase
- NumPy arrays are contiguous memory (cache-friendly)

The function `generate_many_output_samples_Monte_to_numpy()` already existed in the codebase but was never called from the Python-facing entry point. The old entry point `return_samples_Monte_separate_obs_with_QEPG()` used the two-pass path through `generate_many_output_samples_Monte()` + `convert_bitset_row_to_boolean_separate_obs_numpy()`.

---

## 5. Non-Uniform Noise: Poisson + CDF Binary Search

For non-uniform noise, each source i has independent probabilities (px_i, py_i, pz_i) for X, Y, Z errors. We cannot use Floyd's uniform sampling because sources have different firing probabilities. Instead, we use **inverse CDF sampling** (also known as the **inverse transform method**).

### 5.1 The problem: weighted random selection

We have N noise sources with total firing probabilities:

    p_total_i = px_i + py_i + pz_i,    i = 0, 1, ..., N-1

We need to select a random source where source i is chosen with probability proportional to p_total_i. That is:

    P(select source i) = p_total_i / P_all

where P_all = sum of all p_total_i.

### 5.2 Building the CDF

The **Cumulative Distribution Function** (CDF) is a monotonically non-decreasing array where:

    cdf[0] = 0
    cdf[1] = p_total_0
    cdf[2] = p_total_0 + p_total_1
    cdf[3] = p_total_0 + p_total_1 + p_total_2
    ...
    cdf[N] = p_total_0 + p_total_1 + ... + p_total_{N-1} = P_all

The key property: **cdf[i] is the total probability mass of sources 0 through i-1.** So the interval [cdf[i], cdf[i+1]) has length p_total_i.

**Concrete example.** Suppose N = 4 sources with probabilities:

    p_total_0 = 0.002   (source 0)
    p_total_1 = 0.005   (source 1)
    p_total_2 = 0.001   (source 2)
    p_total_3 = 0.003   (source 3)

Then:

    cdf = [0.000, 0.002, 0.007, 0.008, 0.011]
    P_all = 0.011

The number line [0, 0.011) is partitioned into intervals:

    [0.000, 0.002)  →  source 0  (width 0.002)
    [0.002, 0.007)  →  source 1  (width 0.005)
    [0.007, 0.008)  →  source 2  (width 0.001)
    [0.008, 0.011)  →  source 3  (width 0.003)

### 5.3 The inverse CDF sampling algorithm

To select a random source weighted by probability:

1. Draw u ~ Uniform(0, P_all)
2. Find the interval that contains u: find i such that cdf[i] <= u < cdf[i+1]
3. Return source i

**Why this is correct:** The probability that u lands in interval [cdf[i], cdf[i+1]) equals the length of that interval divided by the total length:

    P(source i selected) = (cdf[i+1] - cdf[i]) / P_all = p_total_i / P_all

This exactly matches the desired weighted distribution.

**Example:** If u = 0.004, we look for the interval containing 0.004:

    cdf[0]=0.000 ≤ 0.004 < cdf[2]=0.007  →  source 1 ✓

    (since [0.002, 0.007) is source 1's interval)

### 5.4 Binary search for efficiency

Finding the correct interval is equivalent to finding:

    src = max { j : cdf[j] <= u }

In code, we use `std::upper_bound(cdf.begin(), cdf.end(), u)` which returns an iterator to the first element strictly greater than u, then subtract 1:

```cpp
double u = to_double_01(rng()) * P_all;              // u ~ Uniform(0, P_all)
auto it = std::upper_bound(cdf.begin(), cdf.end(), u); // first cdf[j] > u
std::size_t src = (it - cdf.begin());                 // index of that element
if (src > 0) --src;                                   // back up one: cdf[src] <= u
```

`std::upper_bound` uses binary search on the sorted CDF array, so it runs in O(log N).

**Why upper_bound and not lower_bound?** Because we want the last index where cdf[j] <= u. `upper_bound` finds the first index where cdf[j] > u, so subtracting 1 gives us the last index where cdf[j] <= u. Using `lower_bound` would incorrectly handle the case where u exactly equals a CDF boundary.

### 5.5 Choosing the Pauli type after selecting the source

Once we know source i fired, we need to decide: was it an X, Y, or Z error? The conditional probabilities are:

    P(X | source i fired) = px_i / p_total_i
    P(Y | source i fired) = py_i / p_total_i
    P(Z | source i fired) = pz_i / p_total_i

We precompute thresholds for each source:

    cond_x[i]  = px_i / p_total_i
    cond_xy[i] = (px_i + py_i) / p_total_i

Then a single uniform random r determines the type:

```cpp
double r = to_double_01(rng());
if      (r < cond_x[src])   type = 1;  // X  (probability px/ptotal)
else if (r < cond_xy[src])  type = 2;  // Y  (probability py/ptotal)
else                         type = 3;  // Z  (probability pz/ptotal)
```

### 5.6 Putting it together: the full sparse non-uniform algorithm

```
// Pre-computation (once per call, not per shot):
For each source i:
    p_total_i = px_i + py_i + pz_i
    cdf[i+1] = cdf[i] + p_total_i
    cond_x[i] = px_i / p_total_i
    cond_xy[i] = (px_i + py_i) / p_total_i
P_all = cdf[N]

// Per shot:
k = sample_poisson(P_all)                    // How many errors this shot?
For each of k errors:
    u = uniform(0, P_all)                    // Where on the number line?
    src = binary_search_cdf(u)              // Which source? O(log N)
    type = choose_pauli(src, cond_x, cond_xy) // X, Y, or Z?
    XOR row[src + (type-1)*N] into result   // Propagate through QEPG
```

### 5.7 Correctness proof

We want to show that this algorithm produces the correct marginal probability for each source.

**Claim:** Under the sparse algorithm, the probability that source i contributes at least one error to a given shot is approximately p_total_i (for small p_total_i).

**Proof:**

Let X_i = number of times source i is selected in a given shot.

The total number of selections is k ~ Poisson(P_all). Each selection independently picks source i with probability p_total_i / P_all.

Therefore, conditioned on k:

    X_i | k ~ Binomial(k, p_total_i / P_all)

Unconditionally (averaging over k ~ Poisson(P_all)):

    X_i ~ Poisson(P_all * p_total_i / P_all) = Poisson(p_total_i)

This is because a Poisson-thinned Poisson process produces independent Poisson counts. Formally: if k ~ Poisson(λ) and each of k items independently falls into bin i with probability q_i, then the count in bin i is Poisson(λ * q_i), and counts across bins are independent.

**So X_i ~ Poisson(p_total_i), independent across sources.**

The probability of at least one error at source i:

    P(X_i >= 1) = 1 - exp(-p_total_i) ≈ p_total_i - p_total_i²/2 + ...

Compare to the exact Bernoulli model: P(fire) = p_total_i.

The difference is:

    |P_sparse - P_exact| = |1 - exp(-p) - p| = p²/2 + O(p³)

For p_total_i = 0.001: error = 5 × 10^{-7}. Negligible.

### 5.8 The double-fire issue in GF(2)

In GF(2) (binary field), XOR is the addition operation. If source i is selected twice, we XOR its row twice, which cancels: row ⊕ row = 0.

This means: if X_i = 2, source i has **no effect** on the outcome (as if it never fired). This is different from the Bernoulli model, where a source fires at most once.

**Impact:** The GF(2)-observable probability of source i contributing is:

    P(X_i is odd) = P(X_i=1) + P(X_i=3) + P(X_i=5) + ...

For X_i ~ Poisson(p):

    P(X_i odd) = (1 - exp(-2p)) / 2 = p - 2p²/3 + ...

Compared to exact: P(fire) = p.

The error is O(p²), which is negligible for p < 0.01. This is why the code enforces `use_sparse = (p_max < 0.01)`.

### 5.9 Cost comparison

| | Dense scan | Sparse CDF |
|--|-----------|-----------|
| Per shot | O(N) random draws | O(k) Poisson + O(k log N) CDF lookups |
| Per shot total | O(N) | O(k log N) where k = E[Np] |
| **When sparse wins** | — | **k log N < N**, i.e., **Np log N < N**, i.e., **p < 1/log N** |

For N = 585 (d=3): log₂(585) ≈ 9.2, so sparse wins when p < 0.11. ✓
For N = 6345 (d=11): log₂(6345) ≈ 12.6, so sparse wins when p < 0.08. ✓

At p = 0.001: k ≈ 0.6 to 6.3, so the CDF path does only ~1-6 binary searches vs. scanning 585-6345 locations.

### 5.10 Benchmark: non-uniform MC (single-threaded)

| d | Stim (M/s) | Non-uniform MC (M/s) | Speedup |
|---|----------:|--------------------:|--------:|
| 3 | 7.20 | 15.8 | **2.2x** |
| 5 | 0.93 | 3.4 | **3.6x** |
| 7 | 0.37 | 1.2 | **3.1x** |
| 9 | 0.21 | 0.4 | **1.9x** |
| 11 | 0.08 | 0.2 | **2.3x** |

The non-uniform path is slower than uniform MC (3-9x vs 2-4x over Stim) because of the O(log N) binary search overhead per error, but still faster than Stim at all distances.

---

## 6. Correctness Argument

### Sparse uniform MC vs Dense Bernoulli

The dense Bernoulli model: each location independently fires with probability p. The joint distribution is:

$$P(\text{config}) = \prod_{i \in S} p \cdot \prod_{i \notin S} (1-p)$$

where S is the set of fired locations.

The sparse Poisson+Floyd model: sample k ~ Poisson(Np), then pick k distinct locations uniformly. The joint distribution is:

$$P(\text{config } S) = \text{Poi}(|S|; Np) \cdot \binom{N}{|S|}^{-1}$$

These differ by O(p²) terms. Specifically:
- Bernoulli: P(exactly k fires) = Binomial(N, p) at k
- Poisson+Floyd: P(exactly k fires) = Poisson(Np) at k
- The Poisson approximation has total variation distance ≤ Np² from the Binomial

For p = 0.001: the per-source error is at most p² = 10^{-6}, completely negligible for Monte Carlo estimation.

### Empirical verification

We verified statistical equivalence at multiple error rates:

```
p=0.001: uniform_MC det=0.01342 obs=0.05390 | nonuniform det=0.01342 obs=0.05269
p=0.010: uniform_MC det=0.11769 obs=0.33888 | nonuniform det=0.11771 obs=0.33887
p=0.050: uniform_MC det=0.36087 obs=0.49683 | nonuniform det=0.36106 obs=0.50029
```

Agreement within Monte Carlo noise at all rates. All 472 existing tests pass.

---

## 7. Benchmark Results

### Before optimization (single-threaded)

| d | N (noise) | Stim (M/s) | QEPG-MC (M/s) | Ratio |
|---|----------|----------:|---------------:|------:|
| 3 | 585 | 7.47 | 0.32 | 0.04x (23x slower) |
| 5 | 1485 | 0.94 | 0.059 | 0.06x (16x slower) |
| 7 | 2745 | 0.35 | 0.020 | 0.06x (17x slower) |
| 9 | 4365 | 0.19 | 0.0088 | 0.05x (22x slower) |
| 11 | 6345 | 0.069 | 0.0049 | 0.07x (14x slower) |

### After optimization (single-threaded, OMP_NUM_THREADS=1)

| d | N (noise) | E[k] = Np | Stim (M/s) | QEPG-MC (M/s) | Speedup vs Stim |
|---|----------|----------|----------:|---------------:|------:|
| 3 | 585 | 0.6 | 6.84 | 20.8 | **3.0x** |
| 5 | 1485 | 1.5 | 0.83 | 7.3 | **8.8x** |
| 7 | 2745 | 2.7 | 0.30 | 2.5 | **8.3x** |
| 9 | 4365 | 4.4 | 0.22 | 0.9 | **4.0x** |
| 11 | 6345 | 6.3 | 0.08 | 0.3 | **4.3x** |

### Improvement factor (QEPG-MC old → new)

| d | Old MC (M/s) | New MC (M/s) | Improvement |
|---|----------:|----------:|------:|
| 3 | 0.32 | 20.8 | **65x** |
| 5 | 0.059 | 7.3 | **124x** |
| 7 | 0.020 | 2.5 | **125x** |
| 9 | 0.0088 | 0.9 | **98x** |
| 11 | 0.0049 | 0.3 | **65x** |

---

## 8. Phase 2 Optimizations: AVX2, Cache Locality, Fused Sampling

Building on the sparse Poisson+Floyd algorithm, a second round of optimizations targeted the inner-loop constants: SIMD width, cache access patterns, and per-shot overhead.

### 8.1 AVX2 SIMD (256-bit XOR)

The XOR propagation in `simd::xor_words()` originally used 128-bit SSE2 instructions (`_mm_xor_si128`), processing 2 × uint64_t per iteration. Switching to 256-bit AVX2 (`_mm256_xor_si256`) processes 4 × uint64_t per iteration, halving the number of SIMD instructions for the XOR propagation step.

**Build flags added:**
- Windows: `/arch:AVX2`
- Linux: `-mavx2 -mbmi2`
- macOS: `-mavx2 -mbmi2`

### 8.2 Sorted Row Access + Software Prefetching

When k > 4 errors fire in a shot, their propagation rows may be scattered across the `FlatBitTable`. Sorting errors by row index converts random memory access into sequential access, improving cache line utilization. Software prefetching (`_mm_prefetch`) loads the next row into L1 cache during the current XOR operation.

```cpp
// Sort row indices for sequential cache access (skip for k <= 4)
if (!sorted)
    std::sort(row_idx_scratch_.begin(), row_idx_scratch_.begin() + sample_size);

// XOR with prefetch
for (std::size_t s = 0; s < sample_size; ++s) {
    if (s + 1 < sample_size)
        _mm_prefetch(flat.row_ptr(row_idx_scratch_[s + 1]), _MM_HINT_T0);
    flat.xor_row_into(row_idx_scratch_[s], result_buf);
}
```

The k <= 4 threshold avoids sort overhead when cache locality doesn't matter (few rows to XOR).

### 8.3 Fused Floyd Insertion + XOR Propagation

The original two-pass approach: (1) generate all k error locations into a scratch buffer, (2) iterate scratch buffer and XOR each row. The fused approach XORs each row immediately upon insertion, eliminating the scratch buffer read pass:

```cpp
while (placed < k) {
    std::size_t pos = gen.bounded(error_size);
    if (!collision_bitmap_[pos]) {
        collision_bitmap_[pos] = 1;
        std::size_t type = 1 + gen.bounded(3);
        flat.xor_row_into(pos + (type - 1) * n_noise, result_buf);  // immediate XOR
        ++placed;
    }
}
```

This eliminates one pass over the sample buffer and keeps the hot data (result_buf) in cache.

### 8.4 Byte-Level Bit Unpacking

The `unpack_result_to_numpy` function extracts individual bits from uint64_t words into uint8_t bytes. The original approach extracted 1 bit per loop iteration with a shift-and-mask. The optimized version processes 8 bits at a time by reinterpreting the word as bytes:

```cpp
for (std::size_t b = 0; b < full_bytes; ++b) {
    std::uint8_t byte = byte_ptr[b];
    det_dst[bit_idx + 0] = (byte     ) & 1;
    det_dst[bit_idx + 1] = (byte >> 1) & 1;
    // ... 8 bits unrolled
    bit_idx += 8;
}
```

This reduces loop iterations by 8x and enables better instruction-level parallelism.

### 8.5 Vose's Alias Method for Non-Uniform Source Selection

The non-uniform MC path previously used O(log N) binary search on a CDF array to select which noise source fired. Vose's alias method replaces this with O(1) per-draw sampling after O(N) preprocessing:

**Preprocessing:** Build an alias table where each bin i has:
- `alias_prob[i]`: acceptance probability
- `alias_idx[i]`: redirect index if rejected

**Per-draw:**
```cpp
std::size_t bin = rng.bounded(num_noise);   // uniform bin selection
double u = to_double_01(rng());             // coin flip
std::size_t src = (u < alias_prob[bin]) ? bin : alias_idx[bin];
```

This converts O(k log N) per shot into O(k) per shot for non-uniform sampling.

### 8.6 Phase 2 Benchmark Results

| d | Stim (M/s) | QEPG-MC Phase 1 (M/s) | QEPG-MC Phase 2 (M/s) | Phase 2 vs Stim |
|---|----------:|---------------------:|---------------------:|------:|
| 3 | 6.00 | 20.8 | 27.94 | **4.7x** |
| 5 | 0.60 | 7.3 | 7.89 | **13.2x** |
| 7 | 0.24 | 2.5 | 2.99 | **12.5x** |
| 9 | 0.13 | 0.9 | 1.07 | **7.9x** |
| 11 | 0.06 | 0.3 | 0.36 | **6.1x** |

Phase 2 improvements are most visible at small d (where the XOR propagation and unpacking dominate) and moderate at large d (where the sparse sampling itself dominates).

---

## 9. Files Modified

| File | Change |
|------|--------|
| `QEPG/src/sampler.hpp` | Added `generate_sample_Monte_sparse()`, `generate_and_xor_sparse()` declarations; rewrote `calculate_parity_output_flat()` with sorted row access + prefetching; added `row_idx_scratch_` member; added SIMD intrinsic includes |
| `QEPG/src/sampler.cpp` | Sparse Poisson+Floyd MC sampler; fused sample+XOR in Monte path; Vose's alias table for non-uniform MC; 8-at-a-time bit unpacking; moved `sample_poisson()` and `to_double_01()` to file scope |
| `QEPG/src/LERcalculator.cpp` | Rewrote `return_samples_Monte_separate_obs_with_QEPG()` to use fused NumPy path with uint8 output |
| `QEPG/src/LERcalculator.hpp` | Updated return type to `pair<array_t<uint8_t>, array_t<uint8_t>>` |
| `QEPG/bindings.cpp` | Updated forward declaration to match new return type |
| `setup.py` | Added AVX2 build flags (`/arch:AVX2`, `-mavx2 -mbmi2`) for all platforms |
| `README.md` | Updated performance tables with new MC numbers |
| `plot_benchmark.py` | Added MC speedup labels to bar plot |

---

## 10. Summary

The optimization exploits a fundamental asymmetry: at realistic error rates (p ~ 0.001), the expected number of errors per shot (k ~ Np) is tiny compared to the total number of noise locations (N). The old algorithm paid O(N) per shot to find O(1) errors. The new algorithm pays O(k) by:

1. Sampling k from Poisson(Np) — O(k) via Knuth's algorithm
2. Picking k distinct locations via Floyd's insertion — O(k) expected
3. XOR-propagating k rows — O(k * W) with SIMD

Phase 2 further optimized the inner-loop constants:

4. AVX2 SIMD — 256-bit XOR width (2x over SSE2)
5. Sorted row access + software prefetching — sequential cache access for k > 4
6. Fused Floyd insertion + XOR — eliminates scratch buffer read pass
7. Byte-level bit unpacking — 8x fewer loop iterations
8. Vose's alias method — O(1) non-uniform source selection (replaces O(log N) CDF binary search)

Combined, this achieves 65-125x improvement over the old QEPG-MC path and **4.7-13.2x speedup over Stim's** highly optimized single-threaded detector sampler.
