# ScaLERQEC Code Review Report

**Date:** 2026-03-19
**Reviewer:** Claude Opus 4.6
**Scope:** Full codebase review — C++ backend (QEPG/), Python frontend (src/scalerqec/), tests

---

## Summary

The project implements a sophisticated framework for estimating logical error rates (LER) of quantum error-correcting codes via stratified sampling and Monte Carlo methods, with a high-performance C++ backend for error propagation graph construction and sampling. The codebase is functional and demonstrates strong algorithmic design, but has several correctness bugs, performance bottlenecks, and maintainability issues that should be addressed to achieve production quality.

**Findings by severity:**
- P0 (Bugs): 4
- P1 (Correctness/Performance): 10
- P2 (Maintainability/Minor): 9

---

## P0 — Bugs (Must Fix)

### 1. `add_ZError` writes wrong gate name
**File:** `QEPG/src/clifford.cpp:45`
**Issue:** `"X_ZRROR"` should be `"Z_ERROR"`. This typo means Z errors added via this function will not match any case in backward propagation and will be silently ignored.

### 2. Static RNG seed causes correlated samples
**File:** `QEPG/src/sampler.cpp:203`
**Issue:** `static const uint64_t global_seed` is evaluated once per program lifetime. Combined with thread-ID hash (also constant per thread), repeated calls to `generate_many_output_samples` from the same thread produce identical RNG sequences. This violates Monte Carlo independence.
**Fix:** Add `static std::atomic<uint64_t> call_counter{0}` and XOR `call_counter.fetch_add(1)` into the seed.

### 3. Signed int overflow in backward traversal
**File:** `QEPG/src/QEPG.cpp:160`
**Issue:** `for(int t=gate_size-1;t>=0;t--)` — if `gate_size==0`, `gate_size-1` wraps to `SIZE_MAX`, which cast to `int` is undefined behavior.
**Fix:** `for(size_t t = gate_size; t-- > 0;)`

### 4. Missing `continue` after Reset
**File:** `QEPG/src/QEPG.cpp:214`
**Issue:** The Reset block lacks `continue`, causing fall-through to subsequent gate checks. Correct by accident now, but fragile.

---

## P1 — Correctness & Performance

### 5. String comparison in hot backward traversal loop
**File:** `QEPG/src/QEPG.cpp:163`
**Issue:** Creates `std::string` copy + 7 string comparisons per gate. For d=17 surface codes with ~100K gates, this is a major bottleneck.
**Fix:** Replace `Gate::name` with `enum GateKind` and use `switch`.

### 6. Circuit copied by value into QEPG constructor
**File:** `QEPG/src/QEPG.cpp:24`
**Issue:** Copies entire gate vector, parity groups, measurement mappings.
**Fix:** Take `const&` or `&&`.

### 7. Gate struct uses heap-allocated strings and vectors
**File:** `QEPG/src/clifford.hpp:29`
**Issue:** 2 heap allocations per Gate. For large circuits, this dominates memory.
**Fix:** `struct Gate { GateKind kind; size_t qubits[2]; uint8_t num_qubits; };`

### 8. Linear search in observable list
**File:** `QEPG/src/QEPG.cpp:197`
**Issue:** `std::find` on every measurement. Convert to `unordered_set`.

### 9. No size assertion in DynamicBitset bitwise ops
**File:** `QEPG/src/dynamic_bitset.hpp:143`
**Issue:** `operator^=` reads `rhs.blocks_` without checking sizes match. UB if different.

### 10. FlatBitTable missing null check after allocation
**File:** `QEPG/src/flat_bit_table.hpp:39`
**Issue:** `aligned_alloc` can return nullptr; `memset(nullptr, ...)` is UB.

### 11. Sparse sampler double-XOR cancellation
**File:** `QEPG/src/sampler.cpp:518`
**Issue:** Poisson+CDF path can sample same source twice, causing XOR cancellation (error vanishes). Creates O(p^2) bias.

### 12. Noise accumulation by addition instead of channel composition
**File:** `src/scalerqec/Monte/noise_model_parser.py:329`
**Issue:** `p_combined = p1 + p2` instead of `p1 + p2 - 2*p1*p2`. Significant at p > 0.01.

### 13. DETECTOR prefix matching fails with spaces
**File:** `src/scalerqec/Clifford/stimparser.py:49`
**Issue:** `"DETECTOR("` prefix fails for `"DETECTOR (0,1) rec[-1]"`.

### 14. Missing PAULI_CHANNEL_1/2 handling
**File:** `src/scalerqec/Clifford/stimparser.py:52`
**Issue:** These Stim directives pass through as unknown gates.

---

## P2 — Maintainability

### 15. 5x copy-pasted adaptive batching loop
**File:** `src/scalerqec/Monte/monteLER.py`
**Fix:** Extract `_adaptive_monte_carlo(sample_fn, decode_fn)`.

### 16. Dual Python/C++ QEPG implementations
**Files:** `QEPGpython.py` + `QEPG.cpp`
**Risk:** Must be kept in sync manually when adding gates.

### 17. Wildcard imports in monteLER.py
**Fix:** `from ..Clifford.clifford import CliffordCircuit`

### 18. Dead member variables in QEPG class
**File:** `QEPG/src/QEPG.hpp:262,268-270`
**Issue:** `COLS`, `detectorMatrix_`, `detectorMatrixTranspose_` never used.

### 19. _GATE_NOISE_COUNT duplicates stimparser tables
**File:** `src/scalerqec/Monte/noise_model_parser.py:68`
**Fix:** Derive automatically from decomposition tables.

### 20. Gate dispatch table rebuilt per loop iteration
**File:** `src/scalerqec/Clifford/clifford.py:536`

### 21. Scaler.py is 2130 lines
**Fix:** Split into SamplingStrategy, ModelFitter, Plotter, IterationLogger.

### 22. No CMake build system for C++ backend
**Fix:** Add CMakeLists.txt for cross-platform builds.

### 23. Magic numbers in sparse/dense threshold
**File:** `QEPG/src/sampler.cpp:479`
**Fix:** Document the 0.01 and 0.06 thresholds with bias analysis.

---

## Test Coverage Gaps

1. No round-trip test: C++ QEPG vs Python QEPG on arbitrary circuits
2. No test for CZ gate handling in noise_model_parser
3. No test for multiple observables
4. No stress test for large circuits (d >= 21)
5. No statistical equivalence test at sparse/dense boundary (p=0.01)

---

## Architecture Recommendations

1. **Single source of truth for gate handling:** Make C++ authoritative, deprecate Python QEPG, or generate one from the other.
2. **Shared gate decomposition table:** Unify `_1Q_DECOMPOSITIONS` and `_GATE_NOISE_COUNT` into a single data structure.
3. **Strategy pattern for decoders:** Replace ScalerLDPC/ScalerOptimized inheritance with injected decoder strategy.
4. **Build system:** Add CMakeLists.txt with proper AVX2/SSE2 detection, OpenMP support, and cross-platform CI.
