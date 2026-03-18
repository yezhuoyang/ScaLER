/**
 * @file simd_util.hpp
 * @brief SIMD-accelerated XOR and aligned memory utilities.
 *
 * Provides compile-time SIMD dispatch for XOR operations on uint64_t arrays
 * and cross-platform aligned memory allocation. Supports AVX2, SSE2, ARM NEON,
 * and scalar fallback.
 */

#ifndef SIMD_UTIL_HPP
#define SIMD_UTIL_HPP

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>

// SIMD intrinsic headers
#if defined(__AVX2__)
    #include <immintrin.h>
#elif defined(__SSE2__) || (defined(_MSC_VER) && (defined(_M_X64) || defined(_M_AMD64)))
    #include <emmintrin.h>
    #define QEPG_HAS_SSE2 1
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
    #include <arm_neon.h>
    #define QEPG_HAS_NEON 1
#endif

namespace simd {

// ---------------------------------------------------------------
// Aligned memory allocation (64-byte for cache line alignment)
// ---------------------------------------------------------------

inline std::uint64_t* aligned_alloc_u64(std::size_t count) {
    if (count == 0) return nullptr;
    const std::size_t bytes = count * sizeof(std::uint64_t);

#if defined(_MSC_VER)
    void* ptr = _aligned_malloc(bytes, 64);
#elif defined(__APPLE__) || defined(__linux__)
    void* ptr = nullptr;
    if (posix_memalign(&ptr, 64, bytes) != 0) ptr = nullptr;
#else
    void* ptr = std::aligned_alloc(64, ((bytes + 63) / 64) * 64);
#endif
    return static_cast<std::uint64_t*>(ptr);
}

inline void aligned_free(void* ptr) {
    if (!ptr) return;
#if defined(_MSC_VER)
    _aligned_free(ptr);
#else
    std::free(ptr);
#endif
}

// ---------------------------------------------------------------
// SIMD XOR: dst[i] ^= src[i] for i in [0, n)
// ---------------------------------------------------------------

inline void xor_words(std::uint64_t* dst, const std::uint64_t* src, std::size_t n) noexcept {
#if defined(__AVX2__)
    // AVX2: 256-bit = 4 × uint64_t per instruction
    std::size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256i a = _mm256_load_si256(reinterpret_cast<const __m256i*>(dst + i));
        __m256i b = _mm256_load_si256(reinterpret_cast<const __m256i*>(src + i));
        _mm256_store_si256(reinterpret_cast<__m256i*>(dst + i), _mm256_xor_si256(a, b));
    }
    for (; i < n; ++i) dst[i] ^= src[i];

#elif defined(QEPG_HAS_SSE2)
    // SSE2: 128-bit = 2 × uint64_t per instruction
    std::size_t i = 0;
    for (; i + 2 <= n; i += 2) {
        __m128i a = _mm_load_si128(reinterpret_cast<const __m128i*>(dst + i));
        __m128i b = _mm_load_si128(reinterpret_cast<const __m128i*>(src + i));
        _mm_store_si128(reinterpret_cast<__m128i*>(dst + i), _mm_xor_si128(a, b));
    }
    for (; i < n; ++i) dst[i] ^= src[i];

#elif defined(QEPG_HAS_NEON)
    // ARM NEON: 128-bit = 2 × uint64_t per instruction
    std::size_t i = 0;
    for (; i + 2 <= n; i += 2) {
        uint64x2_t a = vld1q_u64(dst + i);
        uint64x2_t b = vld1q_u64(src + i);
        vst1q_u64(dst + i, veorq_u64(a, b));
    }
    for (; i < n; ++i) dst[i] ^= src[i];

#else
    // Scalar fallback
    for (std::size_t i = 0; i < n; ++i)
        dst[i] ^= src[i];
#endif
}

// ---------------------------------------------------------------
// SIMD zero: memset wrapper for aligned buffers
// ---------------------------------------------------------------

inline void zero_words(std::uint64_t* dst, std::size_t n) noexcept {
    std::memset(dst, 0, n * sizeof(std::uint64_t));
}

} // namespace simd

#endif // SIMD_UTIL_HPP
