/**
 * @file hotspot.cpp
 * @brief Implementation of labeled sampling + decoder-agnostic hotspot aggregation.
 *
 * Two-phase design:
 *
 * Phase 1 (labeled_sample): Extends non-uniform sampling to also track a per-shot
 * category bitmask (uint8). Bit i = 1 if any noise source of Type i fired.
 * The bitmask is computed on-the-fly during sampling at near-zero cost — we
 * already iterate over fired sources to XOR propagation rows, so tracking
 * the category bitmask is just one extra OR per fired source.
 *
 * Phase 3 (hotspot_aggregate): Pure aggregation over the bitmask + logical_error
 * arrays. No sampling, no decoding. Uses OpenMP for large shot counts.
 */

#include "hotspot.hpp"
#include <sstream>
#include <random>
#include <numeric>

namespace hotspot {

/// Convert uint64 to uniform double in [0, 1).
static inline double to_double_01(std::uint64_t v) noexcept {
    return (v >> 11) * (1.0 / (std::uint64_t{1} << 53));
}

/// Poisson sampling for sparse path.
static inline std::size_t sample_poisson(SAMPLE::Xoshiro256pp& rng, double lambda) noexcept {
    if (lambda < 30.0) {
        double L = std::exp(-lambda);
        std::size_t k = 0;
        double p = 1.0;
        do {
            ++k;
            p *= to_double_01(rng());
        } while (p > L);
        return k - 1;
    } else {
        double u1 = to_double_01(rng());
        double u2 = to_double_01(rng());
        if (u1 < 1e-300) u1 = 1e-300;
        constexpr double two_pi = 2.0 * 3.14159265358979323846;
        double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(two_pi * u2);
        double x = lambda + std::sqrt(lambda) * z;
        return x > 0.0 ? static_cast<std::size_t>(x + 0.5) : 0;
    }
}

/// Unpack result_buf (uint64_t words) into numpy byte buffers.
static inline void unpack_result_to_numpy(
    const std::uint64_t* result_buf,
    std::uint8_t* det_dst,
    std::uint8_t& obs_dst,
    std::size_t n_det) noexcept
{
    constexpr std::size_t WORD_BITS = 64;
    const std::size_t full_words = n_det / WORD_BITS;

    for (std::size_t b = 0; b < full_words; ++b) {
        std::uint64_t word = result_buf[b];
        for (std::size_t k = 0; k < WORD_BITS; ++k, word >>= 1)
            det_dst[b * WORD_BITS + k] = static_cast<std::uint8_t>(word & 1);
    }

    const std::size_t rem = n_det % WORD_BITS;
    if (rem) {
        std::uint64_t word = result_buf[full_words];
        for (std::size_t k = 0; k < rem; ++k, word >>= 1)
            det_dst[full_words * WORD_BITS + k] = static_cast<std::uint8_t>(word & 1);
    }

    obs_dst = static_cast<std::uint8_t>(
        (result_buf[n_det / WORD_BITS] >> (n_det % WORD_BITS)) & 1);
}


/// Build configuration key name from bitmask.
static std::string config_key_name(std::uint8_t mask) {
    if (mask == 0) return "(no noise)";

    static const char* names[] = {
        "data_qubit_error",
        "ghost_error",
        "hook_error",
        "measurement_error"
    };

    std::string result;
    for (int i = 0; i < 4; ++i) {
        if (mask & (1 << i)) {
            if (!result.empty()) result += " + ";
            result += names[i];
        }
    }
    return result;
}


// =========================================================================
// Phase 1: labeled_sample — sampling with category bitmask tracking
// =========================================================================

void labeled_sample(
    const QEPG::QEPG& graph,
    const noiselabel::NoiseLabelMap& label_map,
    std::uint8_t* det_buf,
    std::uint8_t* obs_buf,
    std::uint8_t* bitmask_buf,
    std::size_t n_det,
    const double* noise_probs,
    std::size_t num_noise,
    const SAMPLE::CorrelatedPair* corr_pairs,
    std::size_t num_corr_pairs,
    std::size_t num_shots)
{
    const auto& flat = graph.get_parityPropMatrixTransFlat();
    const std::size_t n_noise_flat = flat.n_rows() / 3;
    const std::size_t n_words = flat.words_per_row();

    // Pre-compute per-source total probability and CDF
    std::vector<double> ptotal(num_noise);
    std::vector<double> cdf(num_noise + 1, 0.0);
    for (std::size_t i = 0; i < num_noise; ++i) {
        ptotal[i] = noise_probs[3 * i] + noise_probs[3 * i + 1] + noise_probs[3 * i + 2];
        cdf[i + 1] = cdf[i] + ptotal[i];
    }
    double P_all = cdf[num_noise];

    // Pre-compute conditional Pauli type thresholds
    std::vector<double> cond_x(num_noise, 0.0);
    std::vector<double> cond_xy(num_noise, 0.0);
    for (std::size_t i = 0; i < num_noise; ++i) {
        if (ptotal[i] > 0.0) {
            cond_x[i]  = noise_probs[3 * i] / ptotal[i];
            cond_xy[i] = (noise_probs[3 * i] + noise_probs[3 * i + 1]) / ptotal[i];
        }
    }

    // Decide sparse vs dense
    double p_max = 0.0;
    for (std::size_t i = 0; i < num_noise; ++i) {
        if (ptotal[i] > p_max) p_max = ptotal[i];
    }
    const bool use_sparse = (p_max < 0.01) && (P_all < 0.06 * num_noise);

    // Dense path: uint32 thresholds
    std::vector<std::uint32_t> thresh;
    if (!use_sparse) {
        thresh.resize(3 * num_noise);
        for (std::size_t i = 0; i < num_noise; ++i) {
            double px = noise_probs[3 * i];
            double py = noise_probs[3 * i + 1];
            double pz = noise_probs[3 * i + 2];
            thresh[3 * i]     = static_cast<std::uint32_t>(px * 4294967296.0);
            thresh[3 * i + 1] = static_cast<std::uint32_t>((px + py) * 4294967296.0);
            thresh[3 * i + 2] = static_cast<std::uint32_t>((px + py + pz) * 4294967296.0);
        }
    }

    // Correlated pair thresholds
    std::vector<std::uint32_t> corr_thresh(num_corr_pairs);
    for (std::size_t i = 0; i < num_corr_pairs; ++i) {
        corr_thresh[i] = static_cast<std::uint32_t>(corr_pairs[i].prob * 4294967296.0);
    }

    // Pre-build noise_index -> type bitmask lookup
    std::vector<std::uint8_t> noise_type_bit(num_noise);
    for (std::size_t i = 0; i < num_noise; ++i) {
        noise_type_bit[i] = static_cast<std::uint8_t>(
            1 << static_cast<std::uint8_t>(label_map.get_type(i)));
    }

    // RNG seeding
    static const std::uint64_t global_seed = std::random_device{}();
    static std::atomic<uint64_t> call_counter{0};
    const uint64_t call_id = call_counter.fetch_add(1, std::memory_order_relaxed);

    #pragma omp parallel
    {
        SAMPLE::Xoshiro256pp rng(global_seed
            ^ std::hash<std::thread::id>{}(std::this_thread::get_id())
            ^ (call_id * 0x9E3779B97F4A7C15ULL));

        std::uint64_t* result_buf = simd::aligned_alloc_u64(flat.stride_words());

        #pragma omp for schedule(static)
        for (long long shot = 0; shot < static_cast<long long>(num_shots); ++shot) {
            simd::zero_words(result_buf, n_words);
            std::uint8_t cat_mask = 0;

            if (use_sparse) {
                std::size_t num_errors = sample_poisson(rng, P_all);
                for (std::size_t e = 0; e < num_errors; ++e) {
                    double u = to_double_01(rng()) * P_all;
                    auto it = std::upper_bound(cdf.begin(), cdf.end(), u);
                    std::size_t src = static_cast<std::size_t>(it - cdf.begin());
                    if (src > 0) --src;
                    if (src >= num_noise) src = num_noise - 1;

                    double r = to_double_01(rng());
                    std::size_t type;
                    if (r < cond_x[src])       type = 1;
                    else if (r < cond_xy[src]) type = 2;
                    else                       type = 3;

                    flat.xor_row_into(src + (type - 1) * n_noise_flat, result_buf);
                    cat_mask |= noise_type_bit[src];
                }
            } else {
                for (std::size_t s = 0; s < num_noise; ++s) {
                    std::uint32_t r = static_cast<std::uint32_t>(rng() >> 32);
                    std::uint32_t t_xyz = thresh[3 * s + 2];
                    if (r < t_xyz) {
                        std::size_t type;
                        if (r < thresh[3 * s])          type = 1;
                        else if (r < thresh[3 * s + 1]) type = 2;
                        else                            type = 3;
                        flat.xor_row_into(s + (type - 1) * n_noise_flat, result_buf);
                        cat_mask |= noise_type_bit[s];
                    }
                }
            }

            // Correlated pairs
            for (std::size_t c = 0; c < num_corr_pairs; ++c) {
                std::uint32_t r = static_cast<std::uint32_t>(rng() >> 32);
                if (r < corr_thresh[c]) {
                    std::size_t pidx = rng.bounded(15);
                    std::size_t pa = SAMPLE::TWO_QUBIT_PAULIS[pidx][0];
                    std::size_t pb = SAMPLE::TWO_QUBIT_PAULIS[pidx][1];
                    if (pa > 0) {
                        flat.xor_row_into(
                            corr_pairs[c].source_a + (pa - 1) * n_noise_flat,
                            result_buf);
                        if (corr_pairs[c].source_a < num_noise)
                            cat_mask |= noise_type_bit[corr_pairs[c].source_a];
                    }
                    if (pb > 0) {
                        flat.xor_row_into(
                            corr_pairs[c].source_b + (pb - 1) * n_noise_flat,
                            result_buf);
                        if (corr_pairs[c].source_b < num_noise)
                            cat_mask |= noise_type_bit[corr_pairs[c].source_b];
                    }
                }
            }

            // Unpack detector/observable to numpy buffers
            unpack_result_to_numpy(result_buf,
                det_buf + shot * n_det, obs_buf[shot], n_det);

            // Write category bitmask
            bitmask_buf[shot] = cat_mask;
        }

        simd::aligned_free(result_buf);
    }
}


// =========================================================================
// Phase 3: hotspot_aggregate — pure aggregation, decoder-agnostic
// =========================================================================

HotspotResult hotspot_aggregate(
    const std::uint8_t* bitmask_buf,
    const std::uint8_t* logical_errors,
    std::size_t num_shots,
    const noiselabel::NoiseLabelMap& label_map,
    std::size_t min_config_count)
{
    // Per-category accumulators
    std::array<std::size_t, 4> cat_fired       = {0, 0, 0, 0};
    std::array<std::size_t, 4> cat_fired_error = {0, 0, 0, 0};
    std::unordered_map<std::uint8_t, std::size_t> config_counts;
    std::size_t num_errors = 0;

    // Single-pass aggregation (fast — just array lookups)
    // Use OpenMP reduction for large shot counts
    #pragma omp parallel if(num_shots > 100000)
    {
        std::array<std::size_t, 4> local_cat_fired       = {0, 0, 0, 0};
        std::array<std::size_t, 4> local_cat_fired_error = {0, 0, 0, 0};
        std::unordered_map<std::uint8_t, std::size_t> local_config;
        std::size_t local_errors = 0;

        #pragma omp for schedule(static)
        for (long long i = 0; i < static_cast<long long>(num_shots); ++i) {
            std::uint8_t mask = bitmask_buf[i];
            bool is_error = logical_errors[i] != 0;

            for (int t = 0; t < 4; ++t) {
                if (mask & (1 << t)) {
                    local_cat_fired[t]++;
                    if (is_error) local_cat_fired_error[t]++;
                }
            }

            if (is_error) {
                local_errors++;
                local_config[mask]++;
            }
        }

        #pragma omp critical
        {
            for (int t = 0; t < 4; ++t) {
                cat_fired[t]       += local_cat_fired[t];
                cat_fired_error[t] += local_cat_fired_error[t];
            }
            for (auto& [mask, cnt] : local_config) {
                config_counts[mask] += cnt;
            }
            num_errors += local_errors;
        }
    }

    // Build result
    HotspotResult result;
    result.num_shots = num_shots;
    result.num_logical_errors = num_errors;
    result.num_noise = label_map.size();

    auto cat_counts = label_map.category_counts();

    // Category reports
    for (int t = 0; t < 4; ++t) {
        if (cat_counts[t] == 0) continue;

        CategoryReport cr;
        cr.type = static_cast<std::uint8_t>(t);
        cr.name = noiselabel::noise_type_name(static_cast<noiselabel::NoiseType>(t));
        cr.count = cat_counts[t];

        double n_fired = static_cast<double>(cat_fired[t]);
        double n_fired_error = static_cast<double>(cat_fired_error[t]);
        double n_err = static_cast<double>(num_errors);
        double n_shots = static_cast<double>(num_shots);

        cr.p_fired = n_shots > 0 ? n_fired / n_shots : 0.0;
        cr.p_fired_given_error = n_err > 0 ? n_fired_error / n_err : 0.0;
        cr.p_error_given_fired = n_fired > 0 ? n_fired_error / n_fired : 0.0;
        cr.lift = cr.p_fired > 0 ? cr.p_fired_given_error / cr.p_fired : 1e9;

        result.category_reports.push_back(cr);
    }

    std::sort(result.category_reports.begin(), result.category_reports.end(),
        [](const CategoryReport& a, const CategoryReport& b) { return a.lift > b.lift; });

    // Configuration reports
    for (auto& [mask, cnt] : config_counts) {
        if (cnt < min_config_count) continue;
        ConfigEntry ce;
        ce.configuration = config_key_name(mask);
        ce.count = cnt;
        ce.fraction = num_errors > 0
            ? static_cast<double>(cnt) / static_cast<double>(num_errors)
            : 0.0;
        result.config_reports.push_back(ce);
    }
    std::sort(result.config_reports.begin(), result.config_reports.end(),
        [](const ConfigEntry& a, const ConfigEntry& b) { return a.count > b.count; });

    return result;
}


}  // namespace hotspot
