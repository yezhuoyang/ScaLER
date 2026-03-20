/**
 * @file hotspot.hpp
 * @brief C++ hotspot analysis: labeled sampling + decoder-agnostic aggregation.
 *
 * Two-phase design that separates sampling from decoding:
 *
 * Phase 1 — `labeled_sample()`:
 *   Sample non-uniform noise via QEPG, return detector/observable outcomes
 *   AND a per-shot category bitmask (uint8, bit i = 1 if any Type-i noise fired).
 *   This is the same sampling as `return_samples_nonuniform_to_numpy` but with
 *   the bitmask tracking added at near-zero cost.
 *
 * Phase 2 — User decodes in Python with ANY decoder (pymatching, BPOSD, custom):
 *   logical_error = (observable_outcomes != decoder.decode_batch(detector_outcomes))
 *
 * Phase 3 — `hotspot_aggregate()`:
 *   Takes per-shot category bitmasks + decoded logical_error mask, computes
 *   P(category fired | logical error), lift, configuration breakdown.
 *   Pure aggregation — no sampling, no decoding.
 *
 * Uses OpenMP parallelism, SIMD XOR accumulation, and Xoshiro256++ RNG.
 */

#ifndef HOTSPOT_HPP
#define HOTSPOT_HPP
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>
#include <array>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <atomic>
#include <thread>
#include <mutex>
#include <cmath>

#include "QEPG.hpp"
#include "sampler.hpp"
#include "noiselabel.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace hotspot {

/**
 * @brief Per-category conditional probability report.
 */
struct CategoryReport {
    std::uint8_t type;               ///< NoiseType enum value
    std::string  name;               ///< Category name string
    std::size_t  count;              ///< Number of noise sources with this label
    double       p_fired;            ///< P(any source in category fired) — marginal
    double       p_fired_given_error;///< P(any source in category fired | logical error)
    double       p_error_given_fired;///< P(logical error | any source in category fired)
    double       lift;               ///< p_fired_given_error / p_fired
};

/**
 * @brief Multi-error configuration entry.
 */
struct ConfigEntry {
    std::string configuration;       ///< Sorted " + "-joined category names
    std::size_t count;               ///< Number of logical-error shots with this config
    double      fraction;            ///< Fraction of all logical errors
};

/**
 * @brief Per-source ranking entry.
 */
struct SourceReport {
    std::size_t noise_index;
    std::size_t qubit;
    std::uint8_t type;
    double p_fired;
    double p_fired_given_error;
    double lift;
};

/**
 * @brief Complete hotspot analysis result.
 */
struct HotspotResult {
    std::size_t num_shots;
    std::size_t num_logical_errors;
    std::size_t num_noise;

    std::vector<CategoryReport> category_reports;   ///< Sorted by lift (desc)
    std::vector<ConfigEntry>    config_reports;      ///< Sorted by count (desc)
    std::vector<SourceReport>   source_reports;      ///< Top-N sorted by lift (desc)
};


/**
 * @brief Phase 1: Sample non-uniform noise with per-shot category bitmask tracking.
 *
 * Same noise sampling as `return_samples_nonuniform_to_numpy`, but additionally
 * writes a per-shot category bitmask (uint8) where bit i = 1 if any noise source
 * of Type i fired on that shot. The bitmask overhead is negligible.
 *
 * After this call, the user decodes detector_outcomes with their chosen decoder
 * to determine logical errors, then passes bitmasks + logical_error mask to
 * `hotspot_aggregate()`.
 *
 * @param graph           Pre-compiled QEPG graph.
 * @param label_map       Noise label map (from auto_label).
 * @param det_buf         Output buffer for detector outcomes, shape (num_shots, n_det).
 * @param obs_buf         Output buffer for observable outcomes, shape (num_shots,).
 * @param bitmask_buf     Output buffer for category bitmasks, shape (num_shots,).
 *                        Bit 0 = data_qubit_error, bit 1 = ghost_error,
 *                        bit 2 = hook_error, bit 3 = measurement_error.
 * @param n_det           Number of detectors.
 * @param noise_probs     Flat array of shape (num_noise, 3): [px0, py0, pz0, ...].
 * @param num_noise       Number of independent noise sources.
 * @param corr_pairs      Array of correlated pairs, or nullptr.
 * @param num_corr_pairs  Number of correlated pairs.
 * @param num_shots       Total number of Monte Carlo shots.
 */
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
    std::size_t num_shots);


/**
 * @brief Phase 3: Aggregate hotspot statistics from bitmasks + decoded logical error mask.
 *
 * Pure aggregation — no sampling, no decoding. Takes the per-shot category
 * bitmasks (from labeled_sample) and a logical_error mask (from user's decoder),
 * and computes conditional probabilities and configuration breakdown.
 *
 * @param bitmask_buf     Per-shot category bitmasks, shape (num_shots,). From labeled_sample.
 * @param logical_errors  Per-shot logical error mask, shape (num_shots,). 1 = logical error.
 *                        Computed by user as: (observable != decoder.decode_batch(detectors)).
 * @param num_shots       Number of shots.
 * @param label_map       Noise label map (for category counts in report).
 * @param min_config_count Minimum count for configuration report entries.
 * @return HotspotResult with category_reports and config_reports.
 *         source_reports will be empty (per-source requires per-source fired info).
 */
HotspotResult hotspot_aggregate(
    const std::uint8_t* bitmask_buf,
    const std::uint8_t* logical_errors,
    std::size_t num_shots,
    const noiselabel::NoiseLabelMap& label_map,
    std::size_t min_config_count = 5);


}  // namespace hotspot

#endif  // HOTSPOT_HPP
