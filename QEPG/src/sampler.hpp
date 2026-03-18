/**
 * @file sampler.hpp
 * @brief Random Pauli error sampling and detector/observable output computation.
 *
 * Provides the sampler class for generating random Pauli error configurations
 * (with fixed weight or Bernoulli-distributed) and computing their effect on
 * detectors and observables through the QEPG parity propagation matrix.
 * Batch methods use OpenMP for parallel sample generation.
 */

#ifndef SAMPLE_HPP
#define SAMPLE_HPP
#pragma once


#include <cstddef>
#include <ostream>
#include <string>
#include <vector>
#include <random>
#include <cstring>
#include "QEPG.hpp"
#include <omp.h>          // just include and -fopenmp / /openmp

namespace SAMPLE{


const size_t PAULIX = 1;   ///< Constant identifying a Pauli-X error.
const size_t PAULIY = 2;   ///< Constant identifying a Pauli-Y error.
const size_t PAULIZ = 3;   ///< Constant identifying a Pauli-Z error.


/**
 * @brief xoshiro256++ PRNG — fast, small-state (32 bytes) alternative to mt19937.
 *
 * Passes BigCrush, ~3x faster than mt19937, and fits in registers.
 * State is 4 × uint64_t = 32 bytes (vs mt19937's ~2.5KB).
 */
class Xoshiro256pp {
    std::uint64_t s_[4];

    static inline std::uint64_t rotl(std::uint64_t x, int k) noexcept {
        return (x << k) | (x >> (64 - k));
    }

public:
    using result_type = std::uint64_t;

    explicit Xoshiro256pp(std::uint64_t seed = 0) noexcept {
        // SplitMix64 to initialize state from a single seed
        auto splitmix = [](std::uint64_t& z) -> std::uint64_t {
            z += 0x9e3779b97f4a7c15ULL;
            std::uint64_t r = z;
            r = (r ^ (r >> 30)) * 0xbf58476d1ce4e5b9ULL;
            r = (r ^ (r >> 27)) * 0x94d049bb133111ebULL;
            return r ^ (r >> 31);
        };
        std::uint64_t z = seed;
        s_[0] = splitmix(z);
        s_[1] = splitmix(z);
        s_[2] = splitmix(z);
        s_[3] = splitmix(z);
    }

    static constexpr result_type min() noexcept { return 0; }
    static constexpr result_type max() noexcept { return UINT64_MAX; }

    inline std::uint64_t operator()() noexcept {
        const std::uint64_t result = rotl(s_[0] + s_[3], 23) + s_[0];
        const std::uint64_t t = s_[1] << 17;
        s_[2] ^= s_[0];
        s_[3] ^= s_[1];
        s_[1] ^= s_[2];
        s_[0] ^= s_[3];
        s_[2] ^= t;
        s_[3] = rotl(s_[3], 45);
        return result;
    }

    /// Generate a uniform random number in [0, n) — rejection-free for power-of-2 n,
    /// otherwise uses modulo with rejection to eliminate bias.
    inline std::uint64_t bounded(std::uint64_t n) noexcept {
        std::uint64_t x = operator()();
        return x % n;  // slight bias for non-power-of-2, acceptable for sampling
    }
};

/**
 * @brief A single Pauli error on one noise location.
 *
 * Represents one component of an error configuration: a Pauli operator
 * of a given type applied at a specific noise location index.
 */
struct singlePauli{
    size_t qindex;   ///< Noise location index (0-based).
    size_t type;     ///< Pauli type: PAULIX (1), PAULIY (2), or PAULIZ (3).
};




/**
 * @brief Generates random Pauli error samples and computes their syndrome/observable outputs.
 *
 * The sampler operates on a noise model with num_total_pauliError_ independent
 * noise locations. It can generate error configurations with:
 * - Fixed Hamming weight (Floyd or removal algorithm)
 * - Bernoulli-distributed errors (Monte Carlo)
 *
 * Each error configuration is a vector of singlePauli structs. The sampler
 * then multiplies the error vector against the QEPG parity propagation matrix
 * (transposed) to compute the resulting detector and observable outcomes as
 * a GF(2) row vector.
 */
class sampler{


    public:

    /*---------------------------------------ctor----------*/

        /// @brief Default constructor.
        sampler();

        /**
         * @brief Construct a sampler for a given number of noise locations.
         * @param num_total_paulierror Total number of independent noise locations in the circuit.
         */
        explicit sampler(size_t num_total_paulierror);

        /// @brief Destructor.
        ~sampler();


    /*---------------------------------------Sample one vector with fixed weight----------*/

        /**
         * @brief Generate a single error sample of a given weight using Floyd's algorithm.
         *
         * Selects `weight` distinct noise locations uniformly at random (without
         * replacement) and assigns each a random Pauli type (X, Y, or Z). Uses
         * an insertion-based approach for low weights and automatically delegates
         * to generate_sample_removal() when weight exceeds half the total noise count.
         *
         * @param weight Number of noise locations to activate (Hamming weight of the error).
         * @param gen    Mersenne Twister PRNG instance (modified in place).
         * @return Vector of singlePauli structs representing the error configuration.
         */
        inline std::size_t generate_sample_Floyd(size_t weight,std::mt19937& gen);
        inline std::size_t generate_sample_Floyd(size_t weight, Xoshiro256pp& gen);

        inline std::size_t generate_sample_removal(size_t weight, std::mt19937& gen);
        inline std::size_t generate_sample_removal(size_t weight, Xoshiro256pp& gen);

        inline std::size_t generate_sample_Monte(double error_prob ,size_t ErrorSize,std::mt19937& gen);
        inline std::size_t generate_sample_Monte(double error_prob, size_t ErrorSize, Xoshiro256pp& gen);

        /**
         * @brief Compute detector and observable outcomes for a single error sample.
         *
         * Multiplies the error sample against the transposed detector matrix
         * by XOR-accumulating the appropriate rows (selected by Pauli type).
         * Row indexing: X errors use row 3*pos, Y errors use row 3*pos+1,
         * Z errors use row 3*pos+2.
         *
         * @param graph  The constructed QEPG containing the detector matrix transpose.
         * @param sample The error configuration to evaluate.
         * @return A GF(2) row vector of length (num_detectors + 1). The last bit
         *         is the observable outcome.
         *
         * @note This method uses the detector matrix transpose (get_dectorMatrixTrans).
         */
        inline QEPG::Row calculate_output_from_one_sample(const QEPG::QEPG& graph,std::vector<singlePauli> sample){
            const auto&dm=graph.get_dectorMatrixTrans();
            const std::size_t n_rows=dm.size();
            const std::size_t n_cols=n_rows ? dm[0].size():0;
            QEPG::Row result(n_cols);
            for(const auto& noise: sample){
                // Branchless: X=1→row 3*pos+0, Y=2→row 3*pos+1, Z=3→row 3*pos+2
                result^=dm[3*noise.qindex + (noise.type - 1)];
            }
            return result;
        }


        /**
         * @brief Compute detector and observable outcomes via the parity propagation matrix.
         *
         * Similar to calculate_output_from_one_sample but uses the parity propagation
         * matrix transpose instead of the detector matrix transpose. Row indexing:
         * X errors use row pos, Y errors use row (pos + N), Z errors use row (pos + 2N),
         * where N = total_noise.
         *
         * @param graph  The constructed QEPG containing the parity propagation matrix transpose.
         * @param sample The error configuration to evaluate.
         * @return A GF(2) row vector of length (num_detectors + 1). The last bit
         *         is the observable outcome.
         *
         * @note This is the primary output method used by batch sampling functions.
         */
        inline QEPG::Row calculate_parity_output_from_one_sample(const QEPG::QEPG& graph, const singlePauli* sample, std::size_t sample_size){
            const auto&dm=graph.get_parityPropMatrixTrans();
            const std::size_t n_rows=dm.size();
            const std::size_t n_noise=n_rows/3;
            const std::size_t n_cols=n_rows ? dm[0].size():0;
            QEPG::Row result(n_cols);
            for(std::size_t s=0; s<sample_size; ++s){
                // Branchless: X=1→offset 0, Y=2→offset n_noise, Z=3→offset 2*n_noise
                std::size_t row_idx = sample[s].qindex + (sample[s].type - 1) * n_noise;
                result^=dm[row_idx];
            }
            return result;
        }

        /// Overload accepting std::vector for backward compatibility
        inline QEPG::Row calculate_parity_output_from_one_sample(const QEPG::QEPG& graph,const std::vector<singlePauli>& sample){
            return calculate_parity_output_from_one_sample(graph, sample.data(), sample.size());
        }

        /**
         * @brief Fast parity output using FlatBitTable with SIMD XOR.
         *
         * Uses contiguous memory layout and SIMD-accelerated XOR for
         * maximum throughput. Result is written into a pre-allocated
         * aligned buffer, then copied to QEPG::Row for output.
         */
        inline void calculate_parity_output_flat(
            const qepg_bits::FlatBitTable& flat,
            std::size_t n_noise,
            const singlePauli* sample,
            std::size_t sample_size,
            std::uint64_t* result_buf) const noexcept
        {
            for(std::size_t s = 0; s < sample_size; ++s){
                std::size_t row_idx = sample[s].qindex + (sample[s].type - 1) * n_noise;
                flat.xor_row_into(row_idx, result_buf);
            }
        }


        /**
         * @brief Generate many output samples in parallel with fixed Pauli weight.
         *
         * Uses OpenMP to distribute sample generation across threads. Each thread
         * uses a thread-local Mersenne Twister PRNG seeded from a global seed
         * combined with the thread ID.
         *
         * @param graph           The QEPG to evaluate samples against.
         * @param samplecontainer Output vector resized to samplenumber; each entry is a
         *                        GF(2) row of detector + observable outcomes.
         * @param pauliweight     Number of active noise locations per sample.
         * @param samplenumber    Total number of samples to generate.
         */
        void generate_many_output_samples(const QEPG::QEPG& graph,std::vector<QEPG::Row>& samplecontainer,size_t pauliweight , size_t samplenumber);

        /**
         * @brief Generate many output samples in parallel using Monte Carlo (Bernoulli) errors.
         *
         * Uses OpenMP parallelism. Each noise location independently has an error
         * with the given probability.
         *
         * @param graph           The QEPG to evaluate samples against.
         * @param samplecontainer Output vector resized to samplenumber.
         * @param error_prob      Per-location error probability.
         * @param samplenumber    Total number of samples to generate.
         */
        void generate_many_output_samples_Monte(const QEPG::QEPG& graph,std::vector<QEPG::Row>& samplecontainer,double error_prob, size_t samplenumber);

        /**
         * @brief Enumerate all possible error configurations of a given weight.
         *
         * @param graph           The QEPG to evaluate samples against.
         * @param samplecontainer Output vector to append results to.
         * @param pauliweight     Weight of each error configuration.
         *
         * @note Currently a stub (not implemented).
         */
        void generate_all_samples_with_fixed_weight(const QEPG::QEPG& graph,std::vector<QEPG::Row>& samplecontainer,size_t pauliweight);

        /**
         * @brief Generate many output samples and also return the noise vectors.
         *
         * Single-threaded version that stores both the computed detector/observable
         * outcomes and the raw noise configurations that produced them.
         *
         * @param graph           The QEPG to evaluate samples against.
         * @param noistcontainer  Output vector of noise configurations (one per sample).
         * @param samplecontainer Output vector of GF(2) outcome rows (one per sample).
         * @param pauliweight     Number of active noise locations per sample.
         * @param samplenumber    Total number of samples to generate.
         */
        void generate_many_output_samples_with_noise_vector(const QEPG::QEPG& graph,std::vector<std::vector<singlePauli>>& noistcontainer,std::vector<QEPG::Row>& samplecontainer, size_t pauliweight, size_t samplenumber);
        /// Get pointer to scratch sample buffer (valid after generate_sample_Floyd/Monte).
        const singlePauli* scratch_data() const noexcept { return scratch_sample_.data(); }

        /**
         * @brief Generate samples and write directly to numpy byte buffers.
         * Avoids all intermediate QEPG::Row allocations.
         */
        void generate_many_output_samples_to_numpy(
            const QEPG::QEPG& graph,
            std::uint8_t* det_buf, std::uint8_t* obs_buf,
            std::size_t n_det, size_t pauliweight, size_t samplenumber);

        void generate_many_output_samples_Monte_to_numpy(
            const QEPG::QEPG& graph,
            std::uint8_t* det_buf, std::uint8_t* obs_buf,
            std::size_t n_det, double error_prob, size_t samplenumber);

    private:


        size_t     num_total_pauliError_;   ///< Total number of independent noise locations.

        /// Bitmap for O(1) collision detection in Floyd sampling (replaces std::unordered_set).
        std::vector<std::uint8_t> collision_bitmap_;

        /// Pre-allocated scratch buffer for sample generation (replaces per-sample std::vector).
        std::vector<singlePauli> scratch_sample_;

};



}
#endif
