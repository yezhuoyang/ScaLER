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
#include <unordered_set>
#include "QEPG.hpp"
#include <omp.h>          // just include and -fopenmp / /openmp

namespace SAMPLE{


const size_t PAULIX = 1;   ///< Constant identifying a Pauli-X error.
const size_t PAULIY = 2;   ///< Constant identifying a Pauli-Y error.
const size_t PAULIZ = 3;   ///< Constant identifying a Pauli-Z error.

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
        inline std::vector<singlePauli> generate_sample_Floyd(size_t weight,std::mt19937& gen);

        /**
         * @brief Generate a single error sample using removal-based sampling.
         *
         * Starts with all noise locations active and randomly removes positions
         * until only `weight` remain. More efficient than Floyd's algorithm when
         * the weight is close to the total number of noise locations.
         *
         * @param weight Number of noise locations to keep active.
         * @param gen    Mersenne Twister PRNG instance (modified in place).
         * @return Vector of singlePauli structs representing the error configuration.
         */
        inline std::vector<singlePauli> generate_sample_removal(size_t weight, std::mt19937& gen);

        /**
         * @brief Generate a single error sample using Bernoulli (Monte Carlo) sampling.
         *
         * Each of the first ErrorSize noise locations independently suffers an
         * error with probability error_prob. Activated locations receive a
         * uniformly random Pauli type (X, Y, or Z).
         *
         * @param error_prob Probability that each noise location is activated.
         * @param ErrorSize  Number of noise locations to consider.
         * @param gen        Mersenne Twister PRNG instance (modified in place).
         * @return Vector of singlePauli structs for all activated locations.
         */
        inline std::vector<singlePauli> generate_sample_Monte(double error_prob ,size_t ErrorSize,std::mt19937& gen);

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
            for(singlePauli noise: sample){
                size_t pos=noise.qindex;
                size_t type=noise.type;
                if(type==SAMPLE::PAULIX){
                    result^=dm[3*pos];
                }
                else if(type==SAMPLE::PAULIY){
                    result^=dm[3*pos+1];
                }
                else if(type==SAMPLE::PAULIZ){
                    result^=dm[3*pos+2];
                }
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
        inline QEPG::Row calculate_parity_output_from_one_sample(const QEPG::QEPG& graph,const std::vector<singlePauli>& sample){
            const auto&dm=graph.get_parityPropMatrixTrans();
            const std::size_t n_rows=dm.size();
            const std::size_t n_noise=int(n_rows/3);
            const std::size_t n_cols=n_rows ? dm[0].size():0;
            QEPG::Row result(n_cols);
            for(singlePauli noise: sample){
                size_t pos=noise.qindex;
                size_t type=noise.type;
                if(type==SAMPLE::PAULIX){
                    result^=dm[pos];
                }
                else if(type==SAMPLE::PAULIY){
                    result^=dm[pos+n_noise];
                }
                else if(type==SAMPLE::PAULIZ){
                    result^=dm[pos+n_noise*2];
                }
            }
            return result;
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
    private:


        size_t     num_total_pauliError_;   ///< Total number of independent noise locations.



};



}
#endif
