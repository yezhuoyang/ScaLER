/**
 * @file LERcalculator.hpp
 * @brief High-level sampling API exposed to Python via pybind11.
 *
 * Provides functions that combine circuit compilation, QEPG construction,
 * and error sampling into single calls. These functions are the primary
 * interface used by the Python ScaLERQEC package. They accept STIM program
 * strings (or pre-compiled QEPG graphs) and return detector/observable
 * outcomes as vectors of booleans or NumPy arrays.
 */

#ifndef LER_HPP
#define LER_HPP
#pragma once

#include "QEPG.hpp"
#include "sampler.hpp"
#include <chrono>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>          // <-- defines py::array_t
namespace LERcalculator{

namespace py = pybind11;

/**
 * @brief A single noise event with Pauli type and noise location.
 *
 * Used internally for converting singlePauli structs to a format
 * suitable for Python consumption.
 */
 struct noisesample{
      int type;          ///< Pauli type: 1 (X), 2 (Y), or 3 (Z).
      size_t position;   ///< Noise location index.
 };

 /**
  * @brief Compile a STIM program, build the QEPG, and generate weighted random samples.
  *
  * End-to-end convenience function: parses the STIM string, constructs the
  * backward propagation graph, generates `shots` random error configurations
  * of the given weight, and returns detector + observable outcomes.
  *
  * @param prog_str STIM program string (normalized, one instruction per line).
  * @param weight   Hamming weight of each error sample.
  * @param shots    Number of samples to generate.
  * @return 2D boolean vector of shape (shots, num_detectors + 1). The last
  *         column of each row is the observable outcome.
  */
 std::vector<std::vector<bool>> return_samples(const std::string & prog_str,size_t weight, size_t shots);

 /**
  * @brief Generate weighted random samples using a pre-compiled QEPG graph.
  *
  * Avoids recompilation when the same circuit is sampled multiple times with
  * different weights or shot counts.
  *
  * @param graph  Pre-compiled QEPG graph (from compile_QEPG).
  * @param weight Hamming weight of each error sample.
  * @param shots  Number of samples to generate.
  * @return 2D boolean vector of shape (shots, num_detectors + 1).
  */
 std::vector<std::vector<bool>> return_samples_with_fixed_QEPG(const QEPG::QEPG& graph,size_t weight, size_t shots);

 /**
  * @brief Generate weighted samples with separate NumPy detector/observable outputs.
  *
  * Same as return_samples_with_fixed_QEPG but returns pre-allocated NumPy arrays
  * with OpenMP-parallel unpacking for maximum throughput.
  *
  * @param graph  Pre-compiled QEPG graph.
  * @param weight Hamming weight of each error sample.
  * @param shots  Number of samples to generate.
  * @return A pair of NumPy arrays: (detector_outcomes of shape (shots, num_detectors),
  *         observable_outcomes of shape (shots,)).
  */
 std::pair<py::array_t<std::uint8_t>,py::array_t<std::uint8_t>> return_samples_with_fixed_QEPG_numpy(const QEPG::QEPG& graph,size_t weight, size_t shots);

 /**
  * @brief Generate weighted samples and return both noise vectors and outcomes.
  *
  * In addition to the detector/observable outcomes, also returns the noise
  * configuration (list of (qubit_index, pauli_type) pairs) for each sample.
  * Useful for debugging or analyzing which errors caused which syndromes.
  *
  * @param prog_str STIM program string.
  * @param weight   Hamming weight of each error sample.
  * @param shots    Number of samples to generate.
  * @return A pair: first element is a vector of noise vectors (each noise is a
  *         pair of (position, pauli_type)); second element is the 2D boolean
  *         outcome matrix.
  */
 std::pair<std::vector<std::vector<std::pair<int,int>>> ,std::vector<std::vector<bool>>>  return_samples_with_noise_vector(const std::string & prog_str,size_t weight, size_t shots);


 /**
  * @brief Generate samples for multiple weights in a single call.
  *
  * Compiles the circuit once and samples at each (weight, shot_count) pair.
  *
  * @param prog_str STIM program string.
  * @param weight   Vector of weights to sample at.
  * @param shots    Vector of shot counts, one per weight (must be same length as weight).
  * @return 3D boolean vector: result[i] is the 2D outcome matrix for weight[i].
  */
 std::vector<std::vector<std::vector<bool>>> return_samples_many_weights(const std::string& prog_str,const std::vector<size_t>& weight, const std::vector<size_t>& shots);


 /**
  * @brief Generate samples for multiple weights, returning NumPy arrays.
  *
  * Same as return_samples_many_weights but returns pybind11 NumPy arrays
  * instead of nested vectors for zero-copy Python interop.
  *
  * @param prog_str STIM program string.
  * @param weight   Vector of weights to sample at.
  * @param shots    Vector of shot counts per weight.
  * @return Vector of NumPy bool arrays, one per weight.
  */
 std::vector<py::array_t<bool>> return_samples_many_weights_numpy(const std::string& prog_str,const std::vector<size_t>& weight, const std::vector<size_t>& shots);

 /**
  * @brief Generate samples for multiple weights with separate detector and observable arrays.
  *
  * Compiles the circuit and concatenates all samples into two contiguous NumPy arrays:
  * one for detector outcomes and one for observable outcomes.
  *
  * @param prog_str STIM program string.
  * @param weight   Vector of weights.
  * @param shots    Vector of shot counts per weight.
  * @return A pair of NumPy arrays: (detector_outcomes of shape (total_shots, num_detectors),
  *         observable_outcomes of shape (total_shots,)).
  */
 std::pair<py::array_t<bool>,py::array_t<bool>> return_samples_many_weights_separate_obs(const std::string& prog_str,const std::vector<size_t>& weight, const std::vector<size_t>& shots);


 /**
  * @brief Generate samples for multiple weights with separate outputs, using a pre-compiled QEPG.
  *
  * @param graph  Pre-compiled QEPG graph.
  * @param weight Vector of weights.
  * @param shots  Vector of shot counts per weight.
  * @return A pair of NumPy arrays: (detector_outcomes, observable_outcomes).
  */
 std::pair<py::array_t<bool>,py::array_t<bool>> return_samples_many_weights_separate_obs_with_QEPG(const QEPG::QEPG& graph,const std::vector<size_t>& weight, const std::vector<size_t>& shots);


 /**
  * @brief Generate Monte Carlo (Bernoulli) error samples with separate detector and observable outputs.
  *
  * Each noise location independently has an error with the given probability.
  * Returns detector and observable outcomes as separate NumPy arrays.
  *
  * @param graph      Pre-compiled QEPG graph.
  * @param error_rate Per-location error probability.
  * @param shot       Number of samples to generate.
  * @return A pair of NumPy arrays: (detector_outcomes of shape (shot, num_detectors),
  *         observable_outcomes of shape (shot,)).
  */
 std::pair<py::array_t<std::uint8_t>,py::array_t<std::uint8_t>> return_samples_Monte_separate_obs_with_QEPG(const QEPG::QEPG& graph,const double& error_rate, const size_t& shot);


 /**
  * @brief Enumerate all error configurations of a given weight and return their outcomes.
  *
  * @param prog_str STIM program string.
  * @param weight   Hamming weight for exhaustive enumeration.
  * @return 2D boolean vector of all distinct outcomes.
  *
  * @note Currently delegates to a stub; not fully implemented.
  */
 std::vector<std::vector<bool>> return_all_samples_with_fixed_weights(const std::string& prog_str,const size_t& weight);


 /**
  * @brief Extract the parity propagation matrix (transposed) as a 2D boolean vector.
  *
  * Compiles the circuit, builds the QEPG, and returns the transposed parity
  * propagation matrix. Also prints the circuit and detector matrix to stdout.
  *
  * @param prog_str STIM program string.
  * @return 2D boolean vector representing the transposed parity propagation matrix.
  */
 std::vector<std::vector<bool>> return_detector_matrix(const std::string& prog_str);


 /**
  * @brief Generate samples and return as a NumPy array.
  *
  * Convenience function that compiles the circuit, builds the QEPG, samples,
  * and returns results directly as a NumPy bool array of shape (shots, num_detectors + 1).
  *
  * @param prog_str STIM program string.
  * @param weight   Hamming weight of each error sample.
  * @param shots    Number of samples.
  * @return NumPy bool array of shape (shots, num_detectors + 1).
  */
 py::array_t<bool> return_samples_numpy(const std::string& prog_str,size_t weight, size_t shots);


 /**
  * @brief Compile a STIM program into a reusable QEPG graph object.
  *
  * Parses the STIM string, builds the Clifford circuit, constructs the QEPG
  * via backward_graph_construction(), and returns the result. The returned
  * graph can be passed to functions like return_samples_with_fixed_QEPG()
  * or return_samples_many_weights_separate_obs_with_QEPG() for repeated
  * sampling without recompilation.
  *
  * @param prog_str STIM program string (normalized).
  * @return The compiled QEPG graph, ready for sampling.
  */
 QEPG::QEPG compile_QEPG(const std::string& prog_str);


 /**
  * @brief Generate non-uniform noise samples with separate detector and observable outputs.
  *
  * Each noise source has independent (px, py, pz) probabilities. Supports
  * DEPOLARIZE2 correlated pairs. Uses SIMD XOR, OpenMP parallelism, and
  * Xoshiro256pp RNG for maximum throughput.
  *
  * @param graph           Pre-compiled QEPG graph.
  * @param noise_probs     NumPy array of shape (N, 3) with per-source [px, py, pz].
  * @param corr_sources_a  Source A indices for correlated pairs.
  * @param corr_sources_b  Source B indices for correlated pairs.
  * @param corr_probs      Probabilities for correlated pairs.
  * @param shot            Number of samples to generate.
  * @return A pair of NumPy uint8 arrays: (detector_outcomes, observable_outcomes).
  */
 std::pair<py::array_t<std::uint8_t>, py::array_t<std::uint8_t>>
 return_samples_nonuniform_to_numpy(
     const QEPG::QEPG& graph,
     py::array_t<double> noise_probs,
     py::array_t<std::size_t> corr_sources_a,
     py::array_t<std::size_t> corr_sources_b,
     py::array_t<double> corr_probs,
     std::size_t shot);


}


#endif
