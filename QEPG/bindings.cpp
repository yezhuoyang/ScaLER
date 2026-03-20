/**
 * @file bindings.cpp
 * @brief pybind11 module definition exposing the QEPG C++ backend to Python.
 *
 * Defines the `qepg` Python module with bindings for:
 * - DynamicBitset: wrapper around boost::dynamic_bitset for GF(2) row vectors
 * - CliffordCircuit: STIM program compiler
 * - QEPGGraph: backward error propagation graph
 * - Sampler: random Pauli error generator
 * - Module-level functions for sampling (return_samples, compile_QEPG, etc.)
 *
 * All sampling functions use move semantics (py::return_value_policy::move)
 * to avoid unnecessary copies when returning large arrays to Python.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

#include "src/clifford.hpp"
#include "src/QEPG.hpp"
#include "src/sampler.hpp"
#include "src/LERcalculator.hpp"
#include "src/noiselabel.hpp"
#include "src/hotspot.hpp"

namespace py = pybind11;


namespace SAMPLE {
    class sampler; // Forward declare if binding methods/ctors
}
// Declare other classes if needed for binding their members
namespace clifford {
    class cliffordcircuit;
}
namespace QEPG {
    class QEPG;
}


namespace LERcalculator{
    std::vector<std::vector<bool>> return_samples(const std::string& prog_str, size_t weight, size_t shots);
    std::vector<std::vector<std::vector<bool>>> return_samples_many_weights(const std::string& prog_str,const std::vector<size_t>& weight, const std::vector<size_t>& shots);
    std::vector<std::vector<bool>> return_detector_matrix(const std::string& prog_str);
    std::pair<std::vector<std::vector<std::pair<int,int>>> ,std::vector<std::vector<bool>>>  return_samples_with_noise_vector(const std::string & prog_str,size_t weight, size_t shots);
    std::pair<py::array_t<bool>,py::array_t<bool>> return_samples_many_weights_separate_obs(const std::string& prog_str,const std::vector<size_t>& weight, const std::vector<size_t>& shots);
    py::array_t<bool> return_samples_numpy(const std::string& prog_str,size_t weight, size_t shots);
    std::vector<py::array_t<bool>> return_samples_many_weights_numpy(const std::string& prog_str,const std::vector<size_t>& weight, const std::vector<size_t>& shots);
    QEPG::QEPG compile_QEPG(const std::string& prog_str);
    std::pair<py::array_t<bool>,py::array_t<bool>> return_samples_many_weights_separate_obs_with_QEPG(const QEPG::QEPG& graph,const std::vector<size_t>& weight, const std::vector<size_t>& shots);
    std::vector<std::vector<bool>> return_samples_with_fixed_QEPG(const QEPG::QEPG& graph,size_t weight, size_t shots);
    std::pair<py::array_t<bool>,py::array_t<bool>> return_samples_Monte_separate_obs_with_QEPG(const QEPG::QEPG& graph,const double& error_rate, const size_t& shot);
    std::pair<py::array_t<std::uint8_t>, py::array_t<std::uint8_t>>
    return_samples_nonuniform_to_numpy(
        const QEPG::QEPG& graph,
        py::array_t<double> noise_probs,
        py::array_t<std::size_t> corr_sources_a,
        py::array_t<std::size_t> corr_sources_b,
        py::array_t<double> corr_probs,
        std::size_t shot);
}


/// @brief Define the qepg Python module and all its bindings.
PYBIND11_MODULE(qepg, m) {
    m.doc() = R"pbdoc(
        QEPG: Quantum Error Propagation Graph library.

        A C++ backend for ScaLERQEC that compiles STIM-style quantum error
        correction circuits, builds backward error propagation graphs over
        GF(2), and generates weighted random Pauli error samples with their
        detector/observable outcomes. Uses a custom DynamicBitset (uint64_t
        backed) for efficient bit-level operations and OpenMP for parallel
        sampling.
    )pbdoc";

    /**
     * @brief Python binding for qepg_bits::DynamicBitset.
     *
     * Exposes size(), test(), equality comparison, list conversion,
     * and a truncated string representation for large bitsets.
     */
    py::class_<QEPG::Row>(m, "DynamicBitset",
        "A dynamic-length bitset for GF(2) row vectors.")

        .def("size", &QEPG::Row::size,
             "Return the number of bits in the bitset.")
        .def("test", &QEPG::Row::test, py::arg("pos"),
             "Return True if the bit at the given position is set.")
        .def(py::self == py::self)

        .def("to_list", [](const QEPG::Row& self) {
             std::vector<bool> list(self.size());
             for(size_t i=0; i<self.size(); ++i) list[i] = self[i];
             return list;
         }, "Convert the bitset to a Python list of booleans.")

        .def("__repr__", [](const QEPG::Row& self) {
            std::string s = "<DynamicBitset ";
            if (self.size() > 40) { // Truncate long outputs
                 for(size_t i=0; i<20; ++i) s += (i < self.size() ? (self[i] ? '1' : '0') : '-');
                 s += "...";
                 for(size_t i=self.size()-20; i<self.size(); ++i) s += (self[i] ? '1' : '0');
             } else { // Show full bitstring for small bitsets
                 for(size_t i=0; i<self.size(); ++i) s += (self[i] ? '1' : '0');
             }
            s += ">";
            return s;
        })
        ;

    /**
     * @brief Python binding for clifford::cliffordcircuit.
     *
     * Allows compiling STIM strings and querying circuit metadata from Python.
     */
    py::class_<clifford::cliffordcircuit>(m, "CliffordCircuit",
        "Compiles a normalized STIM program string into an internal Clifford circuit representation.")
        .def(py::init<>(), "Create an empty CliffordCircuit.")
        .def("compile_from_rewrited_stim_string", &clifford::cliffordcircuit::compile_from_rewrited_stim_string,
             py::arg("prog_str"),
             "Parse a normalized STIM program string and populate the circuit.")
        .def("get_num_detector", &clifford::cliffordcircuit::get_num_detector,
             "Return the number of detectors defined in the circuit.")
        .def("get_num_noise", &clifford::cliffordcircuit::get_num_noise,
             "Return the number of depolarizing noise locations.")
        .def("get_num_qubit", &clifford::cliffordcircuit::get_num_qubit,
             "Return the number of qubits in the circuit.")
        ;

    /**
     * @brief Python binding for QEPG::QEPG.
     *
     * Wraps the backward error propagation graph. Construct from a
     * CliffordCircuit, then call backward_graph_construction().
     */
    py::class_<QEPG::QEPG>(m, "QEPGGraph",
        "Quantum Error Propagation Graph. Maps Pauli noise locations to detector/observable outcomes.")
        .def(py::init<clifford::cliffordcircuit, size_t, size_t>(),
             py::arg("circuit"), py::arg("num_detector"), py::arg("num_noise"),
             "Construct a QEPG from a compiled CliffordCircuit.")
        .def("backward_graph_construction", &QEPG::QEPG::backward_graph_construction,
             "Build the parity propagation matrix by backward circuit traversal.")
        .def("get_total_noise", &QEPG::QEPG::get_total_noise,
             "Return the total number of noise locations.")
        .def("get_total_detector", &QEPG::QEPG::get_total_detector,
             "Return the total number of detectors.")
        ;

    /**
     * @brief Python binding for SAMPLE::sampler.
     *
     * Exposes the sampler constructor for use in custom sampling workflows.
     */
    py::class_<SAMPLE::sampler>(m, "Sampler",
        "Random Pauli error sampler with fixed-weight and Monte Carlo modes.")
        .def(py::init<size_t>(), py::arg("num_total_paulierror"),
             "Create a sampler for the given number of noise locations.")
        ;

    // --- Module-level sampling functions ---

    m.def("return_samples", &LERcalculator::return_samples,
          py::arg("prog_str"), py::arg("weight"), py::arg("shots"),
          py::return_value_policy::move,
          R"pbdoc(
              Compile a STIM program, build the QEPG, and generate random error samples.

              Args:
                  prog_str: Normalized STIM program string.
                  weight: Hamming weight of each error sample.
                  shots: Number of samples to generate.

              Returns:
                  2D list of bools, shape (shots, num_detectors + 1).
          )pbdoc");

    m.def("return_samples_with_fixed_QEPG", &LERcalculator::return_samples_with_fixed_QEPG,
          py::arg("graph"), py::arg("weight"), py::arg("shots"),
          py::return_value_policy::move,
          R"pbdoc(
              Generate random error samples using a pre-compiled QEPG graph.

              Args:
                  graph: Pre-compiled QEPGGraph object.
                  weight: Hamming weight of each error sample.
                  shots: Number of samples to generate.

              Returns:
                  2D list of bools, shape (shots, num_detectors + 1).
          )pbdoc");



    m.def("return_samples_Monte_separate_obs_with_QEPG",&LERcalculator::return_samples_Monte_separate_obs_with_QEPG,
        py::arg("graph"), py::arg("error_rate"),py::arg("shot"),
        py::return_value_policy::move,
        R"pbdoc(
            Generate Monte Carlo error samples with separate detector and observable outputs.

            Args:
                graph: Pre-compiled QEPGGraph object.
                error_rate: Per-location error probability.
                shot: Number of samples to generate.

            Returns:
                Tuple of (detector_outcomes, observable_outcomes) as NumPy bool arrays.
        )pbdoc");



    m.def("return_samples_numpy", &LERcalculator::return_samples_numpy,
          py::arg("prog_str"), py::arg("weight"), py::arg("shots"),
          py::return_value_policy::move,
          R"pbdoc(
              Compile, build QEPG, sample, and return results as a NumPy array.

              Args:
                  prog_str: Normalized STIM program string.
                  weight: Hamming weight of each error sample.
                  shots: Number of samples.

              Returns:
                  NumPy bool array of shape (shots, num_detectors + 1).
          )pbdoc");

    m.def("return_samples_many_weights", &LERcalculator::return_samples_many_weights,
        py::arg("prog_str"), py::arg("weight"), py::arg("shots"),
        py::return_value_policy::move,
        R"pbdoc(
            Generate samples at multiple weights in a single compilation pass.

            Args:
                prog_str: Normalized STIM program string.
                weight: List of weights to sample at.
                shots: List of shot counts, one per weight.

            Returns:
                3D list of bools: result[i] is the outcome matrix for weight[i].
        )pbdoc");


    m.def("return_samples_many_weights_numpy", &LERcalculator::return_samples_many_weights_numpy,
        py::arg("prog_str"), py::arg("weight"), py::arg("shots"),
        py::return_value_policy::move,
        R"pbdoc(
            Generate samples at multiple weights, returning NumPy arrays.

            Args:
                prog_str: Normalized STIM program string.
                weight: List of weights.
                shots: List of shot counts per weight.

            Returns:
                List of NumPy bool arrays, one per weight.
        )pbdoc");


    m.def("return_detector_matrix", &LERcalculator::return_detector_matrix,
          py::arg("prog_str"),
          R"pbdoc(
              Extract the transposed parity propagation matrix as a 2D bool list.

              Also prints the circuit and detector matrix to stdout for debugging.

              Args:
                  prog_str: Normalized STIM program string.

              Returns:
                  2D list of bools representing the transposed parity propagation matrix.
          )pbdoc");

    m.def("return_samples_with_noise_vector",
        &LERcalculator::return_samples_with_noise_vector,
        py::arg("prog_str"), py::arg("weight"), py::arg("shots"),
        py::return_value_policy::move,
        R"pbdoc(
            Generate samples and return both noise vectors and outcomes.

            Args:
                prog_str: Normalized STIM program string.
                weight: Hamming weight of each error sample.
                shots: Number of samples.

            Returns:
                Tuple of (noise_vectors, outcomes). Each noise vector is a list
                of (position, pauli_type) pairs.
        )pbdoc");

    m.def("return_samples_many_weights_separate_obs",
        &LERcalculator::return_samples_many_weights_separate_obs,
        py::arg("prog_str"), py::arg("weight"), py::arg("shots"),
        py::return_value_policy::move,
        R"pbdoc(
            Generate samples at multiple weights with separate detector and observable arrays.

            Args:
                prog_str: Normalized STIM program string.
                weight: List of weights.
                shots: List of shot counts per weight.

            Returns:
                Tuple of (detector_outcomes, observable_outcomes) as NumPy bool arrays.
        )pbdoc");

    m.def(
        "compile_QEPG",
        &LERcalculator::compile_QEPG,
        py::arg("prog_str"),
        R"pbdoc(
            Compile a STIM program into a reusable QEPGGraph object.

            Parses the program string, builds the Clifford circuit, runs
            backward_graph_construction(), and returns the resulting graph.
            Use this to avoid recompilation when sampling the same circuit
            with different parameters.

            Args:
                prog_str: Normalized STIM program string.

            Returns:
                QEPGGraph: The compiled error propagation graph.
        )pbdoc",
        py::return_value_policy::move
    );


    m.def("return_samples_many_weights_separate_obs_with_QEPG",
        &LERcalculator::return_samples_many_weights_separate_obs_with_QEPG,
        py::arg("graph"), py::arg("weight"), py::arg("shots"),
        py::return_value_policy::move,
        R"pbdoc(
            Generate samples at multiple weights with a pre-compiled QEPG, returning separate arrays.

            Args:
                graph: Pre-compiled QEPGGraph object.
                weight: List of weights.
                shots: List of shot counts per weight.

            Returns:
                Tuple of (detector_outcomes, observable_outcomes) as NumPy bool arrays.
        )pbdoc");


    m.def("return_samples_nonuniform_to_numpy",
        &LERcalculator::return_samples_nonuniform_to_numpy,
        py::arg("graph"),
        py::arg("noise_probs"),
        py::arg("corr_sources_a"),
        py::arg("corr_sources_b"),
        py::arg("corr_probs"),
        py::arg("shot"),
        py::return_value_policy::move,
        R"pbdoc(
            Generate non-uniform noise samples with per-source (px, py, pz) probabilities.

            Supports DEPOLARIZE2 correlated pairs. Uses SIMD XOR accumulation,
            OpenMP parallelism, and Xoshiro256pp RNG.

            Args:
                graph: Pre-compiled QEPGGraph object.
                noise_probs: NumPy float64 array of shape (N, 3) with [px, py, pz] per source.
                corr_sources_a: NumPy array of source_a indices for correlated pairs.
                corr_sources_b: NumPy array of source_b indices for correlated pairs.
                corr_probs: NumPy array of probabilities for correlated pairs.
                shot: Number of samples to generate.

            Returns:
                Tuple of (detector_outcomes, observable_outcomes) as NumPy uint8 arrays.
        )pbdoc");


    // =========================================================================
    // Noise Label System — positional classifier + hotspot analysis
    // =========================================================================

    py::enum_<noiselabel::NoiseType>(m, "NoiseType",
        "QStab IR error types for noise classification.")
        .value("DATA_QUBIT_ERROR",  noiselabel::NoiseType::DATA_QUBIT_ERROR,
               "Type 0: data qubit error before/after CX phase")
        .value("GHOST_ERROR",       noiselabel::NoiseType::GHOST_ERROR,
               "Type I: data qubit error during CX phase with future CX")
        .value("HOOK_ERROR",        noiselabel::NoiseType::HOOK_ERROR,
               "Type II: ancilla error with remaining CX")
        .value("MEASUREMENT_ERROR", noiselabel::NoiseType::MEASUREMENT_ERROR,
               "Type III: ancilla error after last CX")
        ;

    py::class_<noiselabel::NoiseLabel>(m, "NoiseLabel",
        "Per-noise-source label entry.")
        .def_readonly("noise_index", &noiselabel::NoiseLabel::noise_index)
        .def_readonly("qubit",       &noiselabel::NoiseLabel::qubit)
        .def_readonly("type",        &noiselabel::NoiseLabel::type)
        .def_readonly("is_data",     &noiselabel::NoiseLabel::is_data)
        .def("__repr__", [](const noiselabel::NoiseLabel& self) {
            return std::string("<NoiseLabel idx=") + std::to_string(self.noise_index)
                + " q=" + std::to_string(self.qubit)
                + " type=" + noiselabel::noise_type_name(self.type)
                + (self.is_data ? " data" : " ancilla") + ">";
        })
        ;

    py::class_<noiselabel::NoiseLabelMap>(m, "NoiseLabelMap",
        "Compact noise label map: array of labels indexed by noise source.")
        .def(py::init<>())
        .def(py::init<std::size_t>(), py::arg("num_noise"))
        .def("get_type", &noiselabel::NoiseLabelMap::get_type, py::arg("noise_index"),
             "Get the NoiseType of a noise source.")
        .def("get", &noiselabel::NoiseLabelMap::get, py::arg("noise_index"),
             py::return_value_policy::reference_internal,
             "Get the full NoiseLabel entry.")
        .def("size", &noiselabel::NoiseLabelMap::size,
             "Total number of noise sources.")
        .def("category_counts", &noiselabel::NoiseLabelMap::category_counts,
             "Return [type0_count, type1_count, type2_count, type3_count].")
        .def("get_indices", &noiselabel::NoiseLabelMap::get_indices, py::arg("type"),
             "Get indices of all sources with a given type.")
        .def("__len__", &noiselabel::NoiseLabelMap::size)
        .def("__repr__", [](const noiselabel::NoiseLabelMap& self) {
            auto counts = self.category_counts();
            return std::string("<NoiseLabelMap n=") + std::to_string(self.size())
                + " type0=" + std::to_string(counts[0])
                + " typeI=" + std::to_string(counts[1])
                + " typeII=" + std::to_string(counts[2])
                + " typeIII=" + std::to_string(counts[3]) + ">";
        })
        ;

    m.def("auto_label",
        &noiselabel::auto_label,
        py::arg("circuit"), py::arg("num_data_qubits"),
        py::return_value_policy::move,
        R"pbdoc(
            Classify all noise sources in a compiled CliffordCircuit.

            Uses positional classification relative to the CX schedule:
            - Data qubit noise before first CX → Type 0 (data_qubit_error)
            - Data qubit noise during phase with future CX → Type I (ghost_error)
            - Data qubit noise during phase without future CX → Type 0
            - Ancilla noise with future CX → Type II (hook_error)
            - Ancilla noise without future CX → Type III (measurement_error)

            Args:
                circuit: A compiled CliffordCircuit.
                num_data_qubits: Number of data qubits (indices 0..n-1).

            Returns:
                NoiseLabelMap with per-source labels.
        )pbdoc");


    // --- Hotspot analysis result types ---

    py::class_<hotspot::CategoryReport>(m, "CategoryReport",
        "Conditional probability report for a single noise category.")
        .def_readonly("type",                &hotspot::CategoryReport::type)
        .def_readonly("name",                &hotspot::CategoryReport::name)
        .def_readonly("count",               &hotspot::CategoryReport::count)
        .def_readonly("p_fired",             &hotspot::CategoryReport::p_fired)
        .def_readonly("p_fired_given_error", &hotspot::CategoryReport::p_fired_given_error)
        .def_readonly("p_error_given_fired", &hotspot::CategoryReport::p_error_given_fired)
        .def_readonly("lift",                &hotspot::CategoryReport::lift)
        ;

    py::class_<hotspot::ConfigEntry>(m, "ConfigEntry",
        "Multi-error configuration entry.")
        .def_readonly("configuration", &hotspot::ConfigEntry::configuration)
        .def_readonly("count",         &hotspot::ConfigEntry::count)
        .def_readonly("fraction",      &hotspot::ConfigEntry::fraction)
        ;

    py::class_<hotspot::SourceReport>(m, "SourceReport",
        "Per-source ranking entry.")
        .def_readonly("noise_index",         &hotspot::SourceReport::noise_index)
        .def_readonly("qubit",               &hotspot::SourceReport::qubit)
        .def_readonly("type",                &hotspot::SourceReport::type)
        .def_readonly("p_fired",             &hotspot::SourceReport::p_fired)
        .def_readonly("p_fired_given_error", &hotspot::SourceReport::p_fired_given_error)
        .def_readonly("lift",                &hotspot::SourceReport::lift)
        ;

    py::class_<hotspot::HotspotResult>(m, "HotspotResult",
        "Complete hotspot analysis result.")
        .def_readonly("num_shots",          &hotspot::HotspotResult::num_shots)
        .def_readonly("num_logical_errors", &hotspot::HotspotResult::num_logical_errors)
        .def_readonly("num_noise",          &hotspot::HotspotResult::num_noise)
        .def_readonly("category_reports",   &hotspot::HotspotResult::category_reports)
        .def_readonly("config_reports",     &hotspot::HotspotResult::config_reports)
        .def_readonly("source_reports",     &hotspot::HotspotResult::source_reports)
        ;


    // -----------------------------------------------------------------
    // Phase 1: labeled_sample — sampling with category bitmask tracking
    // -----------------------------------------------------------------
    m.def("labeled_sample",
        [](const QEPG::QEPG& graph,
           const noiselabel::NoiseLabelMap& label_map,
           py::array_t<double> noise_probs_arr,
           py::array_t<std::size_t> corr_sources_a_arr,
           py::array_t<std::size_t> corr_sources_b_arr,
           py::array_t<double> corr_probs_arr,
           std::size_t num_shots)
        {
            // Extract noise_probs
            auto np_info = noise_probs_arr.request();
            if (np_info.ndim != 2 || np_info.shape[1] != 3)
                throw std::runtime_error("noise_probs must have shape (N, 3)");
            const std::size_t num_noise = static_cast<std::size_t>(np_info.shape[0]);
            const double* noise_probs = static_cast<const double*>(np_info.ptr);

            // Extract correlated pairs
            auto ca_info = corr_sources_a_arr.request();
            auto cb_info = corr_sources_b_arr.request();
            auto cp_info = corr_probs_arr.request();
            const std::size_t num_corr = static_cast<std::size_t>(ca_info.shape[0]);

            std::vector<SAMPLE::CorrelatedPair> corr_pairs(num_corr);
            if (num_corr > 0) {
                const std::size_t* ca = static_cast<const std::size_t*>(ca_info.ptr);
                const std::size_t* cb = static_cast<const std::size_t*>(cb_info.ptr);
                const double* cp = static_cast<const double*>(cp_info.ptr);
                for (std::size_t i = 0; i < num_corr; ++i) {
                    corr_pairs[i] = {ca[i], cb[i], cp[i]};
                }
            }

            // Allocate output arrays
            const std::size_t n_det = graph.get_total_detector();
            py::array_t<std::uint8_t> det_result({num_shots, n_det});
            py::array_t<std::uint8_t> obs_result(num_shots);
            py::array_t<std::uint8_t> bitmask_result(num_shots);

            auto det_info = det_result.request();
            auto obs_info = obs_result.request();
            auto bm_info  = bitmask_result.request();

            {
                py::gil_scoped_release release;
                hotspot::labeled_sample(
                    graph, label_map,
                    static_cast<std::uint8_t*>(det_info.ptr),
                    static_cast<std::uint8_t*>(obs_info.ptr),
                    static_cast<std::uint8_t*>(bm_info.ptr),
                    n_det, noise_probs, num_noise,
                    num_corr > 0 ? corr_pairs.data() : nullptr, num_corr,
                    num_shots);
            }

            return py::make_tuple(
                std::move(det_result),
                std::move(obs_result),
                std::move(bitmask_result));
        },
        py::arg("graph"),
        py::arg("label_map"),
        py::arg("noise_probs"),
        py::arg("corr_sources_a"),
        py::arg("corr_sources_b"),
        py::arg("corr_probs"),
        py::arg("num_shots"),
        py::return_value_policy::move,
        R"pbdoc(
            Phase 1: Sample non-uniform noise with per-shot category bitmask tracking.

            Same sampling as return_samples_nonuniform_to_numpy, but additionally
            returns a per-shot category bitmask (uint8) where bit i = 1 if any
            noise source of Type i fired.

            After calling this, decode detector_outcomes with your chosen decoder,
            then pass bitmasks + logical_error mask to hotspot_aggregate().

            Args:
                graph: Pre-compiled QEPGGraph object.
                label_map: NoiseLabelMap classifying each noise source.
                noise_probs: NumPy float64 array of shape (N, 3) with [px, py, pz].
                corr_sources_a: Source A indices for correlated pairs.
                corr_sources_b: Source B indices for correlated pairs.
                corr_probs: Probabilities for correlated pairs.
                num_shots: Number of Monte Carlo shots.

            Returns:
                Tuple of (detector_outcomes, observable_outcomes, category_bitmasks)
                as NumPy uint8 arrays.
        )pbdoc");


    // -----------------------------------------------------------------
    // Phase 3: hotspot_aggregate — decoder-agnostic aggregation
    // -----------------------------------------------------------------
    m.def("hotspot_aggregate",
        [](py::array_t<std::uint8_t> bitmask_arr,
           py::array_t<std::uint8_t> logical_error_arr,
           const noiselabel::NoiseLabelMap& label_map,
           std::size_t min_config_count) -> hotspot::HotspotResult
        {
            auto bm_info = bitmask_arr.request();
            auto le_info = logical_error_arr.request();
            const std::size_t num_shots = static_cast<std::size_t>(bm_info.shape[0]);

            if (le_info.shape[0] != bm_info.shape[0])
                throw std::runtime_error(
                    "bitmask and logical_error arrays must have same length");

            hotspot::HotspotResult result;
            {
                py::gil_scoped_release release;
                result = hotspot::hotspot_aggregate(
                    static_cast<const std::uint8_t*>(bm_info.ptr),
                    static_cast<const std::uint8_t*>(le_info.ptr),
                    num_shots, label_map, min_config_count);
            }
            return result;
        },
        py::arg("bitmasks"),
        py::arg("logical_errors"),
        py::arg("label_map"),
        py::arg("min_config_count") = 5,
        py::return_value_policy::move,
        R"pbdoc(
            Phase 3: Aggregate hotspot statistics from bitmasks + decoded logical errors.

            Pure aggregation — no sampling, no decoding. Takes per-shot category
            bitmasks (from labeled_sample) and a logical_error mask (from YOUR
            decoder), and computes conditional probabilities and configuration
            breakdown.

            Args:
                bitmasks: Per-shot category bitmasks from labeled_sample, shape (N,).
                logical_errors: Per-shot logical error flags, shape (N,). uint8.
                    Compute as: (observable_outcomes != decoder.decode_batch(detectors)).
                label_map: NoiseLabelMap for category count info.
                min_config_count: Minimum count for config entries (default 5).

            Returns:
                HotspotResult with category_reports and config_reports.
        )pbdoc");

}
