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

#include "src/clifford.hpp"
#include "src/QEPG.hpp"
#include "src/sampler.hpp"
#include "src/LERcalculator.hpp"

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

}
