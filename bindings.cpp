// bindings.cpp
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>      // For std::vector, std::string
#include <pybind11/operators.h> // For exposing operators like ==

// Include headers for your C++ types
// Adjust paths if your .hpp files are in a subdirectory (e.g., QEPG/clifford.hpp)
#include "QEPG/clifford.hpp"
#include "QEPG/QEPG.hpp"
#include "QEPG/sampler.hpp"
#include "QEPG/LERcalculator.hpp"
// Include stimparser.hpp if needed by the bindings or types
// #include "QEPG/stimparser.hpp"

// Include Boost headers for types you are binding (like dynamic_bitset)
#include <boost/dynamic_bitset.hpp>

namespace py = pybind11;

// Assuming return_samples is in the SAMPLE namespace as per sampler.hpp
// Declare the function you want to bind if it's not defined in this file
namespace SAMPLE {
    class sampler; // Forward declare if binding methods/ctors
}
// Declare other classes if needed for binding their members
namespace clifford {
    class cliffordcircuit;
}
namespace QEPG {
    class QEPG;
    // using Row = boost::dynamic_bitset<>; // Already in QEPG.hpp
}


namespace LERcalculator{
    std::vector<QEPG::Row> return_samples(std::string prog_str, size_t weight, size_t shots);
}
   

// --- Bindings ---
PYBIND11_MODULE(QEPG, m) { // Use the module name 'QEPG' as seen in your build output
    m.doc() = "Pybind11 bindings for QEPG library";

    // 1. Bind boost::dynamic_bitset (your QEPG::Row)
    // Expose methods to interact with it from Python
    py::class_<boost::dynamic_bitset<>>(m, "DynamicBitset") // Python name for the bitset type
        // You might want a constructor from size
        // .def(py::init<size_t>(), py::arg("size"))
        // Or from a string/list (requires conversion logic)

        .def("size", &boost::dynamic_bitset<>::size, "Get the size of the bitset")
        .def("test", &boost::dynamic_bitset<>::test, py::arg("pos"), "Test if the bit at position pos is set")
        // Add operators if needed (==, !=, etc.)
        .def(py::self == py::self)

        // Add conversion to Python list for easier access
        .def("to_list", [](const boost::dynamic_bitset<>& self) {
             std::vector<bool> list(self.size());
             for(size_t i=0; i<self.size(); ++i) list[i] = self[i];
             return list;
         }, "Convert the bitset to a list of booleans")

        // Add a useful string representation
        .def("__repr__", [](const boost::dynamic_bitset<>& self) {
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

    // 2. Bind clifford::cliffordcircuit class
    py::class_<clifford::cliffordcircuit>(m, "CliffordCircuit")
        .def(py::init<>()) // Default constructor
        // Add other constructors if they exist and you need to use them from Python
        // .def(py::init<size_t>(), py::arg("num_qubits"))
        .def("compile_from_rewrited_stim_string", &clifford::cliffordcircuit::compile_from_rewrited_stim_string,
             py::arg("prog_str"), "Compile circuit from Stim string")
        // Expose getters used in return_samples or potentially useful in Python
        .def("get_num_detector", &clifford::cliffordcircuit::get_num_detector)
        .def("get_num_noise", &clifford::cliffordcircuit::get_num_noise)
        .def("get_num_qubit", &clifford::cliffordcircuit::get_num_qubit)
        // Bind other methods of CliffordCircuit if needed
        ;

    // 3. Bind QEPG::QEPG class
    py::class_<QEPG::QEPG>(m, "QEPGGraph")
        // Bind constructors, making sure argument types match bound C++ types
        .def(py::init<clifford::cliffordcircuit, size_t, size_t>(), // Takes a CliffordCircuit object
             py::arg("circuit"), py::arg("num_detector"), py::arg("num_noise"))
        .def("backward_graph_construction", &QEPG::QEPG::backward_graph_construction)
        // Bind getters for matrices if you want to access them from Python
        // .def("get_detectorMatrix", &QEPG::QEPG::get_detectorMatrix)
        // .def("get_parityPropMatrix", &QEPG::QEPG::get_parityPropMatrix)
        ;

    // 4. Bind SAMPLE::sampler class
    py::class_<SAMPLE::sampler>(m, "Sampler")
        .def(py::init<size_t>(), py::arg("num_total_paulierror"))
        // Bind generate_many_output_samples if you want to call it directly from Python
        // Note: return_samples will call it internally, so not strictly necessary to bind it here
        // .def("generate_many_output_samples", &SAMPLE::sampler::generate_many_output_samples, ...)
        ;

    // 5. Bind the top-level function return_samples
    // It's in the SAMPLE namespace
    m.def("return_samples", &LERcalculator::return_samples, // Use &SAMPLE::return_samples
          py::arg("prog_str"), py::arg("weight"), py::arg("shots"),
          "Function that returns samples based on a circuit and parameters");

}