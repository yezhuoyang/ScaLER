/**
 * @file clifford.hpp
 * @brief Clifford circuit representation for STIM-style quantum error correction programs.
 *
 * Defines the data structures and class for compiling a normalized STIM program string
 * into an internal gate list with associated noise channels, measurements, detectors,
 * and observables. This compiled representation is consumed by the QEPG backward
 * error propagation graph builder.
 */

#ifndef CLIFFORD_HPP
#define CLIFFORD_HPP
#pragma once

#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

namespace clifford{

/**
 * @brief A single circuit instruction (gate, noise, measurement, or reset).
 *
 * Each Gate stores a string name identifying the operation (e.g. "h", "cnot",
 * "DEPOLARIZE1", "M", "R") and a list of qubit indices it acts on.
 */
// TODO: [Review-P1-Perf] Gate uses heap-allocated std::string + std::vector for
// every gate instance. For large circuits (millions of gates), this causes massive
// overhead (2 heap allocs per gate, ~80 bytes each). Replace with:
//   enum class GateKind : uint8_t { Depolarize1, M, R, CNOT, H, P, X, Y, Z, ... };
//   struct Gate { GateKind kind; size_t qubits[2]; uint8_t num_qubits; };
// This reduces Gate from ~80 bytes + 2 heap allocs to ~24 bytes with zero allocs,
// and enables switch-based dispatch in backward_graph_construction().
struct Gate{
    std::string name;              ///< Operation name (e.g. "h", "cnot", "DEPOLARIZE1", "M", "R").
    std::vector<size_t> qubits;    ///< Qubit indices this gate operates on.
};


/**
 * @brief A set of measurement indices that define a detector or observable.
 *
 * A detector is a parity check over a subset of measurement outcomes. This struct
 * stores the list of measurement indices whose XOR gives the detector value.
 */
struct paritygroup{
    std::vector<size_t> indexlist;  ///< Measurement indices belonging to this parity group.
};

/**
 * @brief Inverse mapping from a measurement index to the detector/observable indices it participates in.
 *
 * For a given measurement, this stores all detector (and observable) indices
 * that include that measurement in their parity check.
 */
struct parityIndexgroup{
    std::vector<size_t> indexlist;  ///< Detector/observable indices that reference this measurement.
};


/**
 * @brief Compiles a STIM-style circuit description into an internal gate list.
 *
 * The cliffordcircuit class parses a normalized STIM program string and builds
 * an ordered list of gates, noise channels, measurements, detectors, and
 * observables. Each Clifford gate (H, CX, S, X, Y, Z) is automatically
 * preceded by a DEPOLARIZE1 noise channel on every qubit it touches.
 * Measurements likewise have an associated noise channel.
 *
 * After compilation, the circuit provides accessors for all metadata needed
 * by the QEPG backward propagation algorithm: number of qubits, noise terms,
 * measurements, detectors, and the mappings between them.
 *
 * @note The input STIM string must already be "rewritten" (normalized) so that
 *       each line contains exactly one instruction with space-separated arguments.
 */
class cliffordcircuit{

public:

    /// @brief Default constructor. Creates an empty circuit.
    cliffordcircuit();

    /**
     * @brief Construct a circuit with a pre-specified number of qubits.
     * @param n_qubit Initial number of qubits in the circuit.
     */
    explicit cliffordcircuit(size_t n_qubit);

    /**
     * @brief Set the physical error rate for noise channels.
     * @param p Error probability per noise location.
     */
    void set_error_rate(double p);


    /* single-qubit noise----------------------------------*/

    /**
     * @brief Append an X-error noise channel on a single qubit.
     * @param qindex Target qubit index.
     */
    void add_XError(size_t qindex);

    /**
     * @brief Append a Z-error noise channel on a single qubit.
     * @param qindex Target qubit index.
     */
    void add_ZError(size_t qindex);

    /**
     * @brief Append a single-qubit depolarizing noise channel and increment the noise counter.
     * @param qindex Target qubit index.
     */
    void add_depolarize1(size_t qindex);

    /*Clifford gates---------------------------------------*/

    /**
     * @brief Append a CNOT (CX) gate with depolarizing noise on both qubits.
     * @param qcontrol Control qubit index.
     * @param qtarget  Target qubit index.
     */
    void add_cnot(size_t qcontrol, size_t qtarget);

    /**
     * @brief Append a Hadamard gate with depolarizing noise.
     * @param qindex Target qubit index.
     */
    void add_hadamard(size_t qindex);

    /**
     * @brief Append a Phase (S) gate with depolarizing noise.
     * @param qindex Target qubit index.
     */
    void add_phase(size_t qindex);

    /**
     * @brief Append a Pauli-X gate with depolarizing noise.
     * @param qindex Target qubit index.
     */
    void add_pauliX(size_t qindex);

    /**
     * @brief Append a Pauli-Y gate with depolarizing noise.
     * @param qindex Target qubit index.
     */
    void add_pauliy(size_t qindex);

    /**
     * @brief Append a Pauli-Z gate with depolarizing noise.
     * @param qindex Target qubit index.
     */
    void add_pauliz(size_t qindex);

    /*Reset/Measurement---------------------------------------*/

    /**
     * @brief Append a qubit reset operation (no noise channel added).
     * @param qindex Target qubit index.
     */
    void add_reset(size_t qindex);

    /**
     * @brief Append a measurement with an associated depolarizing noise channel.
     *
     * Increments both the measurement and noise counters. Also creates an
     * empty parityIndexgroup entry for the new measurement in the
     * measurement-to-detector mapping.
     *
     * @param qindex Target qubit index to measure.
     */
    void add_measurement(size_t qindex);

    /*misc -------------------------------------------------*/

    /**
     * @brief Print the full circuit to stdout for debugging.
     *
     * Outputs qubit count, error rate, gate list, detector definitions,
     * measurement-to-detector mappings, observable definition, and
     * measurement-to-gate-index mappings.
     */
    void print_circuit() const;


    /*--Get gate by index-------------------------------------*/

    /**
     * @brief Retrieve a gate by its position in the circuit.
     * @param gateindex Zero-based index into the gate list.
     * @return Const reference to the Gate at that position.
     * @throws std::out_of_range if gateindex is out of bounds.
     */
    const Gate& get_gate(size_t gateindex) const;

    /*Setter/Getter of class members-------------------------------------------*/

    /// @brief Get the number of qubits in the circuit.
    /// @return Number of qubits.
    size_t get_num_qubit() const;

    /**
     * @brief Set the number of qubits explicitly.
     * @param num_qubit New qubit count.
     */
    void set_num_qubit(size_t num_qubit);

    /// @brief Get the total number of gates (including noise channels) in the circuit.
    /// @return Gate count.
    size_t get_gate_num() const;

    /// @brief Get the total number of measurements in the circuit.
    /// @return Measurement count.
    size_t get_num_meas() const;

    /// @brief Get the total number of depolarizing noise locations in the circuit.
    /// @return Noise channel count.
    size_t get_num_noise() const;

    /// @brief Get the number of detectors defined in the circuit.
    /// @return Detector count.
    size_t get_num_detector() const;

    /**
     * @brief Get all detector parity groups.
     * @return Const reference to the vector of paritygroup structs, one per detector.
     */
    const std::vector<paritygroup>& get_detector_parity_group() const;

    /**
     * @brief Get the observable parity group.
     * @return Const reference to the paritygroup defining the logical observable.
     */
    const paritygroup& get_observable_parity_group() const;

    /**
     * @brief Get the detector indices associated with a given measurement.
     * @param mindex Measurement index.
     * @return Const reference to the parityIndexgroup listing detector indices
     *         that include this measurement.
     */
    const parityIndexgroup& get_measure_to_parity_index(const size_t& mindex) const;

    /**
     * @brief Parse a normalized STIM program string and populate the circuit.
     *
     * Processes lines one by one, recognizing the instructions: M, R, H, CX,
     * DETECTOR, and OBSERVABLE_INCLUDE. For DETECTOR and OBSERVABLE lines, it
     * parses rec[...] tokens to build the measurement-to-detector and
     * detector-to-measurement mappings.
     *
     * @param stim_str The full STIM program as a single string (newline-delimited).
     *
     * @note The string must be pre-normalized: one instruction per line, no
     *       REPEAT blocks, no inline noise probabilities. The Python layer
     *       handles this normalization before calling into C++.
     */
    void compile_from_rewrited_stim_string(std::string stim_str);

private:

    size_t     num_qubit_=0;       ///< Number of qubits (auto-expanded as gates are added).
    size_t     num_noise_=0;       ///< Running count of DEPOLARIZE1 noise locations.
    size_t     num_meas_=0;        ///< Running count of measurements.
    size_t     num_detectors_=0;   ///< Running count of detectors.

    double  error_rate_{0.0};      ///< Physical error rate for noise channels.

    std::vector<Gate> circuit_;                        ///< Ordered list of all circuit instructions.
    std::vector<size_t> measureindexList_;              ///< Maps measurement ordinal to gate index in circuit_.

    std::vector<paritygroup> detectors_;                ///< Detector definitions: detector_index -> measurement indices.

    /**
     * @brief Inverse mapping: measurement_index -> list of detector indices.
     *
     * For each measurement, stores which detectors reference it. Built
     * during DETECTOR line parsing.
     */
    std::vector<parityIndexgroup> measure_to_parity_index_;

    paritygroup observable_;                            ///< Observable definition: measurement indices whose parity gives the logical outcome.

};

}

#endif
