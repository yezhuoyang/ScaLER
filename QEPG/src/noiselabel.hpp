/**
 * @file noiselabel.hpp
 * @brief C++ noise label classifier for QStab IR error types.
 *
 * Positional classifier that labels each DEPOLARIZE1 noise source into
 * one of four QStab IR error types:
 *   - Type 0 (data_qubit_error):    data qubit error before/after CX phase
 *   - Type I (ghost_error):         data qubit error during CX phase with future CX
 *   - Type II (hook_error):         ancilla error with remaining CX
 *   - Type III (measurement_error): ancilla error after last CX
 *
 * Classification is POSITIONAL — depends on where the noise occurs
 * relative to the CX schedule, not gate adjacency.
 */

#ifndef NOISELABEL_HPP
#define NOISELABEL_HPP
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>
#include <string>
#include <algorithm>
#include "clifford.hpp"

namespace noiselabel {

/// The four canonical QStab IR error types.
enum class NoiseType : std::uint8_t {
    DATA_QUBIT_ERROR    = 0,  // Type 0 (λ₀)
    GHOST_ERROR         = 1,  // Type I (λ₁)
    HOOK_ERROR          = 2,  // Type II (λ₂)
    MEASUREMENT_ERROR   = 3,  // Type III (λ₃)
    NUM_TYPES           = 4
};

/// Convert NoiseType to string name.
inline const char* noise_type_name(NoiseType t) noexcept {
    switch (t) {
        case NoiseType::DATA_QUBIT_ERROR:  return "data_qubit_error";
        case NoiseType::GHOST_ERROR:       return "ghost_error";
        case NoiseType::HOOK_ERROR:        return "hook_error";
        case NoiseType::MEASUREMENT_ERROR: return "measurement_error";
        default:                           return "unknown";
    }
}

/**
 * @brief Per-noise-source label entry.
 */
struct NoiseLabel {
    std::size_t noise_index;   ///< Index into the noise source array.
    std::size_t qubit;         ///< Qubit this noise acts on.
    NoiseType   type;          ///< Classified error type.
    bool        is_data;       ///< True if noise is on a data qubit.
};


/**
 * @brief Compact noise label map: array of labels indexed by noise source.
 *
 * After classification, `labels[i]` gives the NoiseType for noise source i.
 * Also provides category counts and index lookups by type.
 */
class NoiseLabelMap {
public:
    NoiseLabelMap() = default;

    /// Construct with a given number of noise sources.
    explicit NoiseLabelMap(std::size_t num_noise)
        : labels_(num_noise), num_noise_(num_noise) {}

    /// Set the label for a specific noise index.
    void set(std::size_t noise_index, NoiseType type, std::size_t qubit, bool is_data) {
        if (noise_index >= num_noise_) return;
        labels_[noise_index] = {noise_index, qubit, type, is_data};
    }

    /// Get the type of a noise source.
    NoiseType get_type(std::size_t noise_index) const {
        return labels_[noise_index].type;
    }

    /// Get the full label entry.
    const NoiseLabel& get(std::size_t noise_index) const {
        return labels_[noise_index];
    }

    /// Total number of noise sources.
    std::size_t size() const noexcept { return num_noise_; }

    /// Count sources of each type.
    std::array<std::size_t, 4> category_counts() const {
        std::array<std::size_t, 4> counts = {0, 0, 0, 0};
        for (std::size_t i = 0; i < num_noise_; ++i) {
            counts[static_cast<std::uint8_t>(labels_[i].type)]++;
        }
        return counts;
    }

    /// Get indices of all sources with a given type.
    std::vector<std::size_t> get_indices(NoiseType type) const {
        std::vector<std::size_t> result;
        for (std::size_t i = 0; i < num_noise_; ++i) {
            if (labels_[i].type == type) {
                result.push_back(i);
            }
        }
        return result;
    }

    /// Access underlying label array.
    const std::vector<NoiseLabel>& labels() const noexcept { return labels_; }

private:
    std::vector<NoiseLabel> labels_;
    std::size_t num_noise_ = 0;
};


/**
 * @brief Classify all noise sources in a compiled cliffordcircuit.
 *
 * Uses positional classification:
 *   - Pre-scans for CX(data↔ancilla) positions and ancilla measurement positions
 *   - Data qubit noise before first CX → Type 0
 *   - Data qubit noise during CX phase with future CX → Type I (ghost)
 *   - Data qubit noise during CX phase without future CX → Type 0
 *   - Data qubit noise after all ancilla measurements → Type 0
 *   - Ancilla noise with future CX → Type II (hook)
 *   - Ancilla noise without future CX → Type III (measurement)
 *
 * @param circuit         Compiled cliffordcircuit (gate list with DEPOLARIZE1 entries).
 * @param num_data_qubits Number of data qubits (indices 0..n-1).
 * @return A populated NoiseLabelMap.
 */
inline NoiseLabelMap auto_label(
    const clifford::cliffordcircuit& circuit,
    std::size_t num_data_qubits)
{
    const std::size_t num_noise = circuit.get_num_noise();
    const std::size_t num_gates = circuit.get_gate_num();

    NoiseLabelMap label_map(num_noise);

    // ------------------------------------------------------------------
    // Pre-scan: CX positions per qubit, first_cx_pos, last_ancilla_m_pos
    // ------------------------------------------------------------------
    // CX positions per qubit (only data↔ancilla CX gates)
    // Use flat vectors for cache efficiency
    std::vector<std::vector<std::size_t>> cx_positions(circuit.get_num_qubit());
    std::size_t first_cx_pos = num_gates;  // sentinel
    std::size_t last_ancilla_m_pos = 0;
    bool any_ancilla_m = false;

    for (std::size_t i = 0; i < num_gates; ++i) {
        const auto& gate = circuit.get_gate(i);
        if (gate.name == "cnot" && gate.qubits.size() == 2) {
            bool ctrl_data = gate.qubits[0] < num_data_qubits;
            bool tgt_data  = gate.qubits[1] < num_data_qubits;
            if (ctrl_data != tgt_data) {  // mixed data↔ancilla
                cx_positions[gate.qubits[0]].push_back(i);
                cx_positions[gate.qubits[1]].push_back(i);
                if (i < first_cx_pos) first_cx_pos = i;
            }
        } else if (gate.name == "M" && gate.qubits.size() >= 1
                   && gate.qubits[0] >= num_data_qubits) {
            if (i > last_ancilla_m_pos || !any_ancilla_m) {
                last_ancilla_m_pos = i;
                any_ancilla_m = true;
            }
        }
    }

    // ------------------------------------------------------------------
    // Classify every DEPOLARIZE1 by position
    // ------------------------------------------------------------------
    std::size_t noise_counter = 0;
    for (std::size_t i = 0; i < num_gates; ++i) {
        const auto& gate = circuit.get_gate(i);
        if (gate.name != "DEPOLARIZE1") continue;

        std::size_t qubit = gate.qubits[0];
        bool is_data = qubit < num_data_qubits;

        // Check for future CX for this qubit
        const auto& qcx = cx_positions[qubit];
        bool has_future_cx = std::any_of(qcx.begin(), qcx.end(),
            [i](std::size_t cx_pos) { return cx_pos > i; });

        NoiseType type;
        if (is_data) {
            // Data qubit classification
            if (i < first_cx_pos) {
                type = NoiseType::DATA_QUBIT_ERROR;  // Before any CX
            } else if (any_ancilla_m && i > last_ancilla_m_pos) {
                type = NoiseType::DATA_QUBIT_ERROR;  // After all ancilla M
            } else if (has_future_cx) {
                type = NoiseType::GHOST_ERROR;        // During phase, future CX
            } else {
                type = NoiseType::DATA_QUBIT_ERROR;   // During phase, no future CX
            }
        } else {
            // Ancilla qubit classification
            if (has_future_cx) {
                type = NoiseType::HOOK_ERROR;          // Back-propagation
            } else {
                type = NoiseType::MEASUREMENT_ERROR;   // Only flips M
            }
        }

        label_map.set(noise_counter, type, qubit, is_data);
        ++noise_counter;
    }

    return label_map;
}

}  // namespace noiselabel

#endif  // NOISELABEL_HPP
