/**
 * @file QEPG.hpp
 * @brief Quantum Error Propagation Graph (QEPG) over GF(2).
 *
 * Implements the backward error propagation algorithm that traces how each
 * Pauli noise location (X, Y, Z) propagates through a Clifford circuit to
 * affect detector and observable outcomes. The result is a binary matrix
 * (the parity propagation matrix) stored using Boost dynamic_bitset for
 * efficient GF(2) arithmetic.
 *
 * Uses a custom DynamicBitset class (backed by uint64_t words) instead of
 * Boost, eliminating all non-standard C++ dependencies.
 */

#ifndef QEPG_HPP
#define QEPG_HPP
#pragma once
#include <cstddef>
#include <ostream>
#include <string>
#include <bitset>
#include "dynamic_bitset.hpp"
#include "flat_bit_table.hpp"
#include "clifford.hpp"
#include <iostream>


namespace QEPG{

/// @brief A single row in a GF(2) matrix, represented as a dynamic bitset.
using Row=qepg_bits::DynamicBitset;

/**
 * @brief Multiply two GF(2) matrices represented as vectors of Row bitsets.
 *
 * Computes mat1 * mat2 over GF(2) using transpose-and-AND-popcount.
 * Entry (i,j) of the result is the parity of the bitwise AND of
 * row i of mat1 and column j of mat2.
 *
 * @param mat1 Left matrix with dimensions (R1 x C1).
 * @param mat2 Right matrix with dimensions (C1 x C2).
 * @return Result matrix with dimensions (R1 x C2).
 */
std::vector<Row> bitset_matrix_multiplication(const std::vector<Row>& mat1,const std::vector<Row>& mat2);


/**
 * @brief Compute popcount of the bitwise AND of two bitsets.
 *
 * Equivalent to (a & b).count(). Used for GF(2) inner products:
 * the parity of the result gives the dot product mod 2.
 *
 * @tparam Bitset Any type supporting operator& and count() (e.g. boost::dynamic_bitset).
 * @param a First bitset operand.
 * @param b Second bitset operand.
 * @return Number of bits set in (a & b).
 */
template <typename Bitset>
inline std::size_t and_popcount(const Bitset& a, const Bitset&b)
{
    return (a&b).count();
}


/**
 * @brief Print a single GF(2) row to stdout.
 *
 * @tparam BitRow Type supporting .size() and .test(i) (e.g. boost::dynamic_bitset).
 * @param row   The bit row to print.
 * @param zero  Character to display for a 0 bit (default '0').
 * @param one   Character to display for a 1 bit (default '1').
 */
template<class BitRow>
void inline print_bit_row(const BitRow& row,
                      char zero = '0', char one='1')
{
    const std::size_t cols= row.size();
    for(std::size_t c=0; c<cols;++c){
        std::cout<<(row.test(c)? one: zero);
    }
    std::cout<<"\n";
}


/**
 * @brief Print a single GF(2) row to an arbitrary output stream.
 *
 * @tparam BitRow Type supporting .size() and .test(i).
 * @param row   The bit row to print.
 * @param out   Output stream to write to.
 * @param zero  Character to display for a 0 bit (default '0').
 * @param one   Character to display for a 1 bit (default '1').
 */
template<class BitRow>
inline void print_bit_row(const BitRow& row,
                          std::ostream& out,
                          char zero = '0', char one = '1')
{
    const std::size_t cols = row.size();
    for (std::size_t c = 0; c < cols; ++c) {
        out << (row.test(c) ? one : zero);
    }
    out.put('\n');
}



/**
 * @brief Print a GF(2) matrix (vector of bit rows) to stdout.
 *
 * Validates that all rows have the same width before printing.
 *
 * @tparam BitRow Type supporting .size() and .test(i).
 * @param rows  The matrix rows to print.
 * @param zero  Character to display for a 0 bit (default '0').
 * @param one   Character to display for a 1 bit (default '1').
 */
template<class BitRow>
void print_bit_matrix(const std::vector<BitRow>& rows,
                      char zero = '0', char one='1')
{
    if(rows.empty()) return;

    const std::size_t cols= rows.front().size();
    for(const auto& r: rows){
        if(r.size()!=cols){
            std::cerr<<"[print_bit_matrix] row width mismatach\n";
            return;
        }
        for(std::size_t c=0; c<cols;++c){
            std::cout<<(r.test(c)? one: zero);
        }
        std::cout<<"\n";
    }
}



/**
 * @brief The Quantum Error Propagation Graph.
 *
 * Given a compiled Clifford circuit, constructs a binary matrix that maps
 * every possible single-location Pauli error (X, Y, or Z on each noise
 * location) to the set of detectors and the observable it flips.
 *
 * The core data structure is the parity propagation matrix (transpose),
 * with dimensions 3*num_noise rows (X block, Y block, Z block) by
 * (num_detectors + 1) columns. Column i < num_detectors corresponds to
 * detector i; column num_detectors is the logical observable.
 *
 * The backward_graph_construction() method processes gates in reverse
 * circuit order, propagating Pauli frames through Cliffords to build
 * this matrix without explicit matrix multiplication.
 */
class QEPG{


    public:

        /// @brief Default constructor. Creates an empty graph.
        QEPG();

        /**
         * @brief Construct a QEPG from a compiled Clifford circuit.
         *
         * Initializes the parity propagation matrix transpose with dimensions
         * (3 * total_noise) rows by (total_detectors + 1) columns, all zero.
         *
         * @param othercircuit  The compiled Clifford circuit.
         * @param total_detectors Number of detectors in the circuit.
         * @param total_noise     Number of depolarizing noise locations.
         */
        QEPG(const clifford::cliffordcircuit& othercircuit, size_t total_detectors, size_t total_noise);

        /// @brief Destructor.
        ~QEPG();

        /// Move operations (needed because FlatBitTable is non-copyable).
        QEPG(QEPG&&) = default;
        QEPG& operator=(QEPG&&) = default;
        QEPG(const QEPG&) = delete;
        QEPG& operator=(const QEPG&) = delete;

        /**
         * @brief Build the error propagation graph by backward traversal of the circuit.
         *
         * Processes every gate in reverse order. For each gate type:
         * - DEPOLARIZE1: records the current X/Y/Z propagation state of the
         *   affected qubit into the corresponding rows of parityPropMatrixTranspose_.
         * - M (measurement): sets detector and observable bits in the X and Y
         *   propagation vectors for the measured qubit.
         * - R (reset): clears all propagation vectors for that qubit.
         * - cnot: applies the symplectic CNOT update (X propagates forward on
         *   control, Z propagates backward on target, Y mixes both).
         * - h (Hadamard): swaps X and Z propagation vectors.
         *
         * After this method returns, get_parityPropMatrixTrans() provides the
         * complete noise-to-detector/observable mapping.
         */
        void backward_graph_construction();

        /**
         * @brief Print the detector matrix and parity propagation matrices to stdout.
         * @param zero Character for 0 bits (default '0').
         * @param one  Character for 1 bits (default '1').
         */
        void print_detectorMatrix(char zero = '0', char one='1') const;

        /**
         * @brief Get the detector matrix (noise location to measurement propagation).
         * @return Const reference to the detector matrix rows.
         */
        const std::vector<Row>& get_detectorMatrix() const noexcept;

        /**
         * @brief Get the transpose of the detector matrix.
         * @return Const reference to the transposed detector matrix rows.
         */
        const std::vector<Row>& get_dectorMatrixTrans() const noexcept;

        /**
         * @brief Get the parity propagation matrix.
         *
         * Dimensions: (num_detectors + 1) rows by (3 * num_noise) columns.
         * This is the transpose of the primary working matrix.
         *
         * @return Const reference to the parity propagation matrix rows.
         */
        const std::vector<Row>& get_parityPropMatrix() const noexcept;

        /**
         * @brief Get the transpose of the parity propagation matrix.
         *
         * Dimensions: (3 * num_noise) rows by (num_detectors + 1) columns.
         * Rows [0, N), [N, 2N), [2N, 3N) correspond to X, Y, Z Pauli errors
         * respectively, where N = num_noise. The last column is the observable.
         *
         * This is the primary output of backward_graph_construction().
         *
         * @return Const reference to the transposed parity propagation matrix rows.
         */
        const std::vector<Row>& get_parityPropMatrixTrans() const noexcept;

        /// @brief Get the total number of noise locations.
        /// @return Reference to the noise count.
        const size_t& get_total_noise() const noexcept;

        /// @brief Get the total number of detectors.
        /// @return Reference to the detector count.
        const size_t& get_total_detector() const noexcept;

        /// @brief Get the flat (contiguous, SIMD-aligned) parity propagation matrix transpose.
        const qepg_bits::FlatBitTable& get_parityPropMatrixTransFlat() const noexcept;


    private:

        clifford::cliffordcircuit circuit_;        ///< The compiled Clifford circuit.
        std::size_t total_detectors_ = 0;          ///< Number of detectors.
        std::size_t total_noise_=0;                ///< Number of noise locations.

        std::vector<Row> parityPropMatrix_;        ///< Parity propagation matrix: (num_detectors+1) x (3*num_noise).

        std::vector<Row> parityPropMatrixTranspose_;  ///< Transposed parity propagation matrix: (3*num_noise) x (num_detectors+1).

        // Legacy fields used only by compute_parityPropMatrix() path.
        std::vector<Row> detectorMatrix_;
        std::vector<Row> detectorMatrixTranspose_;

        /// Contiguous, SIMD-aligned version of parityPropMatrixTranspose_.
        qepg_bits::FlatBitTable parityPropMatrixTransFlat_;

        /**
         * @brief Compute the parity propagation matrix via explicit matrix multiplication.
         *
         * Builds the parity group matrix from detector/observable definitions,
         * then multiplies it by the detector matrix. This is the older code path;
         * backward_graph_construction() achieves the same result more efficiently.
         */
        void compute_parityPropMatrix();

        /// Build the FlatBitTable from the vector<Row> parity propagation matrix transpose.
        void build_flat_parity_table();
};
}


#endif
