#ifndef QEPG_HPP
#define QEPG_HPP
#pragma once


#include <cstddef>
#include <ostream>
#include <string>
#include <bitset>  
#include <boost/dynamic_bitset.hpp>
#include "clifford.hpp"
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace QEPG{

using Row=boost::dynamic_bitset<>;
using SpMat = Eigen::SparseMatrix<bool,Eigen::RowMajor>;
using SpVec = Eigen::SparseVector<bool>;

std::vector<Row> bitset_matrix_multiplication(const std::vector<Row>& mat1,const std::vector<Row>& mat2);


template <typename Bitset>
inline std::size_t and_popcount(const Bitset& a, const Bitset&b)
{
    return (a&b).count();
}



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


template<class BitRow>
inline void print_bit_row(const BitRow& row,
                          std::ostream& out,   // default = console
                          char zero = '0', char one = '1')
{
    const std::size_t cols = row.size();
    for (std::size_t c = 0; c < cols; ++c) {
        out << (row.test(c) ? one : zero);
    }
    out.put('\n');
}





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



class QEPG{


    public:
    
        QEPG();
        QEPG(clifford::cliffordcircuit othercircuit, size_t total_detectors, size_t total_noise);
        ~QEPG();

        void backward_graph_construction();

        void backward_graph_construction_Eigen();

        void print_detectorMatrix(char zero = '0', char one='1') const;

        const std::vector<Row>& get_detectorMatrix() const noexcept; 

        const std::vector<Row>& get_dectorMatrixTrans() const noexcept;

        const std::vector<Row>& get_parityPropMatrix() const noexcept; 

        const std::vector<Row>& get_parityPropMatrixTrans() const noexcept; 


        const SpMat& get_detectorMatrix_Eigen() const noexcept; 

        const SpMat& get_dectorMatrixTrans_Eigen() const noexcept;

        const SpMat& get_parityPropMatrix_Eigen() const noexcept; 

        const SpMat& get_parityPropMatrixTrans_Eigen() const noexcept; 


    private:

        clifford::cliffordcircuit circuit_;
        std::size_t total_detectors_ = 0;
        std::size_t total_noise_=0;

        std::size_t COLS = 3*total_noise_;

        std::vector<Row> parityPropMatrix_;    

        SpMat parityPropMatrixEigen_;
        
        std::vector<Row> parityPropMatrixTranspose_;     

        SpMat parityPropMatrixTransposeEigen_;

        std::vector<Row> detectorMatrix_;

        SpMat detectorMatrixEigen_;

        std::vector<Row> detectorMatrixTranspose_;
        
        SpMat detectorMatrixTransposeEigen_;


        void compute_parityPropMatrix();

        void compute_parityPropMatrix_Eigen();
};
}


#endif




