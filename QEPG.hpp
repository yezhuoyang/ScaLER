#ifndef QEPG_HPP
#define QEPG_HPP
#pragma once


#include <cstddef>
#include <ostream>
#include <string>
#include <bitset>  
#include <boost/dynamic_bitset.hpp>
#include "clifford.hpp"

namespace QEPG{

using Row=boost::dynamic_bitset<>;

class QEPG{


    public:
    
        QEPG();
        QEPG(clifford::cliffordcircuit othercircuit, size_t total_detectors, size_t total_noise);
        ~QEPG();

        void backward_graph_construction();

        void print_detectorMatrix(char zero = '0', char one='1') const;

        const std::vector<Row>& get_detectorMatrix() const noexcept; 

        const std::vector<Row>& get_dectorMatrixTrans() const noexcept;

        const std::vector<Row>& get_parityPropMatrix() const noexcept; 

    private:

        clifford::cliffordcircuit circuit;
        std::size_t total_detectors_ = 0;
        std::size_t total_noise_=0;

        std::size_t COLS = 3*total_noise_;

        std::vector<Row> parityPropMatrix_;        

        std::vector<Row> detectorMatrix_;

        std::vector<Row> detectorMatrixTranspose_;
        
        void compute_parityPropMatrix();
};
}


#endif




