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



class QEPG{


    public:
    
        QEPG();
        QEPG(clifford::cliffordcircuit othercircuit, size_t total_detectors, size_t total_noise);
        ~QEPG();

        void backward_graph_construction();


    private:

        clifford::cliffordcircuit circuit;
        std::size_t total_detectors_ = 0;
        std::size_t total_noise_=0;

        std::size_t COLS = 3*total_noise_;
        using Row=boost::dynamic_bitset<>;


        std::vector<Row> detectorMatrix_;

        
};
}


#endif




