#ifndef SAMPLE_HPP
#define SAMPLE_HPP
#pragma once


#include <cstddef>
#include <ostream>
#include <string>
#include <vector>
#include <random>
#include <unordered_set>  
#include "QEPG.hpp"


namespace SAMPLE{



const size_t PAULIX = 1;
const size_t PAULIY = 2;
const size_t PAULIZ = 3;

struct singlePauli{
    size_t qindex;
    size_t type;
};



class sampler{


    public:
    
    /*---------------------------------------ctor----------*/
        sampler();
        explicit sampler(size_t num_total_paulierror);
        ~sampler();


    /*---------------------------------------Sample one vector with fixed weight----------*/        
        std::vector<size_t> sample_fixed_one_two_three(size_t weight);


        /*
        Generate a single sample with weight error by Floyd method. 
        
        */
        std::vector<singlePauli> generate_sample_Floyd(size_t weight);


        QEPG::Row calculate_output_from_one_sample(const QEPG::QEPG& graph,std::vector<singlePauli> sample);


    private:


        size_t     num_total_pauliError_;



};
}


#endif




