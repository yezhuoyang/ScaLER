#ifndef SAMPLE_HPP
#define SAMPLE_HPP
#pragma once


#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

namespace SAMPLE{



class sampler{


    public:
    
    /*---------------------------------------ctor----------*/
        sampler();
        explicit sampler(size_t num_total_paulierror);
        ~sampler();


    /*---------------------------------------Sample one vector with fixed weight----------*/        
        std::vector<size_t> sample_fixed_one_two_three(size_t weight);


    private:

        size_t     num_total_pauliError_;



};
}


#endif




