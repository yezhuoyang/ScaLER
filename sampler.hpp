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
        explicit sampler(int num_total_paulierror);
        ~sampler();

        
    /*---------------------------------------Sample one vector with fixed weight----------*/        
        std::vector<int> sample_fixed_one_two_three(int weight);


    private:

        int     num_total_pauliError_;



};
}


#endif




