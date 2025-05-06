#include "sampler.hpp"



namespace SAMPLE{

/*---------------------------------------ctor----------*/
sampler::sampler()=default;

sampler::sampler(size_t num_total_paulierror):num_total_pauliError_(num_total_paulierror){};

sampler::~sampler()=default;

/*---------------------------------------Sample one vector with fixed weight----------*/        
std::vector<size_t> sampler::sample_fixed_one_two_three(size_t weight){
    
}


std::vector<singlePauli> sampler::generate_sample_Floyd(size_t weight){

    std::vector<singlePauli> result;

    // Seed source for the random number engine
    std::random_device rd;

    // Mersenne Twister engine seeded with rd()
    std::mt19937 gen(rd());

    // Uniform distribution in the range [1, 6]
    std::uniform_int_distribution<> posdistrib(0, num_total_pauliError_-1);
    std::uniform_int_distribution<> typedistrib(1, 3);

    std::unordered_set<size_t> usedpos;


    while(result.size()<weight){
        size_t newpos=(size_t)posdistrib(gen);
        if(usedpos.insert(newpos).second){
            result.emplace_back( singlePauli{ newpos, (size_t)typedistrib(gen)} );  // OK
        }
    }
    return result;
}




}