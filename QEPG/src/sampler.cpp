#include "sampler.hpp"
#include "chrono";


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
    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    // Mersenne Twister engine seeded with rd()
    std::mt19937 gen(seed);

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



void sampler::generate_many_output_samples(const QEPG::QEPG& graph,std::vector<QEPG::Row>& samplecontainer, size_t pauliweight, size_t samplenumber){
    samplecontainer.reserve(samplenumber);
    for(size_t i=0;i<samplenumber;i++){
        std::vector<singlePauli> sample = generate_sample_Floyd(pauliweight);
        samplecontainer.push_back(calculate_parity_output_from_one_sample(graph,sample));
    }
}

/*
Enumerat all possible noise vector with fixed weight
Use recursion
*/
void sampler::generate_all_samples_with_fixed_weight(const QEPG::QEPG& graph,std::vector<QEPG::Row>& samplecontainer,size_t pauliweight){

}


/*
In this implementation, we also return the generated random noise vector 
*/
void sampler::generate_many_output_samples_with_noise_vector(const QEPG::QEPG& graph,std::vector<std::vector<singlePauli>>& noisecontainer,std::vector<QEPG::Row>& samplecontainer, size_t pauliweight, size_t samplenumber){
    samplecontainer.reserve(samplenumber);
    noisecontainer.reserve(samplenumber);
    for(size_t i=0;i<samplenumber;i++){
        std::vector<singlePauli> sample = generate_sample_Floyd(pauliweight);
        noisecontainer.push_back(sample);
        samplecontainer.push_back(calculate_parity_output_from_one_sample(graph,sample));
    }
}





}

