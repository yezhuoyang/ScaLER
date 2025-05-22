#include "sampler.hpp"
#include "chrono"
#include <thread>


namespace SAMPLE{

/*---------------------------------------ctor----------*/
sampler::sampler()=default;

sampler::sampler(size_t num_total_paulierror):num_total_pauliError_(num_total_paulierror){};

sampler::~sampler()=default;

/*---------------------------------------Sample one vector with fixed weight----------*/        


inline std::vector<singlePauli> sampler::generate_sample_Floyd(size_t weight, std::mt19937& gen){

    std::vector<singlePauli> result;
    result.reserve(weight);

    // Uniform distribution in the range [1, 6]
    std::uniform_int_distribution<> posdistrib(0, num_total_pauliError_-1);
    std::uniform_int_distribution<> typedistrib(1, 3);
    std::unordered_set<size_t> usedpos;


    while(result.size()<weight){
        size_t newpos=(size_t)posdistrib(gen);
        if(usedpos.insert(newpos).second){
            result.emplace_back(std::move(singlePauli{ newpos, (size_t)typedistrib(gen)}));  // OK
        }
    }
    return result;
}



void sampler::generate_many_output_samples(const QEPG::QEPG& graph,std::vector<QEPG::Row>& samplecontainer, size_t pauliweight, size_t samplenumber){
    //samplecontainer.reserve(samplenumber);
    samplecontainer.resize(samplenumber);

    static const std::uint64_t global_seed = std::random_device{}();   // log this if you need to replay

    // 2. Parallel region
    #pragma omp parallel
    {
        thread_local std::mt19937 rng{
            static_cast<std::mt19937::result_type>(
                global_seed ^                                    // same run → same base
                std::hash<std::thread::id>{}(std::this_thread::get_id()))  // thread-specific part
        };


        // 3. Work-share the loop
        #pragma omp for schedule(static)
        for (long long i = 0; i < static_cast<long long>(samplenumber); ++i) {
            auto sample = generate_sample_Floyd(pauliweight, rng);
            samplecontainer[i] = calculate_parity_output_from_one_sample(graph, sample);
        }
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
    std::random_device rd;
    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    // Mersenne Twister engine seeded with rd()
    std::mt19937 gen(seed);
    for(size_t i=0;i<samplenumber;i++){
        std::vector<singlePauli> sample = generate_sample_Floyd(pauliweight,gen);
        noisecontainer.push_back(std::move(sample));
        samplecontainer.push_back(std::move(calculate_parity_output_from_one_sample(graph,sample)));
    }
}





}

