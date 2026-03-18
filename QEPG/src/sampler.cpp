/**
 * @file sampler.cpp
 * @brief Implementation of random Pauli error sampling methods.
 *
 * Contains the sampling algorithms (Floyd, removal, Monte Carlo) for
 * generating random error configurations, and the batch methods that
 * use OpenMP to parallelize sample generation across threads.
 */

#include "sampler.hpp"
#include "chrono"
#include <thread>


namespace SAMPLE{

/*---------------------------------------ctor----------*/
sampler::sampler()=default;

sampler::sampler(size_t num_total_paulierror):num_total_pauliError_(num_total_paulierror){};

sampler::~sampler()=default;

/*---------------------------------------Sample one vector with fixed weight----------*/


/**
 * @brief Generate a single error sample using removal-based sampling.
 *
 * Creates a set of all noise location indices, then randomly removes
 * positions until exactly `weight` remain. Each remaining position
 * receives a uniformly random Pauli type (X, Y, or Z).
 *
 * @note More efficient than Floyd's insertion method when weight is
 *       close to num_total_pauliError_ (avoids collision overhead).
 */
inline std::vector<singlePauli> sampler::generate_sample_removal(size_t weight, std::mt19937& gen){
    // Create set of all positions
    std::unordered_set<size_t> remaining_positions;
    for(size_t i = 0; i < num_total_pauliError_; i++){
        remaining_positions.insert(i);
    }

    // Remove (num_total_pauliError_ - weight) random positions
    std::uniform_int_distribution<> posdistrib(0, num_total_pauliError_-1);
    size_t num_to_remove = num_total_pauliError_ - weight;

    // Use same collision-based removal (mirrors addition logic)
    while(remaining_positions.size() > weight){
        size_t pos_to_remove = posdistrib(gen);
        remaining_positions.erase(pos_to_remove);  // erase() is no-op if not found
    }

    // Convert remaining positions to result with random Pauli types
    std::vector<singlePauli> result;
    result.reserve(weight);
    std::uniform_int_distribution<> typedistrib(1, 3);

    for(size_t pos : remaining_positions){
        result.emplace_back(singlePauli{pos, (size_t)typedistrib(gen)});
    }

    return result;
}

/**
 * @brief Generate a single error sample using Floyd's insertion algorithm.
 *
 * Selects `weight` distinct noise locations uniformly at random. Uses
 * a hash set to reject collisions. Automatically delegates to
 * generate_sample_removal() when weight > num_total_pauliError_ / 2
 * to avoid excessive collisions.
 *
 * Each selected position receives a uniformly random Pauli type (1=X, 2=Y, 3=Z).
 */
inline std::vector<singlePauli> sampler::generate_sample_Floyd(size_t weight, std::mt19937& gen){
    // Hybrid strategy: use removal when weight > half of total
    // This avoids collision inefficiency for high weights
    if(weight > num_total_pauliError_ / 2){
        return generate_sample_removal(weight, gen);
    }

    // Original addition strategy for low weights
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


/**
 * @brief Generate a single error sample using Bernoulli (Monte Carlo) sampling.
 *
 * Iterates over the first ErrorSize noise locations. Each location independently
 * has an error with probability error_prob. Activated locations receive a
 * uniformly random Pauli type.
 */
inline std::vector<singlePauli> sampler::generate_sample_Monte(double error_prob,size_t ErrorSize,std::mt19937& gen){
    std::vector<singlePauli> result;
    result.reserve(size_t(error_prob*ErrorSize));
    std::uniform_int_distribution<> typedistrib(1, 3);
    std::bernoulli_distribution dist(error_prob);
    for(size_t pos=0;pos<ErrorSize;++pos){
        if(dist(gen)){
            result.emplace_back(std::move(singlePauli{ pos, (size_t)typedistrib(gen)}));
        }
    }
    return result;
}




/**
 * @brief Generate many samples in parallel with fixed Pauli weight.
 *
 * Uses OpenMP to distribute sample generation across threads. Each thread
 * maintains a thread-local Mersenne Twister PRNG seeded from a global seed
 * XORed with the thread ID hash for reproducibility within a run.
 *
 * Each sample is generated via generate_sample_Floyd() and then its
 * detector/observable outcome is computed via calculate_parity_output_from_one_sample().
 */
void sampler::generate_many_output_samples(const QEPG::QEPG& graph,std::vector<QEPG::Row>& samplecontainer, size_t pauliweight, size_t samplenumber){
    //samplecontainer.reserve(samplenumber);
    samplecontainer.resize(samplenumber);

    static const std::uint64_t global_seed = std::random_device{}();   // log this if you need to replay

    // 2. Parallel region
    #pragma omp parallel
    {
        thread_local std::mt19937 rng{
            static_cast<std::mt19937::result_type>(
                global_seed ^                                    // same run -> same base
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


/**
 * @brief Generate many Monte Carlo samples in parallel.
 *
 * Same parallelization strategy as generate_many_output_samples() but uses
 * Bernoulli sampling (generate_sample_Monte) instead of fixed-weight sampling.
 */
void sampler::generate_many_output_samples_Monte(const QEPG::QEPG& graph,std::vector<QEPG::Row>& samplecontainer,double error_prob, size_t samplenumber){
    //samplecontainer.reserve(samplenumber);
    samplecontainer.resize(samplenumber);
    size_t total_error=graph.get_total_noise();
    static const std::uint64_t global_seed = std::random_device{}();   // log this if you need to replay
    // 2. Parallel region
    #pragma omp parallel
    {
        thread_local std::mt19937 rng{
            static_cast<std::mt19937::result_type>(
                global_seed ^                                    // same run -> same base
                std::hash<std::thread::id>{}(std::this_thread::get_id()))  // thread-specific part
        };
        // 3. Work-share the loop
        #pragma omp for schedule(static)
        for (long long i = 0; i < static_cast<long long>(samplenumber); ++i) {
            auto sample = generate_sample_Monte(error_prob,total_error, rng);
            samplecontainer[i] = calculate_parity_output_from_one_sample(graph, sample);
        }
    }
}


/**
 * @brief Enumerate all error configurations of a given weight (stub).
 *
 * Intended to enumerate all C(N, weight) * 3^weight configurations
 * using recursion. Currently not implemented.
 */
void sampler::generate_all_samples_with_fixed_weight(const QEPG::QEPG& graph,std::vector<QEPG::Row>& samplecontainer,size_t pauliweight){

}


/**
 * @brief Generate many samples and also return the noise vectors (single-threaded).
 *
 * Unlike the parallel batch methods, this version stores both the computed
 * detector/observable outcomes and the raw noise configurations that produced them.
 * Useful for debugging or for algorithms that need to inspect which specific
 * errors caused a given syndrome.
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
        samplecontainer.push_back(std::move(calculate_parity_output_from_one_sample(graph,sample)));
        noisecontainer.push_back(std::move(sample));
    }
}




}
