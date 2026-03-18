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
#include <bit>
#include <thread>


namespace SAMPLE{

/*---------------------------------------ctor----------*/
sampler::sampler()=default;

sampler::sampler(size_t num_total_paulierror)
    : num_total_pauliError_(num_total_paulierror),
      collision_bitmap_(num_total_paulierror, 0),
      scratch_sample_() {
    scratch_sample_.reserve(num_total_paulierror);
};

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
inline std::size_t sampler::generate_sample_removal(size_t weight, std::mt19937& gen){
    // Mark all positions as active using bitmap
    std::memset(collision_bitmap_.data(), 1, num_total_pauliError_);
    size_t remaining = num_total_pauliError_;

    // Remove random positions until only `weight` remain
    std::uniform_int_distribution<size_t> posdistrib(0, num_total_pauliError_-1);
    while(remaining > weight){
        size_t pos = posdistrib(gen);
        if(collision_bitmap_[pos]){
            collision_bitmap_[pos] = 0;
            --remaining;
        }
    }

    // Collect remaining positions into scratch buffer
    scratch_sample_.clear();
    std::uniform_int_distribution<size_t> typedistrib(1, 3);
    for(size_t i = 0; i < num_total_pauliError_; ++i){
        if(collision_bitmap_[i]){
            scratch_sample_.push_back(singlePauli{i, typedistrib(gen)});
        }
    }
    // Clear bitmap for positions we used
    std::memset(collision_bitmap_.data(), 0, num_total_pauliError_);
    return scratch_sample_.size();
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
inline std::size_t sampler::generate_sample_Floyd(size_t weight, std::mt19937& gen){
    // Hybrid strategy: use removal when weight > half of total
    if(weight > num_total_pauliError_ / 2){
        return generate_sample_removal(weight, gen);
    }

    // Addition strategy with bitmap collision detection (no heap allocation)
    scratch_sample_.clear();

    std::uniform_int_distribution<size_t> posdistrib(0, num_total_pauliError_-1);
    std::uniform_int_distribution<size_t> typedistrib(1, 3);

    while(scratch_sample_.size() < weight){
        size_t newpos = posdistrib(gen);
        if(!collision_bitmap_[newpos]){
            collision_bitmap_[newpos] = 1;
            scratch_sample_.push_back(singlePauli{newpos, typedistrib(gen)});
        }
    }
    // Clear only the used positions (O(weight) not O(N))
    for(std::size_t i = 0; i < scratch_sample_.size(); ++i){
        collision_bitmap_[scratch_sample_[i].qindex] = 0;
    }
    return scratch_sample_.size();
}


/**
 * @brief Generate a single error sample using Bernoulli (Monte Carlo) sampling.
 *
 * Iterates over the first ErrorSize noise locations. Each location independently
 * has an error with probability error_prob. Activated locations receive a
 * uniformly random Pauli type.
 */
inline std::size_t sampler::generate_sample_Monte(double error_prob,size_t ErrorSize,std::mt19937& gen){
    scratch_sample_.clear();
    std::uniform_int_distribution<size_t> typedistrib(1, 3);
    std::bernoulli_distribution dist(error_prob);
    for(size_t pos=0;pos<ErrorSize;++pos){
        if(dist(gen)){
            scratch_sample_.push_back(singlePauli{pos, typedistrib(gen)});
        }
    }
    return scratch_sample_.size();
}


/*----------- Xoshiro256++ overloads -----------*/

inline std::size_t sampler::generate_sample_removal(size_t weight, Xoshiro256pp& gen){
    std::memset(collision_bitmap_.data(), 1, num_total_pauliError_);
    size_t remaining = num_total_pauliError_;

    while(remaining > weight){
        size_t pos = gen.bounded(num_total_pauliError_);
        if(collision_bitmap_[pos]){
            collision_bitmap_[pos] = 0;
            --remaining;
        }
    }

    scratch_sample_.clear();
    for(size_t i = 0; i < num_total_pauliError_; ++i){
        if(collision_bitmap_[i]){
            scratch_sample_.push_back(singlePauli{i, 1 + gen.bounded(3)});
        }
    }
    std::memset(collision_bitmap_.data(), 0, num_total_pauliError_);
    return scratch_sample_.size();
}

inline std::size_t sampler::generate_sample_Floyd(size_t weight, Xoshiro256pp& gen){
    if(weight > num_total_pauliError_ / 2){
        return generate_sample_removal(weight, gen);
    }

    scratch_sample_.clear();

    while(scratch_sample_.size() < weight){
        size_t newpos = gen.bounded(num_total_pauliError_);
        if(!collision_bitmap_[newpos]){
            collision_bitmap_[newpos] = 1;
            scratch_sample_.push_back(singlePauli{newpos, 1 + gen.bounded(3)});
        }
    }
    for(std::size_t i = 0; i < scratch_sample_.size(); ++i){
        collision_bitmap_[scratch_sample_[i].qindex] = 0;
    }
    return scratch_sample_.size();
}

inline std::size_t sampler::generate_sample_Monte(double error_prob, size_t ErrorSize, Xoshiro256pp& gen){
    scratch_sample_.clear();
    // Use threshold on upper 32 bits for speed (sufficient precision for error rates)
    const std::uint32_t threshold = static_cast<std::uint32_t>(error_prob * 4294967296.0);
    for(size_t pos = 0; pos < ErrorSize; ++pos){
        // Use upper 32 bits of random value for threshold test
        if(static_cast<std::uint32_t>(gen() >> 32) < threshold){
            scratch_sample_.push_back(singlePauli{pos, 1 + (gen() % 3)});
        }
    }
    return scratch_sample_.size();
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
    samplecontainer.resize(samplenumber);

    static const std::uint64_t global_seed = std::random_device{}();
    const size_t n_err = num_total_pauliError_;
    const auto& flat = graph.get_parityPropMatrixTransFlat();
    const std::size_t n_noise = flat.n_rows() / 3;
    const std::size_t n_words = flat.words_per_row();
    const std::size_t n_cols = flat.n_cols();

    #pragma omp parallel
    {
        sampler local_sampler(n_err);
        Xoshiro256pp rng(global_seed ^ std::hash<std::thread::id>{}(std::this_thread::get_id()));

        // Per-thread aligned result buffer
        std::uint64_t* result_buf = simd::aligned_alloc_u64(flat.stride_words());

        #pragma omp for schedule(static)
        for (long long i = 0; i < static_cast<long long>(samplenumber); ++i) {
            std::size_t n = local_sampler.generate_sample_Floyd(pauliweight, rng);

            // Zero result, accumulate XOR via SIMD, copy to Row
            simd::zero_words(result_buf, n_words);
            local_sampler.calculate_parity_output_flat(
                flat, n_noise, local_sampler.scratch_data(), n, result_buf);

            samplecontainer[i] = QEPG::Row(n_cols, result_buf, n_words);
        }

        simd::aligned_free(result_buf);
    }
}


/**
 * @brief Generate many Monte Carlo samples in parallel.
 *
 * Same parallelization strategy as generate_many_output_samples() but uses
 * Bernoulli sampling (generate_sample_Monte) instead of fixed-weight sampling.
 */
void sampler::generate_many_output_samples_Monte(const QEPG::QEPG& graph,std::vector<QEPG::Row>& samplecontainer,double error_prob, size_t samplenumber){
    samplecontainer.resize(samplenumber);
    size_t total_error=graph.get_total_noise();
    static const std::uint64_t global_seed = std::random_device{}();
    const size_t n_err = num_total_pauliError_;
    const auto& flat = graph.get_parityPropMatrixTransFlat();
    const std::size_t n_noise = flat.n_rows() / 3;
    const std::size_t n_words = flat.words_per_row();
    const std::size_t n_cols = flat.n_cols();

    #pragma omp parallel
    {
        sampler local_sampler(n_err);
        Xoshiro256pp rng(global_seed ^ std::hash<std::thread::id>{}(std::this_thread::get_id()));
        std::uint64_t* result_buf = simd::aligned_alloc_u64(flat.stride_words());

        #pragma omp for schedule(static)
        for (long long i = 0; i < static_cast<long long>(samplenumber); ++i) {
            std::size_t n = local_sampler.generate_sample_Monte(error_prob, total_error, rng);

            simd::zero_words(result_buf, n_words);
            local_sampler.calculate_parity_output_flat(
                flat, n_noise, local_sampler.scratch_data(), n, result_buf);

            samplecontainer[i] = QEPG::Row(n_cols, result_buf, n_words);
        }

        simd::aligned_free(result_buf);
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
    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 gen(static_cast<std::mt19937::result_type>(seed));
    for(size_t i=0;i<samplenumber;i++){
        std::size_t n = generate_sample_Floyd(pauliweight, gen);
        samplecontainer.push_back(calculate_parity_output_from_one_sample(graph, scratch_data(), n));
        noisecontainer.emplace_back(scratch_sample_.begin(), scratch_sample_.begin() + n);
    }
}


/// Helper: unpack result_buf (uint64_t words) into numpy byte buffers
static inline void unpack_result_to_numpy(
    const std::uint64_t* result_buf,
    std::uint8_t* det_dst,
    std::uint8_t& obs_dst,
    std::size_t n_det,
    std::size_t n_words) noexcept
{
    constexpr std::size_t WORD_BITS = 64;
    const std::size_t full_words = n_det / WORD_BITS;

    for(std::size_t b = 0; b < full_words; ++b){
        std::uint64_t word = result_buf[b];
        for(std::size_t k = 0; k < WORD_BITS; ++k, word >>= 1)
            det_dst[b * WORD_BITS + k] = static_cast<std::uint8_t>(word & 1);
    }

    const std::size_t rem = n_det % WORD_BITS;
    if(rem){
        std::uint64_t word = result_buf[full_words];
        for(std::size_t k = 0; k < rem; ++k, word >>= 1)
            det_dst[full_words * WORD_BITS + k] = static_cast<std::uint8_t>(word & 1);
    }

    // Observable is the last bit (at position n_det)
    obs_dst = static_cast<std::uint8_t>(
        (result_buf[n_det / WORD_BITS] >> (n_det % WORD_BITS)) & 1);
}


void sampler::generate_many_output_samples_to_numpy(
    const QEPG::QEPG& graph,
    std::uint8_t* det_buf, std::uint8_t* obs_buf,
    std::size_t n_det, size_t pauliweight, size_t samplenumber)
{
    static const std::uint64_t global_seed = std::random_device{}();
    const size_t n_err = num_total_pauliError_;
    const auto& flat = graph.get_parityPropMatrixTransFlat();
    const std::size_t n_noise = flat.n_rows() / 3;
    const std::size_t n_words = flat.words_per_row();

    #pragma omp parallel
    {
        sampler local_sampler(n_err);
        Xoshiro256pp rng(global_seed ^ std::hash<std::thread::id>{}(std::this_thread::get_id()));
        std::uint64_t* result_buf = simd::aligned_alloc_u64(flat.stride_words());

        #pragma omp for schedule(static)
        for (long long i = 0; i < static_cast<long long>(samplenumber); ++i) {
            std::size_t n = local_sampler.generate_sample_Floyd(pauliweight, rng);

            simd::zero_words(result_buf, n_words);
            local_sampler.calculate_parity_output_flat(
                flat, n_noise, local_sampler.scratch_data(), n, result_buf);

            unpack_result_to_numpy(result_buf,
                det_buf + i * n_det, obs_buf[i], n_det, n_words);
        }

        simd::aligned_free(result_buf);
    }
}


void sampler::generate_many_output_samples_Monte_to_numpy(
    const QEPG::QEPG& graph,
    std::uint8_t* det_buf, std::uint8_t* obs_buf,
    std::size_t n_det, double error_prob, size_t samplenumber)
{
    static const std::uint64_t global_seed = std::random_device{}();
    const size_t n_err = num_total_pauliError_;
    const size_t total_error = graph.get_total_noise();
    const auto& flat = graph.get_parityPropMatrixTransFlat();
    const std::size_t n_noise = flat.n_rows() / 3;
    const std::size_t n_words = flat.words_per_row();

    #pragma omp parallel
    {
        sampler local_sampler(n_err);
        Xoshiro256pp rng(global_seed ^ std::hash<std::thread::id>{}(std::this_thread::get_id()));
        std::uint64_t* result_buf = simd::aligned_alloc_u64(flat.stride_words());

        #pragma omp for schedule(static)
        for (long long i = 0; i < static_cast<long long>(samplenumber); ++i) {
            std::size_t n = local_sampler.generate_sample_Monte(error_prob, total_error, rng);

            simd::zero_words(result_buf, n_words);
            local_sampler.calculate_parity_output_flat(
                flat, n_noise, local_sampler.scratch_data(), n, result_buf);

            unpack_result_to_numpy(result_buf,
                det_buf + i * n_det, obs_buf[i], n_det, n_words);
        }

        simd::aligned_free(result_buf);
    }
}


}
