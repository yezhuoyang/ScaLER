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


inline QEPG::Row sampler::calculate_output_from_one_sample(const QEPG::QEPG& graph,std::vector<singlePauli> sample){

    const auto&dm=graph.get_dectorMatrixTrans();

    const std::size_t n_rows=dm.size();
    const std::size_t n_cols=n_rows ? dm[0].size():0;

    QEPG::Row result(n_cols);

    for(singlePauli noise: sample){
        size_t pos=noise.qindex;
        size_t type=noise.type;
        if(type==SAMPLE::PAULIX){
            result^=dm[3*pos];
        }
        else if(type==SAMPLE::PAULIY){
            result^=dm[3*pos+1];
        }
        else if(type==SAMPLE::PAULIZ){
            result^=dm[3*pos+2];
        }
    }

    return result;
}




inline QEPG::Row sampler::calculate_parity_output_from_one_sample(const QEPG::QEPG& graph,const std::vector<singlePauli>& sample){
    const auto&dm=graph.get_parityPropMatrixTrans();

    const std::size_t n_rows=dm.size();
    const std::size_t n_cols=n_rows ? dm[0].size():0;

    QEPG::Row result(n_cols);

    for(singlePauli noise: sample){
        size_t pos=noise.qindex;
        size_t type=noise.type;
        if(type==SAMPLE::PAULIX){
            result^=dm[3*pos];
        }
        else if(type==SAMPLE::PAULIY){
            result^=dm[3*pos+1];
        }
        else if(type==SAMPLE::PAULIZ){
            result^=dm[3*pos+2];
        }
    }

    return result;
}





void sampler::generate_many_output_samples(const QEPG::QEPG& graph,std::vector<QEPG::Row>& samplecontainer, size_t pauliweight, size_t samplenumber){
    samplecontainer.reserve(samplenumber);
    for(size_t i=0;i<samplenumber;i++){
        sampler sampler_instance;
        std::vector<singlePauli> sample = generate_sample_Floyd(pauliweight);
        samplecontainer.emplace_back(calculate_parity_output_from_one_sample(graph,sample));
    }
}




}