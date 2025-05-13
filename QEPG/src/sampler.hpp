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


        inline QEPG::Row calculate_output_from_one_sample(const QEPG::QEPG& graph,std::vector<singlePauli> sample){
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
        


        inline QEPG::Row calculate_parity_output_from_one_sample(const QEPG::QEPG& graph,const std::vector<singlePauli>& sample){
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


        void generate_many_output_samples(const QEPG::QEPG& graph,std::vector<QEPG::Row>& samplecontainer,size_t pauliweight , size_t samplenumber);


    private:


        size_t     num_total_pauliError_;



};




}
#endif




