#include "LERcalculator.hpp"




namespace LERcalculator{

std::vector<QEPG::Row> return_samples(std::string prog_str,size_t weight, size_t shots){
    clifford::cliffordcircuit c;
    c.compile_from_rewrited_stim_string(prog_str);


    QEPG::QEPG graph(c,c.get_num_detector(),c.get_num_noise());
    graph.backward_graph_construction();

    size_t qubitnum=c.get_num_qubit();

    SAMPLE::sampler sampler(qubitnum);

    std::vector<QEPG::Row> samplecontainer;
    sampler.generate_many_output_samples(graph,samplecontainer,weight,shots);


    // for(QEPG::Row parityresult: samplecontainer){
    //     QEPG::print_bit_row(parityresult);
    // }
    return samplecontainer;
}

}