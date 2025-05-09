#include "LERcalculator.hpp"




namespace LERcalculator{




void convert_bitset_row_to_boolean(std::vector<std::vector<bool>>& result,const std::vector<QEPG::Row>& samplecontainer){
        result.reserve(samplecontainer.size()); // Reserve space

        // Convert each boost::dynamic_bitset<> to std::vector<bool>
        for (const auto& bitset_row : samplecontainer) {
            std::vector<bool> bool_row(bitset_row.size());
            for (size_t i = 0; i < bitset_row.size(); ++i) {
                bool_row[i] = bitset_row[i]; // Access individual bits
            }
            result.push_back(bool_row);
        }
}



 std::vector<std::vector<bool>> return_samples(std::string prog_str,size_t weight, size_t shots){
    clifford::cliffordcircuit c;
    c.compile_from_rewrited_stim_string(prog_str);

    QEPG::QEPG graph(c,c.get_num_detector(),c.get_num_noise());
    graph.backward_graph_construction();

    size_t qubitnum=c.get_num_qubit();

    SAMPLE::sampler sampler(qubitnum);

    std::vector<QEPG::Row> samplecontainer;
    sampler.generate_many_output_samples(graph,samplecontainer,weight,shots);

    std::vector<std::vector<bool>> result;
    convert_bitset_row_to_boolean(result,samplecontainer);

    // for(QEPG::Row parityresult: samplecontainer){
    //     QEPG::print_bit_row(parityresult);
    // }
    return result;
}

}