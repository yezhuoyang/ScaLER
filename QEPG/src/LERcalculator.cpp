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



 std::vector<std::vector<bool>> return_samples(const std::string& prog_str,size_t weight, size_t shots){
    clifford::cliffordcircuit c;
    c.compile_from_rewrited_stim_string(prog_str);

    QEPG::QEPG graph(c,c.get_num_detector(),c.get_num_noise());
    graph.backward_graph_construction();


    SAMPLE::sampler sampler(c.get_num_noise());

    std::vector<QEPG::Row> samplecontainer;

    using clock     = std::chrono::steady_clock;          // monotonic, good for benchmarking
    using microsec  = std::chrono::microseconds;
    auto t0 = clock::now();                               // start timer
    sampler.generate_many_output_samples(graph,samplecontainer,weight,shots);



    auto t1 = clock::now();                               // stop section‑1
    auto compile_us = std::chrono::duration_cast<microsec>(t1 - t0).count();
    std::cout << "[Time to generate these samples:] " << compile_us / 1'000.0 << "ms\n";


    std::vector<std::vector<bool>> result;
    convert_bitset_row_to_boolean(result,samplecontainer);

    // for(QEPG::Row parityresult: samplecontainer){
    //     QEPG::print_bit_row(parityresult);
    // }
    return result;
}



std::vector<std::vector<std::vector<bool>>> return_samples_many_weights(const std::string& prog_str,const std::vector<size_t>& weight, const std::vector<size_t>& shots){
    clifford::cliffordcircuit c;
    c.compile_from_rewrited_stim_string(prog_str);

    QEPG::QEPG graph(c,c.get_num_detector(),c.get_num_noise());
    graph.backward_graph_construction();

    SAMPLE::sampler sampler(c.get_num_noise());


    std::vector<QEPG::Row> samplecontainer;
    std::vector<std::vector<bool>> tmpresult;
    std::vector<std::vector<std::vector<bool>>> result;


    for(size_t i=0;i<weight.size();++i){
        std::cout<<"Weight="<<weight[i]<<"\n";
        samplecontainer.clear();
        tmpresult.clear();
        sampler.generate_many_output_samples(graph,samplecontainer,weight[i],shots[i]);
        convert_bitset_row_to_boolean(tmpresult,samplecontainer);
        result.emplace_back(tmpresult);
    }
    // for(QEPG::Row parityresult: samplecontainer){
    //     QEPG::print_bit_row(parityresult);
    // }
    return result;
}



std::vector<std::vector<bool>> return_detector_matrix(const std::string& prog_str){
    clifford::cliffordcircuit c;
    c.compile_from_rewrited_stim_string(prog_str);

    QEPG::QEPG graph(c,c.get_num_detector(),c.get_num_noise());
    c.print_circuit();
    graph.backward_graph_construction();
    graph.print_detectorMatrix();
    const std::vector<QEPG::Row>& parityMtrans=graph.get_parityPropMatrixTrans();
    const size_t row_size=parityMtrans.size();
    const size_t col_size=parityMtrans[0].size();


    // 2. Allocate the whole target matrix in one go
    std::vector<std::vector<bool>> result(row_size,std::vector<bool>(col_size));

    for(size_t row=0;row<row_size;row++){
        for(size_t column=0;column<col_size;column++){
            result[row][column]=parityMtrans[row][column];
        }
    }
    return result;
}




}