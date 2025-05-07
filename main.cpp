#include <iostream>
#include <vector>
#include <iomanip>



#include "QEPG/clifford.hpp"
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <chrono>

std::string read_file_to_string(const std::string& path)
{
    std::ifstream in(path, std::ios::in | std::ios::binary);
    if(!in) throw std::runtime_error("Cannot open "+ path);

    std::ostringstream buffer;
    buffer<< in.rdbuf();
    return buffer.str();
};



#include "QEPG/QEPG.hpp"
#include "QEPG/sampler.hpp"




void test_bitset_pop_count(){
        std::bitset<8> x("11001110");
        std::bitset<8> y("01100100");
        std::cout<<(x^y)<<"\n";
        size_t pop_count=QEPG::and_popcount(x,y);
        std::cout<<pop_count<<"\n";
}


QEPG::Row make_row(std::initializer_list<int> bits) {
    QEPG::Row r(bits.size());
    std::size_t idx = 0;
    for (int bit : bits) r[idx++] = (bit != 0);
    return r;
}




void test_matrix_multiplication(){
    // Test 1: Identity × A
    std::vector<QEPG::Row> I = {
        make_row({1,0,0}),
        make_row({0,1,0}),
        make_row({0,0,1})
    };
    std::vector<QEPG::Row> A = {
        make_row({1,0,1}),
        make_row({0,1,1}),
        make_row({1,1,0})
    };

    auto C = QEPG::bitset_matrix_multiplication(I, A);
    QEPG::print_bit_matrix(C);
    assert(C == A);                // passes if implementation is correct
}



/*
Benchmard the running time of generating 1 million surface code
sample
*/
void benchmark_surface_million_sample(){

    using clock     = std::chrono::steady_clock;          // monotonic, good for benchmarking
    using microsec  = std::chrono::microseconds;


    auto t0 = clock::now();                               // start timer
    clifford::cliffordcircuit c;
    try{
        std::string stim_str=read_file_to_string("C:/Users/yezhu/OneDrive/Documents/GitHub/Sampling/stimprograms/surface11");
        c.compile_from_rewrited_stim_string(stim_str);
    } catch(const std::exception& e){
        std::cerr<<e.what()<<'\n';
    }
    auto t1 = clock::now();                               // stop section‑1
    auto compile_us = std::chrono::duration_cast<microsec>(t1 - t0).count();
    std::cout << "[Time to read, compile circuit from STIM file:] " << compile_us / 1'000.0 << "ms\n";


    t0 = clock::now();       
    QEPG::QEPG graph(c,c.get_num_detector(),c.get_num_noise());
    graph.backward_graph_construction();
    t1 = clock::now();                               // stop section‑1
    compile_us = std::chrono::duration_cast<microsec>(t1 - t0).count();
    std::cout << "[Time to construct QEPG:] " << compile_us / 1'000.0 << "ms\n";


    size_t qubitnum=c.get_num_qubit();

    SAMPLE::sampler sampler(qubitnum);


    //----------------------------------------------------
    // 1 . Compile circuit + build graph
    //----------------------------------------------------
    t0 = clock::now();                               // start timer

    std::vector<QEPG::Row> samplecontainer;
    sampler.generate_many_output_samples(graph,samplecontainer,10,1000000);

    t1 = clock::now();                               // stop section‑1
    compile_us = std::chrono::duration_cast<microsec>(t1 - t0).count();
    std::cout << "[Time to generate samples:] " << compile_us / 1'000.0 << "ms\n";

    // size_t index=0;
    // for(QEPG::Row parityresult: samplecontainer){
    //     std::cout<<index<<":";
    //     QEPG::print_bit_row(parityresult);
    //     index++;
    // }

}




/*
Generate samples, and store output to file path
*/
void store_outputt_to_file(std::string filepath){

    clifford::cliffordcircuit c;
    try{
        std::string stim_str=read_file_to_string("C:/Users/yezhu/OneDrive/Documents/GitHub/Sampling/stimprograms/surface3");
        c.compile_from_rewrited_stim_string(stim_str);
    } catch(const std::exception& e){
        std::cerr<<e.what()<<'\n';
    }

    QEPG::QEPG graph(c,c.get_num_detector(),c.get_num_noise());
    graph.backward_graph_construction();
    size_t qubitnum=c.get_num_qubit();
    SAMPLE::sampler sampler(qubitnum);


    std::vector<QEPG::Row> samplecontainer;
    sampler.generate_many_output_samples(graph,samplecontainer,2,100);


    /*--------------------------------------------------------------------
      1.  Create / overwrite the output file.
          std::ofstream flushes and closes automatically on destruction.
     -------------------------------------------------------------------*/
     std::ofstream ofs(filepath, std::ios::out | std::ios::trunc);
     if (!ofs) {
         std::cerr << "Cannot open \"" << filepath << "\" for writing.\n";
         return;
     }

    for(QEPG::Row parityresult: samplecontainer){
        QEPG::print_bit_row(parityresult, ofs);   // overload that accepts std::ostream&
    }

    /*--------------------------------------------------------------------
      3.  Check for write errors (optional but helpful).
    -------------------------------------------------------------------*/
     if (!ofs) {
        std::cerr << "I/O error while writing \"" << filepath << "\"\n";
    }

}




int main()
{
    //test_matrix_multiplication();

    //test_bitset_pop_count();
    // clifford::cliffordcircuit* c=new clifford::cliffordcircuit(3);

    // try{
    //     std::string stim_str=read_file_to_string("C:/Users/yezhu/OneDrive/Documents/GitHub/Sampling/stimprograms/surface3");
    //     c->compile_from_rewrited_stim_string(stim_str);
    // } catch(const std::exception& e){
    //     std::cerr<<e.what()<<'\n';
    // }

    // c->print_circuit();
    // delete c;

    
    // std::size_t rows = 4, cols = 8;
    // std::vector<boost::dynamic_bitset<>> M(rows, boost::dynamic_bitset<>(cols));

    // // set a few bits
    // M[0].set(0).set(3);
    // M[1].set(4);
    // M[2].set(1).set(2).set(7);
    // M[3].flip();               // row of all 1s

    // print_bit_matrix(M);       // default '0'/'1'
    // std::cout << '\n';


    //c.print_circuit();

    // QEPG::QEPG graph(c,c.get_num_detector(),c.get_num_noise());
    


    // graph.backward_graph_construction();
    // //graph.print_detectorMatrix('0','1');

    // size_t qubitnum=c.get_num_qubit();

    // SAMPLE::sampler sampler(qubitnum);

   
    // std::vector<SAMPLE::singlePauli> result=sampler.generate_sample_Floyd(50);
    // // for(SAMPLE::singlePauli sample:result){
    // //     std::cout<<"("<<sample.qindex<<","<<sample.type<<") ";
    // // }
    // // std::cout<<"\n";
    // QEPG::Row parityresult=sampler.calculate_parity_output_from_one_sample(graph,result);
    // // std::cout<<"Sample parity outcome:"<<"\n";
    // QEPG::print_bit_row(parityresult);


    // std::vector<SAMPLE::singlePauli> result=sampler.generate_sample_Floyd(10);

    // for(SAMPLE::singlePauli sample:result){
    //     std::cout<<"("<<sample.qindex<<","<<sample.type<<") ";
    // }
    // std::cout<<"\n";

    //benchmark_surface_million_sample();
    store_outputt_to_file("output.txt");
}




