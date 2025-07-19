#include <iostream>
#include <vector>
#include <iomanip>



#include "QEPG/src/clifford.hpp"
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <chrono>
#include <Eigen/Sparse>
//  g++ -O3 -std=c++20 -I/path/to/eigen  main.cpp


std::string read_file_to_string(const std::string& path)
{
    std::ifstream in(path, std::ios::in | std::ios::binary);
    if(!in) throw std::runtime_error("Cannot open "+ path);

    std::ostringstream buffer;
    buffer<< in.rdbuf();
    return buffer.str();
};



#include "QEPG/src/QEPG.hpp"
#include "QEPG/src/sampler.hpp"




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
        std::string stim_str=read_file_to_string("C:/Users/username/OneDrive/Documents/GitHub/Sampling/stimprograms/surface11");
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
        std::string stim_str=read_file_to_string("C:/Users/username/OneDrive/Documents/GitHub/Sampling/stimprograms/surface9");
        c.compile_from_rewrited_stim_string(stim_str);
    } catch(const std::exception& e){
        std::cerr<<e.what()<<'\n';
    }

    QEPG::QEPG graph(c,c.get_num_detector(),c.get_num_noise());
    graph.backward_graph_construction();
    size_t qubitnum=c.get_num_qubit();
    SAMPLE::sampler sampler(qubitnum);


    std::vector<QEPG::Row> samplecontainer;
    sampler.generate_many_output_samples(graph,samplecontainer,8,1000000);


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










/*
Benchmard the running time of generating 1 million surface code
sample
*/
void print_surface_output(){



    clifford::cliffordcircuit c;
    try{
        std::string stim_str=read_file_to_string("C:/Users/username/OneDrive/Documents/GitHub/Sampling/stimprograms/surface11");
        c.compile_from_rewrited_stim_string(stim_str);
    } catch(const std::exception& e){
        std::cerr<<e.what()<<'\n';
    }

    QEPG::QEPG graph(c,c.get_num_detector(),c.get_num_noise());
    graph.backward_graph_construction();

    size_t qubitnum=c.get_num_qubit();

    SAMPLE::sampler sampler(qubitnum);



    std::vector<QEPG::Row> samplecontainer;
    sampler.generate_many_output_samples(graph,samplecontainer,2,1000000);


    for(QEPG::Row parityresult: samplecontainer){
        QEPG::print_bit_row(parityresult);
    }

}






void test_eigen_QEPG(){
 
    clifford::cliffordcircuit c;
    try{
        std::string stim_str=read_file_to_string("C:/Users/username/OneDrive/Documents/GitHub/Sampling/stimprograms/surface15");
        c.compile_from_rewrited_stim_string(stim_str);
    } catch(const std::exception& e){
        std::cerr<<e.what()<<'\n';
    }

    QEPG::QEPG graph(c,c.get_num_detector(),c.get_num_noise());
    graph.backward_graph_construction_Eigen();

    size_t qubitnum=c.get_num_qubit();

}





// using SpMat = Eigen::SparseMatrix<bool, Eigen::RowMajor>;     // NOTE row-major
// using RowIt = Eigen::SparseMatrix<bool, Eigen::RowMajor>::InnerIterator;

// void xor_rows(const SpMat& src, int ia, int ib,
//               SpMat& dst,     int ic)
// {
//     // 1. Reserve enough slots for the result  (union of nnz)
//     std::size_t cap = src.row(ia).nonZeros() + src.row(ib).nonZeros();
//     dst.row(ic).resize(src.cols());          // make sure the row exists
//     dst.row(ic).reserve(cap);                // O(1) inserts

//     RowIt itA(src, ia), itB(src, ib);        // row iterators
//     while (itA || itB)
//     {
//         if (!itB || (itA && itA.index() < itB.index()))
//         {   dst.insertBack(ic, itA.index()) = true; ++itA; }
//         else if (!itA || itB.index() < itA.index())
//         {   dst.insertBack(ic, itB.index()) = true; ++itB; }
//         else                                   // 1⊕1 ⇒ 0  (skip)
//         {   ++itA; ++itB; }
//     }
//     dst.row(ic).finalize();                   // compact & sorted
// }



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


#include "QEPG/src/LERcalculator.hpp";

void print_bool_matrix(const std::vector<std::vector<bool>>& mat)
{
    for (const auto& row : mat) {
        for (bool bit : row) {
            std::cout << (bit ? '1' : '0');   // or std::cout << std::boolalpha << bit;
        }
        std::cout << '\n';
    }
}



 std::vector<std::vector<bool>> return_samples(const std::string& prog_str,size_t weight, size_t shots){
    clifford::cliffordcircuit c;
    c.compile_from_rewrited_stim_string(prog_str);

    QEPG::QEPG graph(c,c.get_num_detector(),c.get_num_noise());
    graph.backward_graph_construction();

    size_t qubitnum=c.get_num_noise();

    SAMPLE::sampler tmpsampler(c.get_num_noise());

    std::vector<QEPG::Row> samplecontainer;

    using clock     = std::chrono::steady_clock;          // monotonic, good for benchmarking
    using microsec  = std::chrono::microseconds;
    auto t0 = clock::now();                               // start timer
    tmpsampler.generate_many_output_samples(graph,samplecontainer,weight,shots);

    std::cout<<"Samples "<<"\n";
    size_t failure=0;
    for(QEPG::Row parityresult:samplecontainer){
        QEPG::print_bit_row(parityresult);
        if(parityresult[1]){
            failure++;
        }
    }
    std::cout<<failure<<"\n";


    auto t1 = clock::now();                               // stop section‑1
    auto compile_us = std::chrono::duration_cast<microsec>(t1 - t0).count();
    std::cout << "[Time to generate these samples    :] " << compile_us / 1'000.0 << "ms\n";


    std::vector<std::vector<bool>> result;
    convert_bitset_row_to_boolean(result,samplecontainer);

    // for(QEPG::Row parityresult: samplecontainer){
    //     QEPG::print_bit_row(parityresult);
    // }
    return result;
}


// ------------------------------------------------------------------
// enumerate_all_vectors ---------------------------------------------------------
//   length  – total number of bits (n)
//   weight  – how many 1-bits (k)
//   visit   – callback invoked with each vector
//             (avoids storing everything in memory)
// ------------------------------------------------------------------
void enumerate_all_vectors(
        std::size_t length,
        std::size_t weight,
        const std::function<void(const std::vector<int>&)>& visit)
{
    if (weight > length) return;          // nothing to generate

    std::vector<int> vec(length, 0);      // current 0/1 vector
    std::vector<std::size_t> idx(weight); // indices of 1-bits

    // initialise first combination 0,1,2,…,k–1
    for (std::size_t i = 0; i < weight; ++i) idx[i] = i;

    auto push_vec = [&]{
        std::fill(vec.begin(), vec.end(), 0);
        for (std::size_t i : idx) vec[i] = 1;
        visit(vec);
    };

    if (weight == 0) {                    // special-case k = 0
        visit(vec);
        return;
    }

    // iterate through all k-combinations of {0,…,n-1}
    while (true) {
        push_vec();                       // emit current vector

        // generate next combination (EKM algorithm)
        std::size_t i = weight;
        while (i-- > 0 && idx[i] == length - weight + i);
        if (i == static_cast<std::size_t>(-1)) break;     // done

        ++idx[i];
        for (std::size_t j = i + 1; j < weight; ++j)
            idx[j] = idx[j - 1] + 1;
    }
}


int main()
{
    //test_matrix_multiplication();

    //test_bitset_pop_count();
    // clifford::cliffordcircuit c(3);

    // std::cout<<"Start compilation!"<<"\n";
    // std::string stim_str;
    // try{
    //     stim_str=read_file_to_string("C:/Users/username/Documents/Sampling/stimprograms/1cnot");
    //     c.compile_from_rewrited_stim_string(stim_str);
    // } catch(const std::exception& e){
    //     std::cerr<<e.what()<<'\n';
    // }

    // QEPG::QEPG graph(c,c.get_num_detector(),c.get_num_noise());
    // graph.backward_graph_construction();

    std::size_t n = 4, k = 2;
    enumerate_all_vectors(n, k, [](const std::vector<int>& v){
        for (int bit : v) std::cout << bit;
        std::cout << '\n';
    });

}




