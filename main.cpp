#include <iostream>
#include <vector>
#include <iomanip>



#include "QEPG/clifford.hpp"
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
        std::string stim_str=read_file_to_string("C:/Users/yezhu/OneDrive/Documents/GitHub/Sampling/stimprograms/surface9");
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
        std::string stim_str=read_file_to_string("C:/Users/yezhu/OneDrive/Documents/GitHub/Sampling/stimprograms/surface11");
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
        std::string stim_str=read_file_to_string("C:/Users/yezhu/OneDrive/Documents/GitHub/Sampling/stimprograms/surface15");
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





int main()
{
    //test_matrix_multiplication();

    //test_bitset_pop_count();
    clifford::cliffordcircuit c(3);

    std::cout<<"Start compilation!"<<"\n";
    try{
        std::string stim_str=read_file_to_string("C:/Users/yezhu/Documents/Sampling/stimprograms/cnot01");
        c.compile_from_rewrited_stim_string(stim_str);
    } catch(const std::exception& e){
        std::cerr<<e.what()<<'\n';
    }

    std::cout<<"Start printing"<<"\n";
    c.print_circuit();


    QEPG::QEPG graph(c,c.get_num_detector(),c.get_num_noise());
    

    graph.backward_graph_construction();
    graph.print_detectorMatrix('0','1');

    //size_t qubitnum=c.get_num_qubit();

    //SAMPLE::sampler sampler(qubitnum);


    // std::size_t rows = 4, cols = 8;
    // std::vector<boost::dynamic_bitset<>> M(rows, boost::dynamic_bitset<>(cols));

    // // set a few bits
    // M[0].set(0).set(3);
    // M[1].set(4);
    // M[2].set(1).set(2).set(7);
    // M[3].flip();               // row of all 1s

    // print_bit_matrix(M);       // default '0'/'1'
    // std::cout << '\n';




   
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
    //store_outputt_to_file("output.txt");
    //print_surface_output();

    // using namespace Eigen;
    // Matrix<bool,Dynamic,Dynamic> A(4,6);
    // A.setZero();
    // A(1,2)=true;
    // A(0,4)=true;
    // A(3,3)=true;

    // RowVectorX<bool> r0=A.row(0);
    // RowVectorX<bool> xor01=(A.row(0).array()^A.row(1).array()).matrix();
    // A.row(2)=xor01;


    // RowVectorX<bool> xor23=(A.col(2).array()^A.col(3).array()).matrix();
    // A.col(5)=xor23;

    // IOFormat cleanFmt(StreamPrecision,DontAlignCols," ","\n","","");


    // std::cout<<"A.row(2)=\n"<<(A.row(2).array().matrix()).cast<int>().format(cleanFmt)<<"\n\n";


    // std::cout<<"A=\n"<<A.cast<int>().format(cleanFmt)<<"\n\n";

    //test_eigen_QEPG();
    //store_outputt_to_file("output");

    // constexpr int    N       = 10'000000000;
    // constexpr double density = 0.0005;     // 0.05 %  ⇒ ~5 M nnz  (≈20 MB)
    // std::mt19937 gen(42);

    // auto A = randomSparse(N, density, gen);
    // auto B = randomSparse(N, density, gen);

    // std::cout << "A nnz = " << A.nonZeros()
    //           << ", B nnz = " << B.nonZeros() << '\n';

    // auto t0 = std::chrono::steady_clock::now();
    // SpMat C = A * B;                                    // sparse × sparse
    // auto t1 = std::chrono::steady_clock::now();

    // std::cout << "C nnz = " << C.nonZeros()
    //           << ", elapsed " << std::chrono::duration<double>(t1-t0).count()
    //           << " s\n";
    


    // using SpMat = Eigen::SparseMatrix<bool,Eigen::RowMajor>;


    // SpMat M(3,8);                     // row-major
    // M.insert(0,1)=true;  M.insert(0,4)=true;
    // M.insert(1,4)=true;  M.insert(1,7)=true;

    // // XOR row-0 and row-1 into row-2
    // Eigen::SparseVector<bool> tmp(M.cols());
    // tmp.reserve(M.row(0).nonZeros() + M.row(1).nonZeros());
    // SpMat::InnerIterator a(M,0), b(M,1);
    // while(a||b){
    //     if(!b || (a && a.index()<b.index()))   { tmp.insertBack(a.index()) = true; ++a; }
    //     else if(!a || b.index()<a.index())     { tmp.insertBack(b.index()) = true; ++b; }
    //     else                                   { ++a; ++b; }              // 1⊕1=0
    // }
    // tmp.finalize();
    // M.row(2) = tmp;                            // overwrite row-2

    // std::cout << "row-2 indices: ";
    // for(SpMat::InnerIterator it(M,2); it; ++it) std::cout<<it.index()<<' ';
    // std::cout << '\n';   // → 1 7
    // test_eigen_QEPG();
}




