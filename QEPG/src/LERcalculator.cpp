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




inline py::array_t<bool> bitset_rows_to_numpy(const std::vector<QEPG::Row>& rows)
{
    using bitset_t  = QEPG::Row;
    using block_t   = bitset_t::block_type;           // usually uint64_t

    const std::size_t n_rows = rows.size();
    const std::size_t n_cols = n_rows ? rows.front().size() : 0;
    if(n_cols==0)
        return py::array_t<bool>({n_rows, n_cols});

    // create a NumPy array of shape (n_rows, n_cols), dtype=bool
    py::array_t<bool> out({n_rows, n_cols});
    // auto buf = out.mutable_unchecked<2>();     // raw, bounds-checked only in debug

    // for (std::size_t r = 0; r < n_rows; ++r) {
    //     const auto& bits = rows[r];            // the current bitset
    //     for (std::size_t c = 0; c < n_cols; ++c) {
    //         buf(r, c) = bits[c];               // copy bit by bit
    //     }
    // }
    // return out;                                // pybind11 returns the array by value
    //Raw buffer pointer(One row=n_cols bytes)
    auto req=out.request();
    auto* base=static_cast<std::uint8_t*>(req.ptr);
    const auto row_stride=n_cols;

    //release the GIL so other python threads can run
    py::gil_scoped_release release;


    //How many bits stored in this block, either 32 or 64
    constexpr std::size_t WORD_BITS = std::numeric_limits<block_t>::digits;



    //-----parallel over rows-------------------------------
    #pragma omp parallel for schedule(static)
    for(long long r=0; r < static_cast<long long>(n_rows);++r)
    {
        const QEPG::Row& bits=rows[static_cast<std::size_t>(r)];
        const std::size_t n_blk=bits.num_blocks();

        /* --- thread-local scratch buffer (namespace scope ==> OK in MSVC) */
        static thread_local std::vector<block_t> tl_buf;
        tl_buf.resize(n_blk);                               // realloc only if needed
        boost::to_block_range(bits, tl_buf.begin());        // fill the buffer

        /*  2. Unpack into the Numpy row (64 bits -> 64 bytes)*/
        std::uint8_t* dst=base+r*row_stride;

        for(std::size_t b=0; b+1 <n_blk;++b){
            std::uint64_t word=static_cast<std::uint64_t>(tl_buf[b]);
            for(int k=0;k<WORD_BITS;++k,word>>=1)
                dst[b*WORD_BITS+k]=static_cast<std::uint8_t>(word&1);
        }

        const std::size_t rem_bits=n_cols&(WORD_BITS-1);
        if(rem_bits){
            std::uint64_t word = static_cast<std::uint64_t>(tl_buf[n_blk-1]);
            for(std::size_t k=0;k<rem_bits;++k, word>>=1)
                dst[(n_blk-1)*WORD_BITS + k]=static_cast<std::uint8_t>(word & 1);
        }
    }
    return out;
}





inline void convert_bitset_row_to_boolean_separate_obs(std::vector<std::vector<bool>>& result,std::vector<bool>& obsresult,const std::vector<QEPG::Row>& samplecontainer){
        result.reserve(samplecontainer.size()); // Reserve space
        obsresult.reserve(samplecontainer.size());
        // Convert each boost::dynamic_bitset<> to std::vector<bool>
        for (const auto& bitset_row : samplecontainer) {
            std::vector<bool> bool_row(bitset_row.size()-1);
            for (size_t i = 0; i < bitset_row.size()-1; ++i) {
                bool_row[i] = bitset_row[i]; // Access individual bits
            }
            result.push_back(bool_row);
            obsresult.push_back(bitset_row[bitset_row.size()-1]);
        }
}


/*
Return a three 
*/
inline void convert_bitset_row_to_boolean_separate_obs_numpy(py::array_t<bool>& detectionresult,py::array_t<bool>& obsresult,const size_t& begin_index,const std::vector<QEPG::Row>& samplecontainer){
        auto detector_buf = detectionresult.mutable_unchecked<2>();     // raw, bounds-checked only in debug
        auto observable_buf = obsresult.mutable_unchecked<1>(); 
        auto row_index=0;
        for (const auto& bitset_row : samplecontainer) {
            std::vector<bool> bool_row(bitset_row.size()-1);
            for (size_t i = 0; i < bitset_row.size()-1; ++i) {
                detector_buf(begin_index+row_index,i)= bitset_row[i]; 
            }
            observable_buf(begin_index+row_index)=bitset_row[bitset_row.size()-1];
            row_index++;
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
    t0 = clock::now();                               // start timer
    convert_bitset_row_to_boolean(result,samplecontainer);
    t1 = clock::now();                               // stop section‑1
    compile_us = std::chrono::duration_cast<microsec>(t1 - t0).count();
    std::cout << "[Time to convert it to bit vector:] " << compile_us / 1'000.0 << "ms\n";


    // for(QEPG::Row parityresult: samplecontainer){
    //     QEPG::print_bit_row(parityresult);
    // }
    return std::move(result);
}




py::array_t<bool> return_samples_numpy(const std::string& prog_str,size_t weight, size_t shots){
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


    py::array_t<bool>  result;
    t0 = clock::now();                               // start timer
    result=bitset_rows_to_numpy(samplecontainer);
    t1 = clock::now();                               // stop section‑1
    compile_us = std::chrono::duration_cast<microsec>(t1 - t0).count();
    std::cout << "[Time to convert it to bit vector:] " << compile_us / 1'000.0 << "ms\n";


    // for(QEPG::Row parityresult: samplecontainer){
    //     QEPG::print_bit_row(parityresult);
    // }
    return std::move(result);    
}






 std::vector<std::vector<bool>> return_all_samples_with_fixed_weights(const std::string& prog_str,const size_t& weight){
    clifford::cliffordcircuit c;
    c.compile_from_rewrited_stim_string(prog_str);

    QEPG::QEPG graph(c,c.get_num_detector(),c.get_num_noise());
    graph.backward_graph_construction();


    SAMPLE::sampler sampler(c.get_num_noise());

    std::vector<QEPG::Row> samplecontainer;

    using clock     = std::chrono::steady_clock;          // monotonic, good for benchmarking
    using microsec  = std::chrono::microseconds;
    auto t0 = clock::now();       
    
    // start timer
    sampler.generate_all_samples_with_fixed_weight(graph,samplecontainer,weight);



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


std::pair<std::vector<std::vector<std::pair<int,int>>> ,std::vector<std::vector<bool>>> 
return_samples_with_noise_vector(const std::string & prog_str,size_t weight, size_t shots){
    clifford::cliffordcircuit c;
    c.compile_from_rewrited_stim_string(prog_str);

    QEPG::QEPG graph(c,c.get_num_detector(),c.get_num_noise());
    graph.backward_graph_construction();


    SAMPLE::sampler sampler(c.get_num_noise());

    std::vector<QEPG::Row> samplecontainer;
    std::vector<std::vector<SAMPLE::singlePauli>> noisecontainer;


    using clock     = std::chrono::steady_clock;          // monotonic, good for benchmarking
    using microsec  = std::chrono::microseconds;
    auto t0 = clock::now();                               // start timer
    //sampler.generate_many_output_samples(graph,samplecontainer,weight,shots);
    sampler.generate_many_output_samples_with_noise_vector(graph,noisecontainer,samplecontainer,weight,shots);


    auto t1 = clock::now();                               // stop section‑1
    auto compile_us = std::chrono::duration_cast<microsec>(t1 - t0).count();
    std::cout << "[Time to generate these samples:] " << compile_us / 1'000.0 << "ms\n";


    std::vector<std::vector<bool>> sampleresult;
    convert_bitset_row_to_boolean(sampleresult,samplecontainer);

    std::vector<std::vector<std::pair<int,int>>> noisegenerated;
    noisegenerated.reserve(shots);
    for(std::vector<SAMPLE::singlePauli> tmpnoisevector: noisecontainer){
        std::vector<std::pair<int,int>> outputnoisevector;
        for(SAMPLE::singlePauli tmpnoise: tmpnoisevector){
              outputnoisevector.push_back(std::pair<int,int>{tmpnoise.qindex,tmpnoise.type});
        }
        noisegenerated.push_back(outputnoisevector);
    }

    return std::pair<std::vector<std::vector<std::pair<int,int>>> ,std::vector<std::vector<bool>>>{std::move(noisegenerated),std::move(sampleresult)};
}



std::vector<std::vector<std::vector<bool>>> return_samples_many_weights(const std::string& prog_str,const std::vector<size_t>& weight, const std::vector<size_t>& shots){
    clifford::cliffordcircuit c;
    c.compile_from_rewrited_stim_string(prog_str);

    QEPG::QEPG graph(c,c.get_num_detector(),c.get_num_noise());
    graph.backward_graph_construction();

    SAMPLE::sampler sampler(c.get_num_noise());

    std::vector<QEPG::Row> samplecontainer;
    std::vector<std::vector<bool>> tmpresult;
    tmpresult.reserve(weight.size());
    std::vector<std::vector<std::vector<bool>>> result;
    result.reserve(weight.size());

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
    return std::move(result);
}




std::vector<py::array_t<bool>> return_samples_many_weights_numpy(const std::string& prog_str,const std::vector<size_t>& weight, const std::vector<size_t>& shots){
    clifford::cliffordcircuit c;
    c.compile_from_rewrited_stim_string(prog_str);

    QEPG::QEPG graph(c,c.get_num_detector(),c.get_num_noise());
    graph.backward_graph_construction();

    SAMPLE::sampler sampler(c.get_num_noise());

    std::vector<QEPG::Row> samplecontainer;
    std::vector<py::array_t<bool>> result;
    result.reserve(weight.size());

    for(size_t i=0;i<weight.size();++i){
        std::cout<<"Weight="<<weight[i]<<"\n";
        samplecontainer.clear();
        py::array_t<bool> tmpresult;
        sampler.generate_many_output_samples(graph,samplecontainer,weight[i],shots[i]);
        tmpresult=bitset_rows_to_numpy(samplecontainer);
        result.emplace_back(std::move(tmpresult));
    }
    // for(QEPG::Row parityresult: samplecontainer){
    //     QEPG::print_bit_row(parityresult);
    // }
    return std::move(result);
}










 std::pair<py::array_t<bool>,py::array_t<bool>> return_samples_many_weights_separate_obs(const std::string& prog_str,const std::vector<size_t>& weight, const std::vector<size_t>& shots){
    clifford::cliffordcircuit c;
    c.compile_from_rewrited_stim_string(prog_str);

    QEPG::QEPG graph(c,c.get_num_detector(),c.get_num_noise());
    graph.backward_graph_construction();

    SAMPLE::sampler sampler(c.get_num_noise());


    std::vector<QEPG::Row> samplecontainer;

    size_t shot_sum=0;
    for(int i=0;i<weight.size();i++){
        shot_sum+=shots[i];
    }

    py::array_t<bool> detectorresult({shot_sum,c.get_num_detector()});
    py::array_t<bool> obsresult(shot_sum);

    auto begin_index=0;
    for(size_t i=0;i<weight.size();++i){
        std::cout<<"Weight="<<weight[i]<<"\n";
        samplecontainer.clear();
        sampler.generate_many_output_samples(graph,samplecontainer,weight[i],shots[i]);
        convert_bitset_row_to_boolean_separate_obs_numpy(detectorresult,obsresult,begin_index,samplecontainer);
        begin_index+=shots[i];
    }
    // for(QEPG::Row parityresult: samplecontainer){
    //     QEPG::print_bit_row(parityresult);
    // }
    return std::pair<py::array_t<bool>,py::array_t<bool>>{std::move(detectorresult),std::move(obsresult)};
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