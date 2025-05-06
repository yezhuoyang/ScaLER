#include <iostream>
#include <vector>
#include <iomanip>




using Matrix = std::vector<std::vector<int>>;
using Row=std::vector<int>;
using size_type=std::size_t;


bool dimensions_match(const Matrix& A, const Matrix& B){
    return !A.empty() && !B.empty() && A.front().size()==B.size();
}

void print_matrix(const Matrix& M, const std::string& name= "M"){
    std::cout<<name<<"("<<M.size()<<" x "<<M.front().size()<<"):\n";
    constexpr int width=6;
    for(const Row& r: M){
        for(int x: r) std::cout << std::setw(width)<<x;
        std::cout << '\n';
    }
    std::cout << "\n";
}


Matrix multiply(const Matrix& A, const Matrix& B){

    if(!dimensions_match(A,B))
        throw std::invalid_argument("incompatible dimensions for A * B");


    size_type R=A.size(),
              C=B.front().size(),
              K=B.size();  

    Matrix Cmat(R,Row(C,0));


    for(size_type i=0; i<R;++i)
        for (size_type k=0; k<K; ++k)
            for (size_type j=0; j<C; ++j)
                Cmat[i][j]+=A[i][k]*B[k][j];

    return Cmat;
}


#include "clifford.hpp"
#include <stdexcept>
#include <fstream>
#include <sstream>

std::string read_file_to_string(const std::string& path)
{
    std::ifstream in(path, std::ios::in | std::ios::binary);
    if(!in) throw std::runtime_error("Cannot open "+ path);

    std::ostringstream buffer;
    buffer<< in.rdbuf();
    return buffer.str();
};



#include "QEPG.hpp"
#include "sampler.hpp"

template<class BitRow>
void print_bit_matrix(const std::vector<BitRow>& rows,
                      char zero = '0', char one='1')
{
    if(rows.empty()) return;

    const std::size_t cols= rows.front().size();
    for(const auto& r: rows){
        if(r.size()!=cols){
            std::cerr<<"[print_bit_matrix] row width mismatach\n";
            return;
        }
        for(std::size_t c=0; c<cols;++c){
            std::cout<<(r.test(c)? one: zero);
        }
        std::cout<<"\n";
    }
}

int main()
{

    
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
    
    // clifford::cliffordcircuit c(3);

    // try{
    //     std::string stim_str=read_file_to_string("C:/Users/yezhu/OneDrive/Documents/GitHub/Sampling/stimprograms/simple");
    //     c.compile_from_rewrited_stim_string(stim_str);
    // } catch(const std::exception& e){
    //     std::cerr<<e.what()<<'\n';
    // }

    // c.print_circuit();

    // QEPG::QEPG graph(c,c.get_num_detector(),c.get_num_noise());


    // graph.backward_graph_construction();

    SAMPLE::sampler sampler(100);


    std::vector<SAMPLE::singlePauli> result=sampler.generate_sample_Floyd(10);

    for(SAMPLE::singlePauli sample:result){
        std::cout<<"("<<sample.qindex<<","<<sample.type<<") ";
    }
    std::cout<<"\n";


}




