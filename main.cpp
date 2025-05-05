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





int main()
{
    clifford::cliffordcircuit c(3);

    try{
        std::string stim_str=read_file_to_string("C:/Users/yezhu/OneDrive/Documents/GitHub/Sampling/stimprograms/surface3");
        c.compile_from_rewrited_stim_string(stim_str);
    } catch(const std::exception& e){
        std::cerr<<e.what()<<'\n';
    }

    c.print_circuit();


}




