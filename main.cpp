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

// Ensure the QEPG.hpp file defines the 'clifford' namespace or class
// and includes the necessary declarations for 'cliffordcircuit'.



int main()
{
    clifford::cliffordcircuit c(3);
    c.set_error_rate(1e-3);
    c.add_hadamard(0);
    c.add_cnot(0, 1);
    c.add_XError(2);
    c.print_circuit();

    const clifford::Gate& g=c.get_gate(2);

    std::cout<<g.name<<"\n";

}




