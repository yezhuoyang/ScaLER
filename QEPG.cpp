#include "QEPG.hpp"



namespace QEPG{

/*------------------------------------------ctor----*/

QEPG::QEPG(){

}

QEPG::QEPG(clifford::cliffordcircuit othercircuit, size_t total_detectors, size_t total_noise):
                    circuit(othercircuit),
                    total_detectors_(total_detectors),
                    total_noise_(total_noise),
                    X_error_(total_detectors,Row(3*total_noise)),
                    Y_error_(total_detectors,Row(3*total_noise)), 
                    Z_error_(total_detectors,Row(3*total_noise))                   
                    {

}


QEPG::~QEPG(){

}



/*--------------------print QEPG graph---------------------------------------*/


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




/*---------------Construction of the QEPG graph-------------------------------*/


void QEPG::backward_graph_construction(){
    
} 








}