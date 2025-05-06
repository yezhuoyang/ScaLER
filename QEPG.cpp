#include "QEPG.hpp"
#include <iostream>

namespace QEPG{

/*------------------------------------------ctor----*/

QEPG::QEPG(){

}

QEPG::QEPG(clifford::cliffordcircuit othercircuit, size_t total_detectors, size_t total_noise):
                    circuit(othercircuit),
                    total_detectors_(total_detectors),
                    total_noise_(total_noise),
                    detectorMatrix_(total_detectors,Row(3*total_noise))                   
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
    size_t gate_size=circuit.get_gate_num();

    std::vector<Row> current_x_prop(circuit.get_num_qubit(),Row(3* total_noise_));
    std::vector<Row> current_y_prop(circuit.get_num_qubit(),Row(3* total_noise_));
    std::vector<Row> current_z_prop(circuit.get_num_qubit(),Row(3* total_noise_));



    for(int t=gate_size-1;t>=0;t--){
        const auto& gate=circuit.get_gate(t);


        std::string name=gate.name;
        std::cout<<name<<"\n";
        /*
        *   First case, when the gate is a depolarization noise
        */
        if(name=="DEPOLARIZE1"){
                std::cout<<"Find depolarization noise!";
        }
        /*
        *   When the gate is a measurement
        */
        if(name=="M"){

        }
        /*
        *   When the gate is a reset
        */
        if(name=="R"){


        }
        /*
        *   When the gate is a CNOT
        */
        if(name=="cnot"){


        }

        if(name=="h"){

        }

    }

} 








}