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
                    detectorMatrix_(othercircuit.get_num_meas(),Row(3*total_noise))                   
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


void QEPG::print_detectorMatrix(char zero, char one){

    std::cout<<"-----------detectorMatrix:----------------------\n";
    for(const auto& row:detectorMatrix_){
        for(std::size_t c=0;c<row.size();++c){
            std::cout<<(row.test(c)? one:zero);
        }
        std::cout<<"\n";
    }
    std::cout<<"-----------detectorMatrix(Transpose):----------------------\n";
    for(const auto& row:detectorMatrixTranspose_){
        for(std::size_t c=0;c<row.size();++c){
            std::cout<<(row.test(c)? one:zero);
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

    size_t total_meas=circuit.get_num_meas();
    size_t current_meas_index=circuit.get_num_meas()-1;
    size_t current_noise_index=total_noise_-1;


    for(int t=gate_size-1;t>=0;t--){
        const auto& gate=circuit.get_gate(t);
        std::string name=gate.name;
        /*
        *   First case, when the gate is a depolarization noise
        */
        if(name=="DEPOLARIZE1"){
                size_t qindex=gate.qubits[0];
                for(size_t j=0;j<total_meas;j++){
                        detectorMatrix_[j].set(current_noise_index,current_x_prop[qindex].test(j));
                        detectorMatrix_[j].set(total_noise_+current_noise_index,current_y_prop[qindex].test(j));
                        detectorMatrix_[j].set(total_noise_*2+current_noise_index,current_z_prop[qindex].test(j));        
                }
                current_noise_index--;
                continue;
        }
        /*
        *   When the gate is a measurement
        */
        if(name=="M"){
            size_t qindex=gate.qubits[0];
            current_x_prop[qindex].set(current_meas_index);
            current_y_prop[qindex].set(current_meas_index);
            current_meas_index--;
            continue;
        }
        /*
        *   When the gate is a reset
        */
        if(name=="R"){
            size_t qindex=gate.qubits[0];
            for(size_t j=0;j<total_meas;j++){
                    current_x_prop[qindex].set(j,false);
                    current_y_prop[qindex].set(j,false);   
                    current_z_prop[qindex].set(j,false);       
            }
        }
        /*
        *   When the gate is a CNOT
        */
        if(name=="cnot"){
            size_t qcontrol=gate.qubits[0];           
            size_t qtarget=gate.qubits[1];
            current_x_prop[qcontrol]^=current_x_prop[qtarget];
            current_z_prop[qtarget]^=current_z_prop[qcontrol];        
            current_y_prop[qcontrol]^=current_x_prop[qtarget];           
            current_y_prop[qtarget]^=current_z_prop[qcontrol];               
            continue;
        }

        if(name=="h"){
            size_t qindex=gate.qubits[0];
            current_x_prop[qindex].swap(current_z_prop[qindex]);   // fast, no copy
        }
    }

    /*
    Compute the transpose of detectorMatrix for future calculation
    */
    const std::size_t n_rows=detectorMatrix_.size();
    const std::size_t n_cols=n_rows ? detectorMatrix_[0].size():0;
    detectorMatrixTranspose_.assign(n_cols, Row(n_rows));
    
    for(std::size_t r=0; r<n_rows; ++r){
        const Row& src= detectorMatrix_[r];
        for(std::size_t c=src.find_first();c!=Row::npos; c=src.find_next(c)){
            detectorMatrixTranspose_[c].set(r);
        }
    }


} 








}