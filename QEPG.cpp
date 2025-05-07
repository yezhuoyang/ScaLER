#include "QEPG.hpp"
#include <iostream>

namespace QEPG{

/*------------------------------------------ctor----*/

QEPG::QEPG(){

}

QEPG::QEPG(clifford::cliffordcircuit othercircuit, size_t total_detectors, size_t total_noise):
                    circuit_(othercircuit),
                    total_detectors_(total_detectors),
                    total_noise_(total_noise),
                    detectorMatrix_(othercircuit.get_num_meas(),Row(3*total_noise))                   
                    {

}


QEPG::~QEPG(){

}



/*--------------------print QEPG graph---------------------------------------*/



void QEPG::print_detectorMatrix(char zero, char one) const{

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
    std::cout<<"-----------ParitygroupMatrix:----------------------\n";
    for(const auto& row:parityPropMatrix_){
        for(std::size_t c=0;c<row.size();++c){
            std::cout<<(row.test(c)? one:zero);
        }
        std::cout<<"\n";
    }
    std::cout<<"-----------ParitygroupMatrixTranspose:----------------------\n";
    for(const auto& row:parityPropMatrixTranspose_){
        for(std::size_t c=0;c<row.size();++c){
            std::cout<<(row.test(c)? one:zero);
        }
        std::cout<<"\n";
    }        
}




/*---------------Construction of the QEPG graph-------------------------------*/


void QEPG::backward_graph_construction(){
    size_t gate_size=circuit_.get_gate_num();
    std::vector<Row> current_x_prop(circuit_.get_num_qubit(),Row(3* total_noise_));
    std::vector<Row> current_y_prop(circuit_.get_num_qubit(),Row(3* total_noise_));
    std::vector<Row> current_z_prop(circuit_.get_num_qubit(),Row(3* total_noise_));

    size_t total_meas=circuit_.get_num_meas();
    size_t current_meas_index=circuit_.get_num_meas()-1;
    size_t current_noise_index=total_noise_-1;


    for(int t=gate_size-1;t>=0;t--){
        const auto& gate=circuit_.get_gate(t);
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
    transpose_matrix(detectorMatrix_,detectorMatrixTranspose_);
    compute_parityPropMatrix();
} 




inline void transpose_matrix(const std::vector<Row>& mat,std::vector<Row>& matTrans){
    const std::size_t n_rows=mat.size();
    const std::size_t n_cols=n_rows ? mat[0].size():0;
    matTrans.assign(n_cols, Row(n_rows));  
    for(std::size_t r=0; r<n_rows; ++r){
        const Row& src= mat[r];
        for(std::size_t c=src.find_first();c!=Row::npos; c=src.find_next(c)){
            matTrans[c].set(r);
        }
    }    
}





const std::vector<Row>& QEPG::get_detectorMatrix() const noexcept{
    return detectorMatrix_;
} 

const std::vector<Row>& QEPG::get_dectorMatrixTrans() const noexcept{
    return detectorMatrixTranspose_;
}

const std::vector<Row>& QEPG::get_parityPropMatrix() const noexcept{
    return parityPropMatrix_;
}







/*
Return the matrix multiplication result of two bitset matrix on Field F2
*/
std::vector<Row> bitset_matrix_multiplication(const std::vector<Row>& mat1,const std::vector<Row>& mat2){
    const size_t row1=mat1.size();
    const size_t col1=row1? mat1[0].size():0;
    const size_t row2=mat2.size();
    const size_t col2=row1? mat2[0].size():0;    
    std::vector<Row> result(row1,Row(col2));
    std::vector<Row> mat2transpose;
    transpose_matrix(mat2,mat2transpose);
    for(size_t i=0;i<row1;i++){
        for(size_t j=0;j<col2;j++){
            result[i][j]=and_popcount(mat1[i],mat2transpose[j])%2? true:false;
        }
    }
    return result;
}


/*
Now we know the propagation of Pauli error to all measurements, we still need to calculate the 
propagation of all pauli error to all detector measurement result
*/
void QEPG::compute_parityPropMatrix(){
    const std::vector<clifford::paritygroup>& detector_parity_group=circuit_.get_detector_parity_group();
    const clifford::paritygroup& observable_group=circuit_.get_observable_parity_group();
    const size_t row_size=detector_parity_group.size()+1;
    const size_t col_size=circuit_.get_num_meas();

    std::vector<Row> paritygroupMatrix(row_size,Row(col_size));

    for(size_t i=0; i<detector_parity_group.size();i++){
        for(size_t index: detector_parity_group[i].indexlist){
            paritygroupMatrix[i][index]=true;
        }
    }
    for(size_t index: observable_group.indexlist){
        paritygroupMatrix[detector_parity_group.size()][index]=true;
    }

    parityPropMatrix_=bitset_matrix_multiplication(paritygroupMatrix,detectorMatrix_);
    transpose_matrix(parityPropMatrix_,parityPropMatrixTranspose_);
}

}