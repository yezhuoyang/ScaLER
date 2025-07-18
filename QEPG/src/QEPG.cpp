#include "QEPG.hpp"
#include <iostream>
#include <chrono>
#include <algorithm>      // <algorithm> also works but <ranges> is canonical

namespace QEPG{

/*------------------------------------------ctor----*/

QEPG::QEPG(){

}

QEPG::QEPG(clifford::cliffordcircuit othercircuit, size_t total_detectors, size_t total_noise):
                    circuit_(othercircuit),
                    total_detectors_(total_detectors),
                    total_noise_(total_noise),
                    parityPropMatrixTranspose_(3* total_noise_,Row(total_detectors+1))             
                    {
                        
}


QEPG::~QEPG(){

}



const size_t& QEPG::get_total_noise() const noexcept{
    return total_noise_;
}

const size_t& QEPG::get_total_detector() const noexcept{
    return total_detectors_;
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
    
    // std::cout<<"-----------ParitygroupMatrix:----------------------\n";
    // for(const auto& row:parityPropMatrix_){
    //     for(std::size_t c=0;c<row.size();++c){
    //         std::cout<<(row.test(c)? one:zero);
    //     }
    //     std::cout<<"\n";
    // }
    
    std::cout<<"-----------ParitygroupMatrixTranspose:----------------------\n";
    for(const auto& row:parityPropMatrixTranspose_){
        for(std::size_t c=0;c<row.size();++c){
            std::cout<<(row.test(c)? one:zero);
        }
        std::cout<<"\n";
    }        
}




/*---------------Construction of the QEPG graph-------------------------------*/





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

void QEPG::backward_graph_construction(){

    // using clock     = std::chrono::steady_clock;          // monotonic, good for benchmarking
    // using microsec  = std::chrono::microseconds;
    // auto t0 = clock::now();                               // start timer

    size_t gate_size=circuit_.get_gate_num();

    /*
    Store the propagation from pauli noise to qubits
    */
    const size_t total_meas=circuit_.get_num_meas();


    /*
    Directly store the propagation from pauli noise to qubits
    */
    const size_t num_detectors=circuit_.get_detector_parity_group().size();
    
    std::vector<Row> current_x_parity_prop(circuit_.get_num_qubit(),Row(num_detectors+1));
    std::vector<Row> current_y_parity_prop(circuit_.get_num_qubit(),Row(num_detectors+1));
    std::vector<Row> current_z_parity_prop(circuit_.get_num_qubit(),Row(num_detectors+1));


    size_t current_meas_index=circuit_.get_num_meas()-1;
    size_t current_noise_index=total_noise_-1;

    const clifford::paritygroup& observable=circuit_.get_observable_parity_group();

    for(int t=gate_size-1;t>=0;t--){

        const auto& gate=circuit_.get_gate(t);
        std::string name=gate.name;

        /*
        *   First case, when the gate is a depolarization noise
        */
        if(name=="DEPOLARIZE1"){
                size_t qindex=gate.qubits[0];
                //TODO: Is there any way to avoid the nested for loop? 
                // detectorMatrixTranspose_[current_noise_index]=current_x_prop[qindex];
                // detectorMatrixTranspose_[total_noise_+current_noise_index]=current_y_prop[qindex];
                // detectorMatrixTranspose_[total_noise_*2+current_noise_index]=current_z_prop[qindex];
                
                // for(size_t j=0;j<total_meas;j++){
                //         detectorMatrix_[j].set(current_noise_index,current_x_prop[qindex].test(j));
                //         detectorMatrix_[j].set(total_noise_+current_noise_index,current_y_prop[qindex].test(j));
                //         detectorMatrix_[j].set(total_noise_*2+current_noise_index,current_z_prop[qindex].test(j));        
                // }
                
                /*
                Uptill now, the fate of this noise is determined
                So in priciple, we can update the parity propagation
                */
                parityPropMatrixTranspose_[current_noise_index]=current_x_parity_prop[qindex];   //current_x_prop(circuit_.get_num_qubit(),Row(3* total_noise_))
                parityPropMatrixTranspose_[total_noise_+current_noise_index]=current_y_parity_prop[qindex];
                parityPropMatrixTranspose_[total_noise_*2+current_noise_index]=current_z_parity_prop[qindex];      
                
                current_noise_index--;
                continue;
        }
        /*
        *   When the gate is a measurement
        */
        if(name=="M"){
            size_t qindex=gate.qubits[0];
            // current_x_prop[qindex].set(current_meas_index);
            // current_y_prop[qindex].set(current_meas_index);

            /*
            Update all affected parity/detector measurement
            */
            const clifford::parityIndexgroup& tmpmeasuregroup=circuit_.get_measure_to_parity_index(current_meas_index);
            for(size_t parityindex: tmpmeasuregroup.indexlist){
                    current_x_parity_prop[qindex].set(parityindex);
                    current_y_parity_prop[qindex].set(parityindex);
            }
            /*
            This measurement will flip the observable
            */
            if(std::find(observable.indexlist.begin(), observable.indexlist.end(), current_meas_index) != observable.indexlist.end()){
                    current_x_parity_prop[qindex].set(num_detectors);
                    current_y_parity_prop[qindex].set(num_detectors);
            }
            current_meas_index--;
            continue;
        }
        /*
        *   When the gate is a reset
        */
        if(name=="R"){
            size_t qindex=gate.qubits[0];

            current_x_parity_prop[qindex].reset();
            current_y_parity_prop[qindex].reset();  
            current_z_parity_prop[qindex].reset();                      

        }
        /*
        *   When the gate is a CNOT
        */
        if(name=="cnot"){
            size_t qcontrol=gate.qubits[0];           
            size_t qtarget=gate.qubits[1];   
            current_x_parity_prop[qcontrol]^=current_x_parity_prop[qtarget];
            current_z_parity_prop[qtarget]^=current_z_parity_prop[qcontrol];
            current_y_parity_prop[qcontrol]^=current_x_parity_prop[qtarget];
            current_y_parity_prop[qtarget]^=current_z_parity_prop[qcontrol];            
            continue;
        }

        if(name=="h"){
            size_t qindex=gate.qubits[0];
            current_x_parity_prop[qindex].swap(current_z_parity_prop[qindex]);
        }
    }
} 


/// XOR the ia-th row of `dst` and ib-th row of src, store into the ia-th row of `dst`
// dst is a non-const reference as it is being modified.
void xorRow(SpMat &dst, int ia,
    const SpMat &src,       int ib) // src remains const as it is not modified
{
        // 1. Reserve enough space: union of non‑zeros in both source rows
        // We can still use innerVector to get the non-zeros count efficiently
        std::size_t cap = dst.innerVector(ia).nonZeros()
                + src.innerVector(ib).nonZeros();

        SpVec tmp(dst.cols());
        tmp.reserve(cap);

        // 2. Merge the two sorted sparse-row iterators (Θ(nnzA + nnzB))
        // Use SpMat::InnerIterator as shown in your working code
        SpMat::InnerIterator a(dst, ia); // Iterator for the ia-th row of dst
        SpMat::InnerIterator b(src, ib); // Iterator for the ib-th row of src

        while(a || b) { // Continue as long as either iterator is valid
            if (!b || (a && a.index() < b.index())) {
                // Element only in 'a' (dst row) or 'a' is before 'b'
                // This element is part of the result (1 XOR 0 = 1)
                tmp.insertBack(a.index()) = true;
                ++a; // Move to the next element in 'a'
            }
            else if (!a || b.index() < a.index()) {
                // Element only in 'b' (src row) or 'b' is before 'a'
                // This element is part of the result (0 XOR 1 = 1)
                tmp.insertBack(b.index()) = true;
                ++b; // Move to the next element in 'b'
            }
            else {
                // Elements at the same index in both 'a' and 'b'
                // Their XOR is 1 XOR 1 = 0, so we skip both
                ++a; // Move to the next element in 'a'
                ++b; // Move to the next element in 'b'
            }
        }

        // 3. Finalize and assign back to dst’s ia-th row
        tmp.finalize();
        // Assign the temporary sparse vector back to the row.
        // This requires a mutable reference to the row.
        dst.innerVector(ia) = tmp;
}



void QEPG::backward_graph_construction_Eigen(){

    using clock     = std::chrono::steady_clock;          // monotonic, good for benchmarking
    using microsec  = std::chrono::microseconds;
    auto t0 = clock::now();                               // start timer

    detectorMatrixEigen_.setZero();
    size_t gate_size=circuit_.get_gate_num();
    SpMat current_x_prop(circuit_.get_num_qubit(),3*total_noise_);
    SpMat current_y_prop(circuit_.get_num_qubit(),3*total_noise_);
    SpMat current_z_prop(circuit_.get_num_qubit(),3*total_noise_);


    current_x_prop.setZero();
    current_y_prop.setZero();
    current_z_prop.setZero();

    size_t total_meas=circuit_.get_num_meas();
    size_t current_meas_index=circuit_.get_num_meas()-1;
    size_t current_noise_index=total_noise_-1;
    
    SpMat detectorMatrixTransposeEigen_(3*total_noise_,circuit_.get_num_meas());
    detectorMatrixTransposeEigen_.setZero();

    for(int t=gate_size-1;t>=0;t--){

        const auto& gate=circuit_.get_gate(t);
        std::string name=gate.name;
        //std::cout<<t<<" "<<name<<"\n";
        /*
        *   First case, when the gate is a depolarization noise
        */
        if(name=="DEPOLARIZE1"){
                size_t qindex=gate.qubits[0];

                detectorMatrixTransposeEigen_.row(current_noise_index)=current_x_prop.row(qindex);
                detectorMatrixTransposeEigen_.row(total_noise_+current_noise_index)=current_y_prop.row(qindex);
                detectorMatrixTransposeEigen_.row(total_noise_*2+current_noise_index)=current_z_prop.row(qindex);

                current_noise_index--;
                continue;
        }
        /*
        *   When the gate is a measurement
        */
        if(name=="M"){
            size_t qindex=gate.qubits[0];
            current_x_prop.insert(qindex,current_meas_index)=true;
            current_y_prop.insert(qindex,current_meas_index)=true;
            current_meas_index--;
            continue;
        }
        /*
        *   When the gate is a reset
        */
        if(name=="R"){
            size_t qindex=gate.qubits[0];
            Eigen::SparseVector<bool> empty_rowx(current_x_prop.cols()); // nnz = 0
            current_x_prop.row(qindex) = empty_rowx;                          
            Eigen::SparseVector<bool> empty_rowy(current_y_prop.cols()); // nnz = 0
            current_y_prop.row(qindex) = empty_rowy;    
            Eigen::SparseVector<bool> empty_rowz(current_z_prop.cols()); // nnz = 0
            current_z_prop.row(qindex) = empty_rowz;    
        }
        /*
        *   When the gate is a CNOT
        */
        if(name=="cnot"){
            size_t qcontrol=gate.qubits[0];           
            size_t qtarget=gate.qubits[1];
        // XOR control ← control ⊕ target
        
        xorRow(current_x_prop, qcontrol,
            current_x_prop, qtarget);

        // XOR target ← target ⊕ control
        xorRow(current_z_prop, qtarget,
                current_z_prop, qcontrol);

        // update Y-propagations similarly
        xorRow(current_y_prop, qcontrol,
                current_x_prop, qtarget);

        xorRow(current_y_prop, qtarget,
                current_z_prop, qcontrol);
        continue;
        }
        if(name=="h"){
            size_t qindex=gate.qubits[0];
            SpVec rowx=current_x_prop.innerVector(qindex);
            SpVec rowz=current_z_prop.innerVector(qindex);
            rowx.swap(rowz);
        }
    }
}


/*
In this implementation, we directly construct the parity group matrix
by the backward algorithm. This implementation avoid large sparse matrix multiplication.
*/
void QEPG::backward_parity_matrix_construction_Eigen(){

    using clock     = std::chrono::steady_clock;          // monotonic, good for benchmarking
    using microsec  = std::chrono::microseconds;
    auto t0 = clock::now();                               // start timer

    const std::vector<clifford::paritygroup>& detector_parity_group=circuit_.get_detector_parity_group();
    const clifford::paritygroup& observable_group=circuit_.get_observable_parity_group();
    
    const size_t row_size=3*circuit_.get_num_noise();
    const size_t col_size=detector_parity_group.size()+1;

    SpMat parityPropMatrixEigen_(row_size,col_size);


    SpMat current_x_prop(circuit_.get_num_qubit(),3*total_noise_);
    SpMat current_y_prop(circuit_.get_num_qubit(),3*total_noise_);
    SpMat current_z_prop(circuit_.get_num_qubit(),3*total_noise_);


    /*
    Problem now: Current_xyz_prop are essential for resursive algorithm to work, so 
    we cannot directly change the first dimension to be detector_parity_group.size()+1.   
    What can we do then? There must be similar but different induction/resursive rule
    */

    current_x_prop.setZero();
    current_y_prop.setZero();
    current_z_prop.setZero();    
    


    size_t total_meas=circuit_.get_num_meas();
    size_t current_meas_index=circuit_.get_num_meas()-1;
    size_t current_noise_index=total_noise_-1;
    
    SpMat detectorMatrixTransposeEigen_(3*total_noise_,circuit_.get_num_meas());
    detectorMatrixTransposeEigen_.setZero();


    size_t gate_size=circuit_.get_gate_num();
    for(int t=gate_size-1;t>=0;t--){

        const auto& gate=circuit_.get_gate(t);
        std::string name=gate.name;
        //std::cout<<t<<" "<<name<<"\n";
        /*
        *   First case, when the gate is a depolarization noise
        */
        if(name=="DEPOLARIZE1"){
                size_t qindex=gate.qubits[0];

                current_noise_index--;
                continue;
        }
        /*
        *   When the gate is a measurement
        */
        if(name=="M"){
            size_t qindex=gate.qubits[0];
            current_x_prop.insert(qindex,current_meas_index)=true;
            current_y_prop.insert(qindex,current_meas_index)=true;
            current_meas_index--;
            continue;
        }
        /*
        *   When the gate is a reset
        */
        if(name=="R"){
            size_t qindex=gate.qubits[0];
            Eigen::SparseVector<bool> empty_rowx(current_x_prop.cols()); // nnz = 0
            current_x_prop.row(qindex) = empty_rowx;                          
            Eigen::SparseVector<bool> empty_rowy(current_y_prop.cols()); // nnz = 0
            current_y_prop.row(qindex) = empty_rowy;    
            Eigen::SparseVector<bool> empty_rowz(current_z_prop.cols()); // nnz = 0
            current_z_prop.row(qindex) = empty_rowz;    
        }
        /*
        *   When the gate is a CNOT
        */
        if(name=="cnot"){
            size_t qcontrol=gate.qubits[0];           
            size_t qtarget=gate.qubits[1];
        // XOR control ← control ⊕ target
        
        xorRow(current_x_prop, qcontrol,
            current_x_prop, qtarget);

        // XOR target ← target ⊕ control
        xorRow(current_z_prop, qtarget,
                current_z_prop, qcontrol);

        // update Y-propagations similarly
        xorRow(current_y_prop, qcontrol,
                current_x_prop, qtarget);

        xorRow(current_y_prop, qtarget,
                current_z_prop, qcontrol);
        continue;
        }
        if(name=="h"){
            size_t qindex=gate.qubits[0];
            SpVec rowx=current_x_prop.innerVector(qindex);
            SpVec rowz=current_z_prop.innerVector(qindex);
            rowx.swap(rowz);
        }
    }
    // auto t1 = clock::now();                               // stop section‑1
    // auto compile_us = std::chrono::duration_cast<microsec>(t1 - t0).count();
    // std::cout << "[Backword detector matrix construction:] " << compile_us / 1'000.0 << "ms\n";



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



const std::vector<Row>& QEPG::get_parityPropMatrixTrans() const noexcept{
    return parityPropMatrixTranspose_;
}



const SpMat& QEPG::get_detectorMatrix_Eigen() const noexcept{
    return detectorMatrixEigen_;
} 

const SpMat& QEPG::get_dectorMatrixTrans_Eigen() const noexcept{
    return detectorMatrixTransposeEigen_;
}

const SpMat& QEPG::get_parityPropMatrix_Eigen() const noexcept{
    return parityPropMatrixEigen_;
}

const SpMat& QEPG::get_parityPropMatrixTrans_Eigen() const noexcept{
    return parityPropMatrixTransposeEigen_;
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

    using clock     = std::chrono::steady_clock;          // monotonic, good for benchmarking
    using microsec  = std::chrono::microseconds;
    auto t0 = clock::now();                               // start timer

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
    // auto t1 = clock::now();                               // stop section‑1
    // auto compile_us = std::chrono::duration_cast<microsec>(t1 - t0).count();
    // std::cout << "[Set up parity group:] " << compile_us / 1'000.0 << "ms\n";



    // std::cout << "[dectector matrix size:] " << detectorMatrix_.size()<<","<<detectorMatrix_[0].size()<< "\n";
    // std::cout << "[dectector matrix transpose size:] " <<detectorMatrixTranspose_.size()<<","<<detectorMatrixTranspose_[0].size()<< "\n";
    // t0 = clock::now();                               // start timer
    parityPropMatrix_=bitset_matrix_multiplication(paritygroupMatrix,detectorMatrix_);
    // t1 = clock::now();                               // stop section‑1
    // compile_us = std::chrono::duration_cast<microsec>(t1 - t0).count();
    // std::cout << "[bitset_matrix_multiplication(paritygroupMatrix,detectorMatrix_):] " << compile_us / 1'000.0 << "ms\n";
    // t0 = clock::now();                               // start timer
    transpose_matrix(parityPropMatrix_,parityPropMatrixTranspose_);
    // t1 = clock::now();                               // stop section‑1
    // compile_us = std::chrono::duration_cast<microsec>(t1 - t0).count();
    // std::cout << "[transpose_matrix(parityPropMatrix_,parityPropMatrixTranspose_):] " << compile_us / 1'000.0 << "ms\n";
    
}



void QEPG::compute_parityPropMatrix_Eigen(){
    // using clock     = std::chrono::steady_clock;          // monotonic, good for benchmarking
    // using microsec  = std::chrono::microseconds;
    // auto t0 = clock::now();                               // start timer

    // const std::vector<clifford::paritygroup>& detector_parity_group=circuit_.get_detector_parity_group();
    // const clifford::paritygroup& observable_group=circuit_.get_observable_parity_group();
    // const size_t row_size=detector_parity_group.size()+1;
    // const size_t col_size=circuit_.get_num_meas();

    // Eigen::Matrix<uint8_t,Eigen::Dynamic,Eigen::Dynamic> paritygroupMatrix(row_size,col_size);

    // for(size_t i=0; i<detector_parity_group.size();i++){
    //     for(size_t index: detector_parity_group[i].indexlist){
    //         paritygroupMatrix(i,index)=1;
    //     }
    // }
    // for(size_t index: observable_group.indexlist){
    //     paritygroupMatrix(detector_parity_group.size(),index)=1;
    // }
    // auto t1 = clock::now();                               // stop section‑1
    // auto compile_us = std::chrono::duration_cast<microsec>(t1 - t0).count();
    // std::cout << "[Set up parity group:] " << compile_us / 1'000.0 << "ms\n";



    // Eigen::Matrix<uint8_t,Eigen::Dynamic,Eigen::Dynamic>  tmpdetectorMatrixEigen  =
    //         detectorMatrixEigen_.cast<uint8_t>();                 // 0/1 → 0/1 as bytes


    // t0 = clock::now();                               // start timer
    // parityPropMatrixEigen_ =
    // ( (paritygroupMatrix * tmpdetectorMatrixEigen)      // uint8 product
    //   .array()                                          //  ⟶ Array<uint8>
    //   .template cast<bool>() )                          //  ⟶ Array<bool>
    // .matrix();                                          //  ⟶ Matrix<bool>
    // t1 = clock::now();                               // stop section‑1
    // compile_us = std::chrono::duration_cast<microsec>(t1 - t0).count();
    // std::cout << "[Eigen3_matrix_multiplication(paritygroupMatrix,detectorMatrix_):] " << compile_us / 1'000.0 << "ms\n";
    // t0 = clock::now();                               // start timer
    // parityPropMatrixTransposeEigen_=parityPropMatrixEigen_.transpose();
    // t1 = clock::now();                               // stop section‑1
    // compile_us = std::chrono::duration_cast<microsec>(t1 - t0).count();
    // std::cout << "[transpose_matrix(parityPropMatrix_,parityPropMatrixTranspose_):] " << compile_us / 1'000.0 << "ms\n";
}



}