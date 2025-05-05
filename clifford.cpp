#include "clifford.hpp"
#include <iomanip>
#include <iostream>


namespace clifford{



/*-------------------------------------------ctor */

cliffordcircuit::cliffordcircuit() = default;

cliffordcircuit::cliffordcircuit(int n_qubit)
        : num_qubit_(n_qubit){}

/*-------------------------------------------configuration*/


void cliffordcircuit::set_error_rate(double p){error_rate_=p;}


/*--------------------------------------------Add quantum noise*/


void cliffordcircuit::add_XError(int qindex) {circuit_.push_back({"X_ERROR", {qindex}});}
void cliffordcircuit::add_ZError(int qindex) {circuit_.push_back({"X_ZRROR", {qindex}});}
void cliffordcircuit::add_depolarize1(int qindex) {circuit_.push_back({"DEPOLARIZE1", {qindex}});}



/*--------------------------------------------1 qubit gate*/


void cliffordcircuit::add_hadamard(int qindex){circuit_.push_back({"h", {qindex}});}
void cliffordcircuit::add_phase(int qindex){circuit_.push_back({"p", {qindex}});}
void cliffordcircuit::add_pauliX(int qindex){circuit_.push_back({"x", {qindex}});}
void cliffordcircuit::add_pauliy(int qindex){circuit_.push_back({"y", {qindex}});}
void cliffordcircuit::add_pauliz(int qindex){circuit_.push_back({"z", {qindex}});}


/*--------------------------------------------2 qubit gate---------------------------*/


void cliffordcircuit::add_cnot(int qcontrol, int qtarget){circuit_.push_back({"cnot", {qcontrol,qtarget}});}


/*--------------------------------------------Reset/Measurement gadget---------------*/

void cliffordcircuit::add_reset(int qindex){circuit_.push_back({"R", {qindex}});}
void cliffordcircuit::add_measurement(int qindex){circuit_.push_back({"M", {qindex}});}


/*---------------------------------------------visualizatation-----------------------*/


void cliffordcircuit::print_circuit() const{
    std::cout<<"\n--- Clifford circuit ("<<num_qubit_<<" qubits, p= "<< std::setprecision(4)<< error_rate_<<")-------\n";

    for(std::size_t i=0; i< circuit_.size(); ++i){
        std::cout<<std::setw(4)<<i << ":"<<std::left<<std::setw(10)<<circuit_[i].name <<" ";
        for(auto q: circuit_[i].qubits) std::cout<<"q"<<q<<" ";
        std::cout << '\n';
    }
    std::cout<<"-------------------------------------------------------------\n";
}



/*--Get gate by index-------------------------------------*/
const Gate& cliffordcircuit::get_gate(int gateindex) const{return circuit_.at(gateindex);}

/*Get member-------------------------------------------*/
int cliffordcircuit::get_num_qubit() const{return num_qubit_;}
int cliffordcircuit::get_gate_num() const{return circuit_.size();}



}