#include "clifford.hpp"
#include <iomanip>
#include <iostream>
#include <string_view>

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

void cliffordcircuit::set_num_qubit(int num_qubit) {num_qubit_=num_qubit;}




/*Helper function------------------------------------------------------------*/


template <typename Callback>
void for_each_line(std::string_view sv, Callback&& callback){
    while(!sv.empty()){
        auto pos = sv.find_first_of("\n\r");
        std::string_view line=sv.substr(0,pos);
        callback(line);

        if(pos == sv.npos) break;
        sv.remove_prefix(pos+1);

        if(!sv.empty()&&sv.front()=='\n'&& line.back()=='\r')
        sv.remove_prefix(1);
    }
}






/*compile from stim string---------------------------------------------------*/
void cliffordcircuit::compile_from_rewrited_stim_string(std::string stim_str){
    
    for_each_line(stim_str, [](std::string_view line){
        std::cout<<line<<"\n";
    });


}



}