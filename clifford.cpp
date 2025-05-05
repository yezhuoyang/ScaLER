#include "clifford.hpp"
#include <iomanip>
#include <iostream>
#include <string_view>
#include <charconv>   // std::from_chars

namespace clifford{



/*-------------------------------------------ctor */

cliffordcircuit::cliffordcircuit() = default;

cliffordcircuit::cliffordcircuit(size_t n_qubit)
        : num_qubit_(n_qubit){}

/*-------------------------------------------configuration*/


void cliffordcircuit::set_error_rate(double p){error_rate_=p;}


/*--------------------------------------------Add quantum noise*/


void cliffordcircuit::add_XError(size_t qindex) {
    circuit_.push_back({"X_ERROR", {qindex}});
    num_qubit_=std::max(num_qubit_,qindex+1);
}

void cliffordcircuit::add_ZError(size_t qindex) {
    circuit_.push_back({"X_ZRROR", {qindex}});
    num_qubit_=std::max(num_qubit_,qindex+1);
}

void cliffordcircuit::add_depolarize1(size_t qindex) {
    circuit_.push_back({"DEPOLARIZE1", {qindex}});
    num_qubit_=std::max(num_qubit_,qindex+1);
}



/*--------------------------------------------1 qubit gate*/


void cliffordcircuit::add_hadamard(size_t qindex){
    circuit_.push_back({"h", {qindex}});
    num_qubit_=std::max(num_qubit_,qindex+1);
}

void cliffordcircuit::add_phase(size_t qindex){
    circuit_.push_back({"p", {qindex}});
    num_qubit_=std::max(num_qubit_,qindex+1);
}

void cliffordcircuit::add_pauliX(size_t qindex){
    circuit_.push_back({"x", {qindex}});
    num_qubit_=std::max(num_qubit_,qindex+1);
}

void cliffordcircuit::add_pauliy(size_t qindex){
    circuit_.push_back({"y", {qindex}});
    num_qubit_=std::max(num_qubit_,qindex+1);
}

void cliffordcircuit::add_pauliz(size_t qindex){
    circuit_.push_back({"z", {qindex}});
    num_qubit_=std::max(num_qubit_,qindex+1);
}


/*--------------------------------------------2 qubit gate---------------------------*/


void cliffordcircuit::add_cnot(size_t qcontrol, size_t qtarget){
    circuit_.push_back({"cnot", {qcontrol,qtarget}});
    num_qubit_=std::max(num_qubit_,qcontrol+1);
    num_qubit_=std::max(num_qubit_,qtarget+1);    
}


/*--------------------------------------------Reset/Measurement gadget---------------*/

void cliffordcircuit::add_reset(size_t qindex){
    circuit_.push_back({"R", {qindex}});
    num_qubit_=std::max(num_qubit_,qindex+1);
}

void cliffordcircuit::add_measurement(size_t qindex){
    circuit_.push_back({"M", {qindex}});
    num_qubit_=std::max(num_qubit_,qindex+1);
}


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
const Gate& cliffordcircuit::get_gate(size_t gateindex) const{return circuit_.at(gateindex);}

/*Get member-------------------------------------------*/
int cliffordcircuit::get_num_qubit() const{return num_qubit_;}
int cliffordcircuit::get_gate_num() const{return circuit_.size();}

void cliffordcircuit::set_num_qubit(size_t num_qubit) {num_qubit_=num_qubit;}

int cliffordcircuit::get_num_meas() const{
    return num_meas_;
}

int cliffordcircuit::get_num_noise() const{
    return num_noise_;
}


/*Helper functions for parsing------------------------------------------------------------*/


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


inline std::string_view next_token(std::string_view& sv){
    //First, we discard all leading spaces/tabs
    const auto first=sv.find_first_not_of(" \t");
    if(first==std::string_view::npos){
        sv={};
        return {};
    }
    sv.remove_prefix(first);

    //Token is substring [0, pos)
    const auto pos=sv.find_first_of(" \t");
    std::string_view token=sv.substr(0,pos);
    sv.remove_prefix(token.size());
    return token;
}


inline int to_int(std::string_view tok){
    int value{};
    std::from_chars(tok.data(),tok.data()+tok.size(),value);
    return value;
}


inline std::size_t to_size_t(std::string_view tok){
    std::size_t value{};
    const char* begin = tok.data();
    const char* end= begin + tok.size();

    auto [ptr, ec] = std::from_chars(begin,end,value,10);

    if (ec == std::errc::invalid_argument || ptr != end)
        throw std::invalid_argument{"token is not a non‑negative integer"};
    if (ec == std::errc::result_out_of_range)
        throw std::out_of_range{"integer value exceeds std::size_t range"};

    return value;
}






/*compile from stim string---------------------------------------------------*/
void cliffordcircuit::compile_from_rewrited_stim_string(std::string stim_str){
    
    for_each_line(stim_str, [this](std::string_view line){
        std::string_view rest=line;
        std::string_view op=next_token(rest);

        if(op=="M"){
            size_t qindex=to_size_t(next_token(rest));
            add_measurement(qindex);
        }
        else if(op=="R"){
            size_t qindex=to_size_t(next_token(rest));
            add_reset(qindex);           
        }
        else if(op=="H"){
            size_t qindex=to_size_t(next_token(rest));
            add_hadamard(qindex);   
        }       
        else if(op=="CX"){
            size_t qcontrol=to_size_t(next_token(rest)); 
            size_t qtarget=to_size_t(next_token(rest));
            add_cnot(qcontrol,qtarget);
        }
        else if(op=="DETECTOR"){
           
        }
        else if(op=="OBSERVABLE"){

        }
    });


}



}