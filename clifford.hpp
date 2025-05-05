#ifndef CLIFFORD_HPP
#define CLIFFORD_HPP
#pragma once

#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

namespace clifford{

struct Gate{
    std::string name;
    std::vector<size_t> qubits;
};


struct paritygroup{
    std::vector<size_t> indexlist;
};


class cliffordcircuit{

public:

    cliffordcircuit();
    explicit cliffordcircuit(size_t n_qubit);
    void set_error_rate(double p);


    /* single-qubit noise----------------------------------*/
    void add_XError(size_t qindex);
    void add_ZError(size_t qindex);
    void add_depolarize1(size_t qindex);

    /*Clifford gates---------------------------------------*/
    void add_cnot(size_t qcontrol, size_t qtarget);
    void add_hadamard(size_t qindex);
    void add_phase(size_t qindex);
    void add_pauliX(size_t qindex);
    void add_pauliy(size_t qindex);
    void add_pauliz(size_t qindex);

    /*Reset/Measurement---------------------------------------*/
    void add_reset(size_t qindex);
    void add_measurement(size_t qindex);

    /*misc -------------------------------------------------*/
    void print_circuit() const;


    /*--Get gate by index-------------------------------------*/
    const Gate& get_gate(size_t gateindex) const;

    /*Setter/Getter of class members-------------------------------------------*/
    int get_num_qubit() const;
    void set_num_qubit(size_t num_qubit);
    int get_gate_num() const;
    int get_num_meas() const;
    int get_num_noise() const;

    /*compile from stim string---------------------------------------------------*/
    void compile_from_rewrited_stim_string(std::string stim_str);

private:



    size_t     num_qubit_=0;
    size_t     num_noise_=0;
    size_t     num_meas_=0;

    double  error_rate_{0.0};

    std::vector<Gate> circuit_;
    std::vector<size_t> measureindexList;

    std::vector<paritygroup> detectors_;
    paritygroup observable_;

};

}





#endif