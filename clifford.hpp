#ifndef CLIFFORD_HPP
#define CLIFFORD_HPP
#pragma once

#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

namespace clifford{


class cliffordcircuit{

public:

    cliffordcircuit();
    explicit cliffordcircuit(int n_qubit);
    void set_error_rate(double p);


    /* single-qubit noise----------------------------------*/
    void add_XError(int qindex);
    void add_ZError(int qindex);
    void add_depolarize1(int qindex);

    /*Clifford gates---------------------------------------*/
    void add_cnot(int qcontrol, int qtarget);
    void add_hadamard(int qindex);
    void add_phase(int qindex);
    void add_pauliX(int qindex);
    void add_pauliy(int qindex);
    void add_pauliz(int qindex);

    /*Reset/Measurement---------------------------------------*/
    void add_reset(int qindex);
    void add_measurement(int qindex);

    /*misc -------------------------------------------------*/
    void print_circuit() const;


private:


    struct Gate{
        std::string name;
        std::vector<int> qubits;
    };


    int     num_qubit_;
    double  error_rate_{0.0};


    std::vector<Gate> circuit_;




};

}





#endif