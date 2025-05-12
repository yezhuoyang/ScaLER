#ifndef LER_HPP
#define LER_HPP
#pragma once

#include "QEPG.hpp"
#include "sampler.hpp"
#include <chrono>


namespace LERcalculator{

 std::vector<std::vector<bool>> return_samples(const std::string & prog_str,size_t weight, size_t shots);


 std::vector<std::vector<std::vector<bool>>> return_samples_many_weights(const std::string& prog_str,const std::vector<size_t>& weight, const std::vector<size_t>& shots);



}


#endif