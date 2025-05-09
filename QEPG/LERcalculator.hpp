#ifndef LER_HPP
#define LER_HPP
#pragma once

#include "QEPG.hpp"
#include "sampler.hpp"



namespace LERcalculator{

 std::vector<std::vector<bool>> return_samples(std::string prog_str,size_t weight, size_t shots);


}


#endif