$Env:VCPKG_ROOT = 'C:\vcpkg'
g++ -std=c++20 -I"C:\local\boost_1_87_0" -I"$Env:VCPKG_ROOT\installed\x64-windows\include" -Wall -Wextra -pedantic -O3 main.cpp  QEPG/src/LERcalculator.cpp QEPG/src/clifford.cpp QEPG/src/QEPG.cpp QEPG/src/sampler.cpp -o demo
./demo