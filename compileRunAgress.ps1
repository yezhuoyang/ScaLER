$Env:VCPKG_ROOT = 'C:\vcpkg'
g++ -std=c++20 -I"C:\local\boost_1_87_0" -I"$Env:VCPKG_ROOT\installed\x64-windows\include"  -O3 main.cpp QEPG/clifford.cpp QEPG/QEPG.cpp QEPG/sampler.cpp -o demo
./demo