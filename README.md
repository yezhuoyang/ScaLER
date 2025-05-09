# Sampling
This is the compiler implemented in C++ that sampling and calculate logical error rate. 



# Environment and compilation


The code is developed in C++17. 


We are using dynamic bitset from boost, to install it, run:

```bash
choco install   boost-msvc-14.3 -y
```

The boost header file will be stored under the path "C:\local\boost_1_87_0\boost". Add this path into VScode cpp include path. 

The matrix size turn out to be the bottleneck. 
We also use vcpkg and install the Eigen3 library for matrix operations:


https://chatgpt.com/share/681cbfa0-f8d8-8005-9878-a8798ff9a88a



We use pybind11 to convert the samples from C++ objects to python objects. To install using vcpkg, run the following command:

```bash
vcpkg install pybind11
```

Also need to add the path of the python header file

```bash
py -c "from sysconfig import get_paths as gp; print(gp()['include'])"
```


# Plan of development


[ ] Add pybinding, so python can directly call the C++ function
[ ] Check the correctness of QEPG at small scale
[ ] In the construction of detectormatrix, ditectly construct the parityprop matrix
[ ] Check the correctness of QEPG at larger scale: Compare the simulation result with STIM
[ ] Optimize the speed by parallism