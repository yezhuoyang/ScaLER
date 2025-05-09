# Sampling
This is the compiler implemented in C++ that sampling and calculate logical error rate. 



# Environment and compilation


> [!IMPORTANT]
> Be careful about setting up environment

The code is developed in C++17. 


We are using dynamic bitset from boost, to install it, run:

```bash
choco install   boost-msvc-14.3 -y
```

The boost header file will be stored under the path "C:\local\boost_1_87_0\boost". Add this path into VScode cpp include path. 

<<<<<<< HEAD
The matrix size turn out to be the bottleneck. 
=======


> [!TIP]
> Avoid using dense matrix!


>>>>>>> bb98f20e0fdc7bfccc5a908c8dd02c3bbbe337fe
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


- [ ] **Integrate pybind11**: Enable Python to directly call C++ functions  
- [ ] **Verify QEPG correctness on small-scale examples**  
- [ ] **Construct the parity propagation matrix directly** during detector matrix construction  
- [ ] **Validate QEPG at larger scale** by comparing simulation results with STIM  
- [ ] **Optimize performance** using multithreading or other parallelism strategies  
