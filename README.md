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



> [!TIP]
> Avoid using dense matrix!


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

Run the following command to build the QEPG package with pybinding:

```bash
py setup.py build_ext --inplace
```

Run the following command to clear the previously compiled output:

```bind
py setup.py clean --all    
```


# How to compile run python script


To compile QEPG python package by pybind11:

```bash
(Under QEPG folder)./compilepybind.ps1
```

The python code is divided into different modules. For example, to run the test_by_stim.py file under test module, stay at the root folder and execute:

```bash
(Under Sampling folder)py -m test.test_by_stim   
```

# How to run tests


# How to run benchmark




# Plan of development


- [X] **Integrate pybind11**: Enable Python to directly call C++ functions  
- [ ] **Optimize python interface**: Optimize the interface which convert C++ generated samples to python samples
- [X] **Verify QEPG correctness on small-scale examples**  
- [X] **Construct the parity propagation matrix directly** during detector matrix construction  
- [X] **Validate QEPG at larger scale** by comparing simulation results with STIM  
- [ ] **Optimize performance** using multithreading or other parallelism strategies  
- [X] **LER symbolic calculator** Calculate logical error rate(of small circuit) by dynamic algorithm
- [ ] **LER calculation benchmark(Exact)** Set up a benchmark of small circuit with known logical error rate
- [ ] **LER calculation benchmark(Compare with STIM)** Set up a benchmark of surface code/repetition code
- [ ] **QEPG generating speed benchmark** Set up a benchmark comparing the time used to generate the QEPG graph
- [ ] **Sampling rate benchmark** Set up a benchmark comparing the time used to get 1 million samples
- [ ] **Truncate weight of symbolic calculator**  Add a maximum weight in the dynamic algorithm, since we only care about small weight(W=1,2,3)