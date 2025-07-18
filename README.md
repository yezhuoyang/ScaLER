# ScaLER

Scalable testing of logical error rate for QEC circuit. 



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

```bash
py setup.py clean --all    
```


We also need to convert C++ object to python object directly. So "Python.h" needs to be added to the search path. Typically, it is under:


```bash
C:\Users\yezhu\miniconda3\include
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

# How to run tests?


All test script are kept under test/ folder. You can test the correct ness of our QEPG implementation, test with ground truth for small scale circuit. 





# How to run benchmark circuit and reproduce the table?

All our benchmark circuit are stored undered stimprograms/ folder. To reproduce one circuit, for example, for surface code with distance 7, execute the following python script:



```python
from ScaLER import stratified_Scurve_LERcalc
d=7
p = 0.001
repeat=5
sample_budget = 100_000_0000
t = (d - 1) // 2
stim_path = f"your/path/stimprograms/surface/surface7"
figname = f"Surface{d}"
titlename = f"Surface{d}"
output_filename = f"Surface{d}.txt"
testinstance = stratified_Scurve_LERcalc(p, sampleBudget=sample_budget, k_range=5, num_subspace=6, beta=4)
testinstance.set_t(t)
testinstance.set_sample_bound(
    MIN_NUM_LE_EVENT=100,
    SAMPLE_GAP=100,
    MAX_SAMPLE_GAP=5000,
    MAX_SUBSPACE_SAMPLE=50000
)
with open(output_filename, "w") as f:
    with redirect_stdout(f):
        testinstance.calculate_LER_from_file(stim_path, p, 0, figname, titlename, repeat)
```



# Symbolic calculator of LER

In this part, I explain how to get the ground truth of logical error rate by Symbolic calculator. Reader can reproduce Table 2



```python
from ScaLER import symbolicLER

tmp=symbolicLER(0.001)
filepath="your/file/path/to/circuit"
print(tmp.calculate_LER_from_file(filepath,0.001))
p=0.001

num_noise=tmp._num_noise

# for weight in range(1,11):
#     print("LER in the subspace {} is {}".format(weight,tmp.evaluate_LER_subspace(p,weight)))        


for weight in range(1,12):
    print("SubspaceLER {} is {}".format(weight,tmp.subspace_LER(weight)))     
```




# Use Monte to test any circuit

In this part, I explain how to test any circuit with the widely use random fault injection method.



```python
from ScaLER import stimLERcalc


p=0.001
filepath="your/file/path/to/circuit"
dlist=[15]
repeat=5
stim_path = base_dir+rel+str(d)
# 3) build your output filename:
out_fname =  result_dir+str(p)+"-"+str(code_type)+str(d)+"-resultMonte.txt"     # e.g. "surface3-result.txt"
# 4) redirect prints for just this file:

with open(out_fname, "w") as outf, redirect_stdout(outf):
    print(f"---- Processing {stim_path} ----")

    calculator=stimLERcalc(10)
    # pass the string path into your function:
    ler = calculator.calculate_LER_from_my_random_sampler(500000000, filepath, p,repeat)
```




# Use ScaLER to test any circuit

In this part, I explain how to use ScaLER to test and input circuit. I will explain how to change hyper parameters. In a python script, run the folloing code:


```python
from ScaLER import stratified_Scurve_LERcalc
p = 0.001
repeat=5
sample_budget = 100_000_0000
t = (d - 1) // 2
stim_path = f"C:/Users/yezhu/Documents/Sampling/stimprograms/surface/surface{d}"
figname = f"Surface{d}"
titlename = f"Surface{d}"
output_filename = f"Surface{d}.txt"
testinstance = stratified_Scurve_LERcalc(p, sampleBudget=sample_budget, k_range=5, num_subspace=6, beta=4)
testinstance.set_t(t)
testinstance.set_sample_bound(
    MIN_NUM_LE_EVENT=100,
    SAMPLE_GAP=100,
    MAX_SAMPLE_GAP=5000,
    MAX_SUBSPACE_SAMPLE=50000
)
with open(output_filename, "w") as f:
    with redirect_stdout(f):
        testinstance.calculate_LER_from_file(stim_path, p, 0, figname, titlename, repeat)
```

Under the same directory, you will see a output figure which shows the fitted curve, and a output.txt file which shows all results. 







