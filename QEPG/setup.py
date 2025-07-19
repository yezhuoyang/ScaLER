import setuptools
import pybind11
from pybind11.setup_helpers import Pybind11Extension
import sys



# pick the right “O3”‑style flag
extra_compile_args = ['/std:c++20', '/EHsc']
if sys.platform == 'win32':
    extra_compile_args.append('/O2')    # MSVC “optimize for speed”
    extra_compile_args.append('/openmp:llvm')   
    lflags = [
        "/DEBUG"               # <-- carry PDB into final .pyd
    ]
else:
    extra_compile_args.append('-O3')    # GCC/Clang “optimize even more”


ext_modules=[
    setuptools.Extension(
        'QEPG',
        ['bindings.cpp','src/QEPG.cpp','src/clifford.cpp','src/sampler.cpp','src/LERcalculator.cpp'],
        include_dirs=[
            pybind11.get_include(),
            '.',
            'C:/local/boost_1_87_0/',
            'C:/Users/username/OneDrive/Documents/GitHub/vcpkg/installed/x64-windows/include',
            'C:/vcpkg/installed/x64-windows/include',
            'C:/Users/username/AppData/Local/Programs/Python/Python311/Include',
            'C:/Users/username/miniconda3/Include',
        ],
        extra_compile_args=extra_compile_args,
        extra_link_args=lflags,
        language='c++'
    ),
]


# Configure setuptools setup
setuptools.setup(
    name='QEPG', # A package name
    version='0.0.1', # Version number
    author="anonymous",
    author_email="anonymous",
    description="QEPG",
    ext_modules=ext_modules,
    zip_safe=False, # Important: Set to False for C++ extensions
)