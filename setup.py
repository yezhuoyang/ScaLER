import os
import sys
import setuptools
import pybind11
from pybind11.setup_helpers import Pybind11Extension, build_ext

# Base include dirs: pybind11 + your source tree
include_dirs = [
    pybind11.get_include(),
    "QEPG",
    "QEPG/src",
]

extra_compile_args = []
extra_link_args = []


if sys.platform == "win32":
    # --- Windows flags ---
    extra_compile_args = ["/std:c++20", "/EHsc", "/O2", "/openmp:llvm"]
    extra_link_args = ["/DEBUG"]

elif sys.platform == "darwin":
    # --- macOS flags (Apple Clang does not support -fopenmp directly) ---
    homebrew_prefix = os.environ.get("HOMEBREW_PREFIX", "/opt/homebrew")
    omp_inc = os.path.join(homebrew_prefix, "opt", "libomp", "include")
    omp_lib = os.path.join(homebrew_prefix, "opt", "libomp", "lib")
    extra_compile_args = ["-std=c++20", "-O3", "-Xpreprocessor", "-fopenmp"]
    extra_link_args = ["-lomp"]
    if os.path.isdir(omp_inc):
        extra_compile_args.append(f"-I{omp_inc}")
    if os.path.isdir(omp_lib):
        extra_link_args.append(f"-L{omp_lib}")

else:
    # --- Linux flags ---
    extra_compile_args = ["-std=c++20", "-O3", "-fopenmp"]
    extra_link_args = ["-fopenmp"]


ext_modules = [
    Pybind11Extension(
        "scalerqec.qepg",
        [
            "QEPG/bindings.cpp",
            "QEPG/src/QEPG.cpp",
            "QEPG/src/clifford.cpp",
            "QEPG/src/sampler.cpp",
            "QEPG/src/LERcalculator.cpp",
        ],
        include_dirs=include_dirs,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
        language="c++",
    )
]


setuptools.setup(
    name="scalerqec",
    version="1.0.0",
    description="Scalable Quantum Error Correction testing Tools for logical error rate and software correctness",
    author="John Ye",
    packages=setuptools.find_packages(where="src"),
    package_dir={"": "src"},
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
)
