"""
ScaLERQEC -- Scalable Logical Error Rate Estimation for Quantum Error Correction
================================================================================

ScaLERQEC is a Python toolkit for estimating the logical error rate (LER) of
stabilizer codes under realistic noise models.  It provides multiple estimation
strategies so users can trade off accuracy, speed, and code distance:

* **Stratified sampling** with S-curve fitting (the ScaLER algorithm)
* **Exact symbolic** computation via dynamic programming
* **Monte Carlo** estimation via stim detector sampling

Key classes
-----------
:class:`StratifiedLERcalc`
    Stratified fault-injection LER estimator.
:class:`StratifiedScurveLERcalc`
    Stratified estimator with S-curve fitting for improved accuracy.
:class:`SymbolicLERcalc`
    Exact symbolic LER computed as a polynomial in the physical error rate.
:class:`MonteLERcalc`
    Standard Monte Carlo LER estimator using stim detector sampling.
:class:`CliffordCircuit`
    Clifford circuit representation with STIM import and QEPG construction.
:class:`StabCode`
    Stabilizer code definition bundling a code with its noise model.

Quick example
-------------
::

    from scalerqec import StabCode, StratifiedScurveLERcalc

    code = StabCode.from_stim("path/to/surface_code.stim")
    calc = StratifiedScurveLERcalc(code)
    result = calc.run()
    print(f"Estimated LER: {result}")
"""

# Expose C++ backend as scaler.qepg
from . import qepg

# Re-export high-level components for easy access
from .Stratified.stratifiedLER import StratifiedLERcalc
from .Stratified.stratifiedScurveLER import StratifiedScurveLERcalc
from .Symbolic.symbolicLER import SymbolicLERcalc
from .Monte.monteLER import MonteLERcalc
from .Clifford.clifford import CliffordCircuit
from .QEC.qeccircuit import StabCode
from .Clifford.noiselabel import NoiseLabel, NoiseLabelMap
from .QEC.autolabel import auto_label_from_circuit, auto_label_from_stabcode
from .Analysis.hotspot import HotspotAnalyzer, HotspotResult, CategoryReport
from .Analysis.visualize import labeled_circuit_latex


__all__ = [
    "StratifiedLERcalc",
    "StratifiedScurveLERcalc",
    "SymbolicLERcalc",
    "MonteLERcalc",
    "CliffordCircuit",
    "qepg",
    "StabCode",
    "NoiseLabel",
    "NoiseLabelMap",
    "auto_label_from_circuit",
    "auto_label_from_stabcode",
    "HotspotAnalyzer",
    "HotspotResult",
    "CategoryReport",
    "labeled_circuit_latex",
]
