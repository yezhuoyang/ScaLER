"""
S-curve model implementations for ScaLER.

This package provides multiple S-curve models for fitting logical error rate data:
- OurScurveModel: The model from Definition 1 in the paper
- IBMScurveModel: IBM's min-fail enclosure model from Definition 2

Usage:
    from scalerqec.Stratified.models import ModelType, ModelFactory, OurScurveModel

    # Create a specific model
    model = ModelFactory.create(ModelType.OUR_MODEL, t=3, gamma=0.05)

    # Or create all models for comparison
    models = ModelFactory.create_all(t=3, gamma=0.05)
"""

from scalerqec.Stratified.models.base import ScurveModelBase
from scalerqec.Stratified.models.our_model import OurScurveModel
from scalerqec.Stratified.models.ibm_model import IBMScurveModel
from scalerqec.Stratified.models.factory import ModelType, ModelFactory

__all__ = [
    "ScurveModelBase",
    "OurScurveModel",
    "IBMScurveModel",
    "ModelType",
    "ModelFactory",
]
