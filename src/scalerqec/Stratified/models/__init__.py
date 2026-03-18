"""S-curve model implementations for the ScaLER algorithm.

This subpackage provides pluggable S-curve models that describe how the
conditional logical error probability P_L(w) varies with fault weight *w*.
Two concrete models are included:

- :class:`OurScurveModel` (Definition 1): a logistic S-curve with a
  pole term at the fault-tolerant threshold,
  ``P_L(w) = 0.5 / (1 + exp(-(w-mu)/alpha + beta/sqrt(w-t)))``.
- :class:`IBMScurveModel` (Definition 2): IBM's "min-fail enclosure"
  model with power-law onset,
  ``P_L(w) = 0.5 * [1 - exp(-2*f0 * (w/w0)^gamma_ibm)]``.

All models derive from :class:`ScurveModelBase` and can be created via
:class:`ModelFactory`::

    from scalerqec.Stratified.models import ModelType, ModelFactory

    model = ModelFactory.create(ModelType.OUR_MODEL, t=3, gamma=0.05)
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
