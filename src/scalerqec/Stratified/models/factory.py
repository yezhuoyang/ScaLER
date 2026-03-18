from __future__ import annotations

"""Factory pattern for S-curve model instantiation.

This module provides:

- :class:`ModelType` -- an enum identifying the available model types.
- :class:`ModelFactory` -- a class-based factory for creating
  :class:`~scalerqec.Stratified.models.base.ScurveModelBase` instances
  by type.  Custom models can be registered at runtime via
  :meth:`ModelFactory.register`.
"""

from enum import Enum
from typing import Dict, Type

from scalerqec.Stratified.models.base import ScurveModelBase
from scalerqec.Stratified.models.our_model import OurScurveModel
from scalerqec.Stratified.models.ibm_model import IBMScurveModel


class ModelType(Enum):
    """Enumeration of available S-curve model types.

    Attributes:
        OUR_MODEL: Definition 1 model with pole term
            (:class:`~scalerqec.Stratified.models.our_model.OurScurveModel`).
        IBM_MODEL: Definition 2 IBM min-fail enclosure model
            (:class:`~scalerqec.Stratified.models.ibm_model.IBMScurveModel`).
    """

    OUR_MODEL = "our_model"
    IBM_MODEL = "ibm_model"


class ModelFactory:
    """Factory for creating :class:`ScurveModelBase` instances by type.

    The factory maintains an internal registry mapping
    :class:`ModelType` values to concrete model classes.  New model
    types can be added at runtime with :meth:`register`.

    Example::

        model = ModelFactory.create(ModelType.OUR_MODEL, t=3, gamma=0.05)
        all_models = ModelFactory.create_all(t=3, gamma=0.05)
        available = ModelFactory.available_models()
    """

    _registry: Dict[ModelType, Type[ScurveModelBase]] = {
        ModelType.OUR_MODEL: OurScurveModel,
        ModelType.IBM_MODEL: IBMScurveModel,
    }

    @classmethod
    def register(
        cls, model_type: ModelType, model_class: Type[ScurveModelBase]
    ) -> None:
        """
        Register a new model type.

        Args:
            model_type: The enum value for this model
            model_class: The model class to instantiate
        """
        cls._registry[model_type] = model_class

    @classmethod
    def create(cls, model_type: ModelType, **kwargs) -> ScurveModelBase:
        """
        Create an instance of the specified model type.

        Args:
            model_type: The type of model to create
            **kwargs: Arguments passed to the model constructor (e.g., t, gamma)

        Returns:
            An instance of the specified model

        Raises:
            ValueError: If model_type is not registered
        """
        if model_type not in cls._registry:
            raise ValueError(f"Unknown model type: {model_type}")
        return cls._registry[model_type](**kwargs)

    @classmethod
    def create_all(cls, **kwargs) -> Dict[ModelType, ScurveModelBase]:
        """
        Create instances of all registered models.

        Args:
            **kwargs: Arguments passed to each model constructor

        Returns:
            Dictionary mapping ModelType to model instance
        """
        return {mt: cls.create(mt, **kwargs) for mt in cls._registry}

    @classmethod
    def available_models(cls) -> list[ModelType]:
        """
        Get list of available model types.

        Returns:
            List of registered ModelType values
        """
        return list(cls._registry.keys())

    @classmethod
    def get_model_class(cls, model_type: ModelType) -> Type[ScurveModelBase]:
        """
        Get the class for a model type without instantiating.

        Args:
            model_type: The type of model

        Returns:
            The model class

        Raises:
            ValueError: If model_type is not registered
        """
        if model_type not in cls._registry:
            raise ValueError(f"Unknown model type: {model_type}")
        return cls._registry[model_type]
