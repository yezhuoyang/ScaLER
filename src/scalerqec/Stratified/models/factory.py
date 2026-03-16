"""
Factory for creating S-curve model instances.

This module provides:
- ModelType enum for identifying different model types
- ModelFactory class for creating model instances
"""

from enum import Enum
from typing import Dict, Type

from scalerqec.Stratified.models.base import ScurveModelBase
from scalerqec.Stratified.models.our_model import OurScurveModel
from scalerqec.Stratified.models.ibm_model import IBMScurveModel


class ModelType(Enum):
    """Enum for S-curve model types."""

    OUR_MODEL = "our_model"
    IBM_MODEL = "ibm_model"


class ModelFactory:
    """
    Factory for creating S-curve model instances.

    Usage:
        # Create a specific model
        model = ModelFactory.create(ModelType.OUR_MODEL, t=3, gamma=0.05)

        # Create all models for comparison
        models = ModelFactory.create_all(t=3, gamma=0.05)

        # Get list of available models
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
