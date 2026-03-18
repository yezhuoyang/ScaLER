"""
Tests for S-curve model implementations.

Tests the model abstraction layer including:
- OurScurveModel (Definition 1 from paper)
- IBMScurveModel (Definition 2 from paper)
- ModelFactory
"""

import pytest
import numpy as np
from typing import List

from scalerqec.Stratified.models import (
    ScurveModelBase,
    OurScurveModel,
    IBMScurveModel,
    ModelType,
    ModelFactory,
)


class TestOurScurveModel:
    """Tests for Our S-curve model (Definition 1)."""

    def test_initialization(self):
        """Test model initialization with default and custom parameters."""
        # Default initialization
        model = OurScurveModel()
        assert model.t == 0
        assert model.gamma == 0.05
        assert model.name == "OurModel"

        # Custom initialization
        model = OurScurveModel(t=5, gamma=0.1)
        assert model.t == 5
        assert model.gamma == 0.1

    def test_param_names(self):
        """Test that param_names returns expected parameters."""
        model = OurScurveModel()
        assert model.param_names == ("alpha", "mu", "beta")

    def test_predict_bounds(self):
        """Test that predictions are in valid range (0, 0.5)."""
        model = OurScurveModel(t=3, gamma=0.05)
        model.set_params(alpha=10.0, mu=50.0, beta=5.0)

        # Test for various weights > t
        for w in range(4, 100):
            p = model.predict(w)
            assert 0 <= p <= 0.5, f"P_L({w}) = {p} is out of bounds"

    def test_predict_fault_tolerant_region(self):
        """Test that predictions are 0 for w <= t."""
        model = OurScurveModel(t=5, gamma=0.05)
        model.set_params(alpha=10.0, mu=50.0, beta=5.0)

        for w in range(0, 6):  # w <= t = 5
            p = model.predict(w)
            assert p == 0.0, f"P_L({w}) = {p} should be 0 in fault-tolerant region"

    def test_predict_array_input(self):
        """Test prediction with array input."""
        model = OurScurveModel(t=2, gamma=0.05)
        model.set_params(alpha=10.0, mu=30.0, beta=5.0)

        weights = np.array([3, 5, 10, 20, 50])
        predictions = model.predict(weights)

        assert isinstance(predictions, np.ndarray)
        assert len(predictions) == len(weights)
        assert all(0 <= p <= 0.5 for p in predictions)

    def test_transform_inverse_roundtrip(self):
        """Test that transform(inverse_transform(y)) ≈ y."""
        model = OurScurveModel()

        y_values = [0.5, 1.0, 2.0, 5.0, 10.0]
        for y in y_values:
            p = model.inverse_transform(y)
            y_back = model.transform(p)
            assert np.isclose(y, y_back, rtol=1e-10), f"Roundtrip failed for y={y}"

    def test_linear_prediction_matches_transform(self):
        """Test that linear_prediction matches transform(predict)."""
        model = OurScurveModel(t=3, gamma=0.05)
        model.set_params(alpha=10.0, mu=50.0, beta=5.0)

        for w in range(5, 50):
            p = model.predict(w)
            if 0 < p < 0.5:
                y_from_predict = model.transform(p)
                y_from_linear = model.linear_prediction(w)
                assert np.isclose(y_from_predict, y_from_linear, rtol=1e-6), \
                    f"Mismatch at w={w}: transform(predict)={y_from_predict}, linear_prediction={y_from_linear}"

    def test_sweet_spot_in_valid_range(self):
        """Test that sweet spot is > t."""
        model = OurScurveModel(t=5, gamma=0.05)
        model.set_params(alpha=10.0, mu=50.0, beta=5.0)

        # Set internal params (a, b, c)
        model.set_internal_params(a=-0.1, b=5.0, c=5.0)

        w_sweet = model.calculate_sweet_spot()
        assert w_sweet > 5, f"Sweet spot {w_sweet} should be > t=5"

    def test_sweet_spot_varies_with_gamma(self):
        """Test that different gamma values give different sweet spots."""
        t = 5
        a, b, c = -0.1, 5.0, 5.0

        sweet_spots = []
        for gamma in [0.01, 0.05, 0.1, 0.2]:
            model = OurScurveModel(t=t, gamma=gamma)
            model.set_internal_params(a=a, b=b, c=c)
            sweet_spots.append(model.calculate_sweet_spot())

        # Sweet spots should generally decrease as gamma increases
        # (higher gamma means less emphasis on curvature)
        assert len(set(sweet_spots)) > 1, "Sweet spots should vary with gamma"

    def test_get_set_params(self):
        """Test parameter getting and setting."""
        model = OurScurveModel()

        model.set_params(alpha=15.0, mu=60.0, beta=8.0)
        params = model.get_params()

        assert params["alpha"] == 15.0
        assert params["mu"] == 60.0
        assert params["beta"] == 8.0

    def test_internal_params(self):
        """Test internal (a, b, c) parameter access."""
        model = OurScurveModel(t=3)
        model.set_internal_params(a=-0.05, b=3.0, c=4.0)

        internal = model.get_internal_params()
        assert internal["a"] == -0.05
        assert internal["b"] == 3.0
        assert internal["c"] == 4.0


class TestIBMScurveModel:
    """Tests for IBM S-curve model (Definition 2)."""

    def test_initialization(self):
        """Test model initialization."""
        model = IBMScurveModel()
        assert model.t == 0
        assert model.gamma == 0.05
        assert model.name == "IBMModel"

        model = IBMScurveModel(t=5, gamma=0.1)
        assert model.t == 5
        assert model.gamma == 0.1

    def test_param_names(self):
        """Test that param_names returns expected parameters."""
        model = IBMScurveModel()
        assert model.param_names == ("f0", "w0", "gamma_ibm")

    def test_predict_bounds(self):
        """Test that predictions are in valid range [0, 0.5)."""
        model = IBMScurveModel(t=3, gamma=0.05)
        model.set_params(f0=0.5, w0=10.0, gamma_ibm=1.5)

        for w in range(1, 100):
            p = model.predict(w)
            assert 0 <= p <= 0.5, f"P_L({w}) = {p} is out of bounds"

    def test_predict_array_input(self):
        """Test prediction with array input."""
        model = IBMScurveModel(t=2, gamma=0.05)
        model.set_params(f0=0.5, w0=10.0, gamma_ibm=1.5)

        weights = np.array([3, 5, 10, 20, 50])
        predictions = model.predict(weights)

        assert isinstance(predictions, np.ndarray)
        assert len(predictions) == len(weights)
        assert all(0 <= p <= 0.5 for p in predictions)

    def test_monotonicity(self):
        """Test that P_L increases with weight."""
        model = IBMScurveModel(t=0, gamma=0.05)
        model.set_params(f0=0.5, w0=10.0, gamma_ibm=1.5)

        weights = list(range(1, 50))
        predictions = [model.predict(w) for w in weights]

        for i in range(len(predictions) - 1):
            assert predictions[i] <= predictions[i + 1], \
                f"P_L should be monotonically increasing: P_L({weights[i]})={predictions[i]} > P_L({weights[i+1]})={predictions[i+1]}"

    def test_sweet_spot_positive(self):
        """Test that sweet spot is positive."""
        model = IBMScurveModel(t=3, gamma=0.05)
        model.set_params(f0=0.5, w0=10.0, gamma_ibm=1.5)

        w_sweet = model.calculate_sweet_spot()
        assert w_sweet > 0, f"Sweet spot {w_sweet} should be positive"

    def test_onset_weight(self):
        """Test get_onset_weight method."""
        model = IBMScurveModel(t=5)
        model.set_params(f0=0.5, w0=15.0, gamma_ibm=1.5)

        assert model.get_onset_weight() == 15.0


class TestModelFactory:
    """Tests for ModelFactory."""

    def test_create_our_model(self):
        """Test creating OurScurveModel."""
        model = ModelFactory.create(ModelType.OUR_MODEL, t=5, gamma=0.1)

        assert isinstance(model, OurScurveModel)
        assert model.t == 5
        assert model.gamma == 0.1

    def test_create_ibm_model(self):
        """Test creating IBMScurveModel."""
        model = ModelFactory.create(ModelType.IBM_MODEL, t=3, gamma=0.05)

        assert isinstance(model, IBMScurveModel)
        assert model.t == 3
        assert model.gamma == 0.05

    def test_create_all(self):
        """Test creating all models."""
        models = ModelFactory.create_all(t=5, gamma=0.1)

        assert ModelType.OUR_MODEL in models
        assert ModelType.IBM_MODEL in models
        assert isinstance(models[ModelType.OUR_MODEL], OurScurveModel)
        assert isinstance(models[ModelType.IBM_MODEL], IBMScurveModel)

    def test_available_models(self):
        """Test getting list of available models."""
        available = ModelFactory.available_models()

        assert ModelType.OUR_MODEL in available
        assert ModelType.IBM_MODEL in available

    def test_invalid_model_type(self):
        """Test that invalid model type raises error."""
        with pytest.raises(ValueError):
            ModelFactory.create("invalid_type")  # type: ignore


class TestModelFitting:
    """Tests for model fitting functionality."""

    @pytest.fixture
    def synthetic_data(self) -> tuple[List[int], List[float], List[int], List[int]]:
        """Generate synthetic S-curve data for testing."""
        t = 3
        # Generate weights
        weights = list(range(4, 30))

        # Generate P_L values following approximate S-curve
        p_values = []
        for w in weights:
            # Approximate S-curve: P_L increases with w
            x = (w - 10) / 5  # Center around w=10
            p = 0.5 / (1 + np.exp(-x))
            # Add small noise
            p = min(0.49, max(0.01, p + np.random.normal(0, 0.01)))
            p_values.append(p)

        # Generate sample counts and LE counts
        sample_counts = [1000] * len(weights)
        le_counts = [int(p * 1000) for p in p_values]

        return weights, p_values, sample_counts, le_counts

    def test_our_model_fit(self, synthetic_data):
        """Test that OurScurveModel can fit data."""
        weights, p_values, sample_counts, le_counts = synthetic_data

        model = OurScurveModel(t=3, gamma=0.05)
        model.fit(weights, p_values, sample_counts, le_counts)

        # Check that R² is reasonable
        assert model.r_squared > 0.5, f"R² = {model.r_squared} is too low"

        # Check that parameters were set
        params = model.get_params()
        assert "alpha" in params
        assert "mu" in params
        assert "beta" in params

    def test_ibm_model_fit(self, synthetic_data):
        """Test that IBMScurveModel can fit data."""
        weights, p_values, sample_counts, le_counts = synthetic_data

        model = IBMScurveModel(t=3, gamma=0.05)
        model.fit(weights, p_values, sample_counts, le_counts)

        # Check that R² is reasonable (IBM model may not fit as well for this data)
        assert model.r_squared >= 0.0

        # Check that parameters were set
        params = model.get_params()
        assert "f0" in params
        assert "w0" in params
        assert "gamma_ibm" in params

    def test_fit_updates_r_squared(self, synthetic_data):
        """Test that fitting updates R² score."""
        weights, p_values, sample_counts, le_counts = synthetic_data

        model = OurScurveModel(t=3, gamma=0.05)
        assert model.r_squared == 0.0  # Initial value

        model.fit(weights, p_values, sample_counts, le_counts)
        assert model.r_squared != 0.0  # Updated after fit

    def test_fit_with_insufficient_data(self):
        """Test fitting with insufficient data."""
        model = OurScurveModel(t=3, gamma=0.05)

        # Only one data point
        model.fit([5], [0.1], [100], [10])

        # Should not crash, R² should remain 0
        assert model.r_squared == 0.0


class TestStatisticalEstimators:
    """Tests for statistical estimator methods in base class."""

    def test_sigma_estimator_valid_input(self):
        """Test sigma estimator with valid input."""
        # n=1000 samples, m=100 errors
        sigma = ScurveModelBase._sigma_estimator(1000, 100)
        assert sigma > 0
        assert np.isfinite(sigma)

    def test_sigma_estimator_edge_cases(self):
        """Test sigma estimator with edge cases."""
        # Zero errors
        sigma = ScurveModelBase._sigma_estimator(1000, 0)
        assert sigma == 1.0  # Default fallback

        # Equal errors to samples (invalid for this formula)
        sigma = ScurveModelBase._sigma_estimator(100, 50)
        assert sigma == 1.0  # n <= 2*m

    def test_bias_estimator_valid_input(self):
        """Test bias estimator with valid input."""
        bias = ScurveModelBase._bias_estimator(1000, 100)
        assert np.isfinite(bias)

    def test_bias_estimator_edge_cases(self):
        """Test bias estimator with edge cases."""
        # Zero errors
        bias = ScurveModelBase._bias_estimator(1000, 0)
        assert bias == 0.0

    def test_r_squared_score(self):
        """Test R² calculation."""
        y_true = [1.0, 2.0, 3.0, 4.0, 5.0]
        y_pred = [1.1, 2.0, 2.9, 4.1, 4.9]

        r2 = ScurveModelBase._r_squared_score(y_true, y_pred)
        assert 0 <= r2 <= 1
        assert r2 > 0.9  # Should be high for this good fit

    def test_r_squared_perfect_fit(self):
        """Test R² for perfect fit."""
        y_true = [1.0, 2.0, 3.0, 4.0, 5.0]
        y_pred = y_true.copy()

        r2 = ScurveModelBase._r_squared_score(y_true, y_pred)
        assert r2 == 1.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
