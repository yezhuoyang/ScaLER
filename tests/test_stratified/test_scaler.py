"""
Test suite for the Scaler class implementation.
"""
import os
import sys
import pytest
import numpy as np

# Add the src directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))

from scalerqec.Stratified.Scaler import Scaler


class TestScalerBasic:
    """Basic tests for Scaler functionality."""

    def test_scaler_initialization(self):
        """Test that Scaler initializes with correct default values."""
        scaler = Scaler(error_rate=0.001, time_budget=30)
        assert scaler._error_rate == 0.001
        assert scaler._time_budget == 30.0
        assert scaler._remaining_time_budget == 30.0
        assert scaler._sampling_rate == 0.0
        assert scaler._LER == 0.0

    def test_scaler_parse_surface3(self):
        """Test parsing a small surface code circuit."""
        filepath = "c:/Users/yezhu/GitRepos/ScaLERQEC/stimprograms/surface/surface3"
        scaler = Scaler(error_rate=0.001, time_budget=30)
        scaler.parse_from_file(filepath)

        assert scaler._num_noise > 0
        assert scaler._num_detector > 0
        assert scaler._matcher is not None
        assert scaler._QEPG_graph is not None

    def test_calc_logical_error_rate_with_fixed_w(self):
        """Test calculating logical error rate for a fixed weight."""
        filepath = "c:/Users/yezhu/GitRepos/ScaLERQEC/stimprograms/surface/surface3"
        scaler = Scaler(error_rate=0.001, time_budget=30)
        scaler.parse_from_file(filepath)

        # Test with a few different weights
        ler1 = scaler.calc_logical_error_rate_with_fixed_w(shots=1000, w=1)
        ler2 = scaler.calc_logical_error_rate_with_fixed_w(shots=1000, w=5)

        # Error rate should be between 0 and 1
        assert 0 <= ler1 <= 1
        assert 0 <= ler2 <= 1
        # Higher weight should generally have higher error rate (though not guaranteed with sampling)
        # We just check they're reasonable values


class TestScalerSampling:
    """Tests for sampling functionality."""

    def test_measure_sample_rates(self):
        """Test that sample rate measurement works."""
        filepath = "c:/Users/yezhu/GitRepos/ScaLERQEC/stimprograms/surface/surface3"
        scaler = Scaler(error_rate=0.001, time_budget=30)
        scaler.parse_from_file(filepath)

        initial_budget = scaler._remaining_time_budget
        scaler.measure_sample_rates()

        # Sample rate should be positive
        assert scaler._sampling_rate > 0
        # Time budget should have decreased
        assert scaler._remaining_time_budget < initial_budget

    def test_sampling_step(self):
        """Test a single sampling step."""
        filepath = "c:/Users/yezhu/GitRepos/ScaLERQEC/stimprograms/surface/surface3"
        scaler = Scaler(error_rate=0.001, time_budget=30)
        scaler.parse_from_file(filepath)
        scaler.measure_sample_rates()

        wlist = [5, 10]
        slist = [100, 100]

        initial_budget = scaler._remaining_time_budget
        elapsed = scaler._sampling_step(wlist, slist)

        # Should have elapsed some time
        assert elapsed > 0
        # Subspace stats should be updated
        assert 5 in scaler._subspace_sample_used
        assert 10 in scaler._subspace_sample_used
        assert scaler._subspace_sample_used[5] == 100
        assert scaler._subspace_sample_used[10] == 100


class TestScalerFitting:
    """Tests for S-curve fitting functionality."""

    def test_fit_log_S_model_with_data(self):
        """Test S-curve fitting with some sample data."""
        filepath = "c:/Users/yezhu/GitRepos/ScaLERQEC/stimprograms/surface/surface3"
        scaler = Scaler(error_rate=0.001, time_budget=60)
        scaler._circuit_level_code_distance = 3
        scaler._t = 1
        scaler.parse_from_file(filepath)
        scaler.measure_sample_rates()

        # Generate some sample data
        scaler.determine_lower_w()
        scaler.determine_saturated_w()

        # Sample a few weights
        wlist = [scaler._has_logical_errorw + i * 2
                 for i in range(5)
                 if scaler._has_logical_errorw + i * 2 <= scaler._saturatew]
        slist = [1000] * len(wlist)

        if wlist:
            scaler._sampling_step(wlist, slist)
            scaler.fit_log_S_model(savefigure=False)

            # Check that parameters were updated
            assert scaler._a != 0.0 or scaler._b != 0.0
            assert scaler._sweet_spot is not None


class TestScalerIntegration:
    """Integration tests for the full workflow."""

    @pytest.mark.slow
    def test_calculate_LER_surface3(self):
        """Test full LER calculation on surface3."""
        filepath = "c:/Users/yezhu/GitRepos/ScaLERQEC/stimprograms/surface/surface3"
        scaler = Scaler(error_rate=0.001, time_budget=60)

        ler = scaler.calculate_LER_from_file(
            filepath=filepath,
            pvalue=0.001,
            codedistance=3,
            figname="test_surface3_",
            titlename="Test Surface 3",
            repeat=1
        )

        # LER should be a reasonable value
        assert ler is not None
        assert 0 < ler < 1
        assert scaler._R_square_score >= 0

    @pytest.mark.slow
    def test_calculate_LER_surface5(self):
        """Test full LER calculation on surface5."""
        filepath = "c:/Users/yezhu/GitRepos/ScaLERQEC/stimprograms/surface/surface5"
        scaler = Scaler(error_rate=0.001, time_budget=120)

        ler = scaler.calculate_LER_from_file(
            filepath=filepath,
            pvalue=0.001,
            codedistance=5,
            figname="test_surface5_",
            titlename="Test Surface 5",
            repeat=1
        )

        # LER should be a reasonable value
        assert ler is not None
        assert 0 < ler < 1


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v", "-s"])
