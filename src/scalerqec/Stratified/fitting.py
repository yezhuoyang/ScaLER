"""Goodness-of-fit statistics for S-curve model evaluation.

This module provides the coefficient of determination (R^2) used to
assess how well a fitted S-curve model explains the measured P_L(w)
data in transformed (log-logit) space.
"""

from __future__ import annotations

from typing import Sequence
import numpy as np
import numpy.typing as npt

FloatArray = npt.NDArray[np.float64]


def r_squared(
    y_true: Sequence[float], y_pred: Sequence[float], clip: bool = False
) -> float:
    """Compute the coefficient of determination (R^2).

    R^2 measures the fraction of variance in the observed values that
    is explained by the model predictions.  A value of 1.0 indicates a
    perfect fit; values below 0 indicate the model is worse than a
    horizontal line at the mean.

    Args:
        y_true: Observed (true) values.
        y_pred: Model-predicted values, same length as *y_true*.
        clip: If ``True``, negative R^2 values are clipped to 0 so the
            returned score lies in [0, 1].

    Returns:
        The R^2 statistic.

    Raises:
        ValueError: If *y_true* and *y_pred* have different shapes.
    """
    yt: FloatArray = np.asarray(y_true, dtype=np.float64)
    yp: FloatArray = np.asarray(y_pred, dtype=np.float64)

    if yt.shape != yp.shape:
        raise ValueError("y_true and y_pred must have the same shape")

    ss_res = np.sum((yt - yp) ** 2)  # residual sum of squares
    ss_tot = np.sum((yt - yt.mean()) ** 2)  # total sum of squares

    # Handle the degenerate case where variance is zero
    if ss_tot == 0.0:
        return 1.0 if ss_res == 0.0 else 0.0

    r2 = 1.0 - ss_res / ss_tot
    return max(0.0, r2) if clip else r2
