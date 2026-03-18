"""Formatting utilities for presenting numerical results with uncertainty.

Provides functions for displaying scientific results in a compact,
human-readable notation suitable for logging and reports.
"""

import numpy as np


def format_with_uncertainty(value: float, std: float) -> str:
    r"""Format a value with its uncertainty in scientific notation.

    Produces a compact string of the form ``coeff(+std_coeff)*10^exp``,
    e.g., ``1.23(+0.45)*10^-3``.

    For zero values, returns ``0(+<std in scientific notation>)``.

    Args:
        value (float): The central value to format.
        std (float): The standard deviation (uncertainty).

    Returns:
        str: Formatted string with coefficient, uncertainty, and
        power-of-ten exponent.

    Example::

        >>> format_with_uncertainty(0.00123, 0.00045)
        '1.23(+0.45)*10^-3'
    """
    if value == 0:
        return f"0(+{std:.2e})"
    exponent = int(np.floor(np.log10(abs(value))))
    coeff = value / (10**exponent)
    std_coeff = std / (10**exponent)
    return f"{coeff:.2f}(+{std_coeff:.2f})*10^{exponent}"
