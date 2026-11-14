"""Exercise the constrained CPU ABI with distinct bounds and custom coordinates."""

import numpy as np
import pycpufit.cpufit as cf


def test_distinct_bounds_per_fit():
    x = np.linspace(2, 6, 20, dtype=np.float32)
    data = np.tile(3 * x, (3, 1))
    starts = np.tile(np.array([0, 0.5], np.float32), (3, 1))
    bounds = np.array([[0, 0, 0, 1], [0, 0, 0, 2], [0, 0, 0, 0.75]], np.float32)
    result = cf.fit_constrained(
        data,
        None,
        cf.ModelID.LINEAR_1D,
        starts,
        bounds,
        np.full(2, cf.ConstraintType.LOWER_UPPER, np.int32),
        parameters_to_fit=np.array([0, 1], np.int32),
        user_info=x,
    )
    np.testing.assert_allclose(result[0][:, 0], 0, atol=1e-6)
    np.testing.assert_allclose(result[0][:, 1], bounds[:, 3], atol=2e-5)
