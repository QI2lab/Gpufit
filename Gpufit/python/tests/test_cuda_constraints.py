"""Numerical acceptance through the compiled native optimizer on an actual device."""

import os
import subprocess
import sys

import numpy as np
import pygpufit.gpufit as gf
import pytest

pytestmark = pytest.mark.gpu


@pytest.mark.parametrize("model_id", [13, 14, 15, 16])
def test_custom_coordinates_across_internal_chunks(model_id, request):
    """Every chunk uses its own coordinates and the full batch's size table."""
    if _chunk_subprocess(request, "2"):
        return
    z, y, x = np.meshgrid(
        np.linspace(-2, 2, 7), np.linspace(-1, 1, 9), np.linspace(-1, 1, 9), indexing="ij"
    )
    coords = np.stack([x.ravel(), y.ravel(), z.ravel()])
    n = x.size
    centers = np.array([[0, 0, 0], [3, -2, 1], [-4, 5, -3]], np.float32)
    params = np.array(
        [
            [100, *c, 0.4, 0.8, 5]
            if model_id in (13, 14)
            else [100, *c, 0.4, 0.5, 0.8, 5, *([0, 0, 0] if model_id == 15 else [])]
            for c in centers
        ],
        np.float32,
    )
    start = params.copy()
    start[:, 1:4] += 0.05
    relative = coords / np.array([0.4, 0.4 if model_id in (13, 14) else 0.5, 0.8])[:, None]
    image = 5 + 100 * np.exp(-0.5 * np.sum(relative**2, axis=0))
    if model_id == 14:
        lorentz = 1 + relative[2] ** 2
        image = 5 + 100 * np.exp(-0.5 * np.sum(relative[:2] ** 2, axis=0) / lorentz) / lorentz
    info = np.concatenate(
        [
            (coords[None] + centers[:, :, None]).ravel(),
            np.full(3, n),
            np.zeros(0 if model_id == 14 else (2 if model_id == 13 else 3)),
        ]
    ).astype(np.float32)
    # Fix widths/angles to avoid Euler angle degeneracy at zero rotation.
    mask = np.zeros(params.shape[1], np.int32)
    mask[:4] = 1
    result = gf.fit(
        np.tile(image, (3, 1)).astype(np.float32),
        None,
        model_id,
        start,
        user_info=info,
        parameters_to_fit=mask,
        tolerance=1e-6,
    )
    np.testing.assert_array_equal(result[1], 0)
    np.testing.assert_allclose(result[0], params, atol=3e-4, rtol=3e-4)


def test_invalid_chunk_limit_rejected(request):
    if _chunk_subprocess(request, "invalid"):
        return
    with pytest.raises(RuntimeError, match="GPUFIT_MAX_CHUNK_SIZE"):
        gf.fit(np.ones((1, 3), np.float32), None, gf.ModelID.LINEAR_1D, np.ones((1, 2), np.float32))


def _chunk_subprocess(request, limit):
    # The statically linked Windows CRT reads its environment when the DLL loads.
    # Set the option before process startup rather than monkeypatching afterward.
    if os.environ.get("GPUFIT_CHUNK_TEST_CHILD") == "1":
        return False
    environment = dict(os.environ, GPUFIT_MAX_CHUNK_SIZE=limit, GPUFIT_CHUNK_TEST_CHILD="1")
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", request.node.nodeid, "--tb=short"],
        env=environment,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return True


def test_device_binary_identity():
    assert gf.cuda_available(), gf.get_last_error()
    identity = gf.get_capabilities()
    assert identity["use_cublas"] is True
    assert identity["precision"] == "float32"
    assert identity["models"]["GAUSS_3D_ROT_ARB"] == {"id": 15, "parameters": 11}


def test_distinct_bounds_per_fit_across_chunks(request):
    if _chunk_subprocess(request, "2"):
        return
    x = np.arange(20, dtype=np.float32)
    data = np.tile(3 * x, (3, 1))
    starts = np.tile(np.array([0, 0.5], np.float32), (3, 1))
    bounds = np.array([[0, 0, 0, 1], [0, 0, 0, 2], [0, 0, 0, 0.75]], np.float32)
    result = gf.fit_constrained(
        data,
        None,
        gf.ModelID.LINEAR_1D,
        starts,
        bounds,
        np.full(2, 3, np.int32),
        parameters_to_fit=np.array([0, 1], np.int32),
        user_info=x,
    )
    np.testing.assert_allclose(result[0][:, 1], bounds[:, 3], atol=2e-5)


def test_actual_bounds_fixed_parameters_and_weighted_objective():
    x = np.arange(30, dtype=np.float32)
    data = (2 + 3 * x)[None, :]
    start = np.array([[2, 0.5]], np.float32)
    bounds = np.array([[2, 2, 0, 1]], np.float32)
    weights = np.linspace(0.25, 2, x.size, dtype=np.float32)[None, :]
    result = gf.fit_constrained(
        data,
        weights,
        gf.ModelID.LINEAR_1D,
        start,
        bounds,
        np.array([3, 3], np.int32),
        parameters_to_fit=np.array([0, 1], np.int32),
        user_info=x,
        max_number_iterations=100,
    )
    np.testing.assert_allclose(result[0], [[2, 1]], atol=2e-5)
    predicted = result[0][0, 0] + result[0][0, 1] * x
    np.testing.assert_allclose(
        result[2][0], np.sum(weights[0] * (predicted - data[0]) ** 2), rtol=1e-5
    )


@pytest.mark.parametrize("model_id", [13, 14, 15, 16])
@pytest.mark.parametrize("constrained", [False, True])
def test_arbitrary_coordinate_3d_recovery(model_id, constrained):
    z, y, x = np.meshgrid(
        np.linspace(-2, 2, 9), np.linspace(-1.5, 1.5, 13), np.linspace(-1.5, 1.5, 13), indexing="ij"
    )
    center = [0.13, -0.17, 0.24]  # XYZ parameter convention
    sigma = [0.42, 0.42, 0.83] if model_id in (13, 14) else [0.42, 0.53, 0.83]
    image = 5 + 100 * np.exp(
        -0.5
        * (
            ((x - center[0]) / sigma[0]) ** 2
            + ((y - center[1]) / sigma[1]) ** 2
            + ((z - center[2]) / sigma[2]) ** 2
        )
    )
    if model_id == 14:
        lorentz = 1 + ((z - center[2]) / sigma[2]) ** 2
        image = (
            5
            + 100
            * np.exp(
                -0.5
                * (((x - center[0]) / sigma[0]) ** 2 + ((y - center[1]) / sigma[1]) ** 2)
                / lorentz
            )
            / lorentz
        )
    angles = [0.2, 0.3, 0.4]
    if model_id == 15:

        def rz(angle):
            c, s = np.cos(angle), np.sin(angle)
            return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

        c, s = np.cos(angles[1]), np.sin(angles[1])
        ry = np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
        inverse = (rz(angles[0]) @ ry @ rz(angles[2])).T
        relative = inverse @ (
            np.stack([x.ravel(), y.ravel(), z.ravel()]) - np.array(center)[:, None]
        )
        image = (
            5 + 100 * np.exp(-0.5 * np.sum((relative / np.array(sigma)[:, None]) ** 2, axis=0))
        ).reshape(x.shape)
    true = np.array(
        [100, *center, sigma[0], sigma[2], 5]
        if model_id in (13, 14)
        else [100, *center, *sigma, 5, *(angles if model_id == 15 else [])],
        np.float32,
    )
    start = true.copy()
    start[:4] += [-10, 0.05, -0.04, 0.06]
    start = start[None, :]
    info = np.concatenate(
        [
            np.stack([x.ravel(), y.ravel(), z.ravel()]).ravel(),
            [x.size],
            np.zeros(0 if model_id == 14 else (2 if model_id == 13 else 3)),
        ]
    ).astype(np.float32)
    data = image.ravel()[None, :].astype(np.float32)
    kwargs = {"user_info": info, "tolerance": 1e-6, "max_number_iterations": 100}
    if constrained:
        lower = true - np.maximum(1, abs(true))
        upper = true + np.maximum(1, abs(true))
        lower[4 : 6 if model_id in (13, 14) else 7] = 0.05
        bounds = np.stack([lower, upper], axis=-1).reshape(1, -1).astype(np.float32)
        result = gf.fit_constrained(
            data, None, model_id, start, bounds, np.full(true.size, 3, np.int32), **kwargs
        )
    else:
        result = gf.fit(data, None, model_id, start, **kwargs)
    assert result[1][0] == 0
    np.testing.assert_allclose(result[0][0], true, rtol=2e-4, atol=2e-4)
