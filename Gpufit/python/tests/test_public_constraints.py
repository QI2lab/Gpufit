"""Public ABI regression tests; numerical tests require the compiled QI2lab fork."""

import ctypes
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest


class NativeFunction:
    def __init__(self, result=0):
        self.result = result
        self.calls = []

    def __call__(self, *args):
        self.calls.append(args)
        return self.result


@pytest.fixture
def binding(monkeypatch):
    from types import SimpleNamespace

    library = SimpleNamespace(
        **{
            name: NativeFunction()
            for name in (
                "gpufit",
                "gpufit_constrained",
                "gpufit_cuda_available",
                "gpufit_get_cuda_version",
                "gpufit_get_last_error",
            )
        }
    )
    monkeypatch.setattr(ctypes.cdll, "LoadLibrary", lambda path: library)
    path = Path(__file__).parents[1] / "pygpufit" / "gpufit.py"
    spec = importlib.util.spec_from_file_location("binding_under_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    library.gpufit_get_capabilities = NativeFunction(
        json.dumps(
            {
                "abi_version": 3,
                "precision": "float32",
                "models": module.MICROSCOPY_MODELS,
                "constraints_layout": "fit-major-interleaved-lower-upper",
            }
        ).encode()
    )
    return module, library


def arrays():
    return np.ones((2, 5), np.float32), np.ones((2, 2), np.float32)


def test_constrained_pointer_order_and_forwarding(binding):
    gf, native = binding
    data, start = arrays()
    weights = np.full_like(data, 2)
    constraints = np.tile([0, 3, -1, 4], (2, 1)).astype(np.float32)
    types = np.array([1, 3], np.int32)
    mask = np.array([0, 1], np.int32)
    info = np.array([3, 7], np.float32)
    result = gf.fit_constrained(
        data,
        weights,
        gf.ModelID.LINEAR_1D,
        start,
        constraints,
        types,
        tolerance=0.003,
        max_number_iterations=87,
        parameters_to_fit=mask,
        estimator_id=gf.EstimatorID.LSE,
        user_info=info,
    )
    args = native.gpufit_constrained.calls[-1]
    assert len(args) == 18
    assert np.ctypeslib.as_array(args[6], shape=(8,)).tolist() == constraints.ravel().tolist()
    assert np.ctypeslib.as_array(args[7], shape=(2,)).tolist() == [1, 3]
    assert args[8].value == pytest.approx(0.003)
    assert args[9].value == 87
    assert np.ctypeslib.as_array(args[10], shape=(2,)).tolist() == [0, 1]
    assert args[12].value == info.nbytes
    assert result[0].shape == (2, 2)
    assert result[1].dtype == np.int32


def test_unbounded_keeps_original_native_dispatch(binding):
    gf, native = binding
    gf.fit(*arrays()[:1], None, gf.ModelID.LINEAR_1D, arrays()[1])
    assert len(native.gpufit.calls) == 1
    assert not native.gpufit_constrained.calls


@pytest.mark.parametrize(
    "change", ["shape", "dtype", "types", "mask", "reversed", "infeasible", "nan"]
)
def test_malformed_inputs_rejected_before_native_call(binding, change):
    gf, native = binding
    data, start = arrays()
    bounds = np.tile([0, 3, 0, 3], (2, 1)).astype(np.float32)
    types = np.array([3, 3], np.int32)
    mask = np.array([1, 1], np.int32)
    if change == "shape":
        bounds = bounds[:, :3]
    if change == "dtype":
        bounds = bounds.astype(np.float64)
    if change == "types":
        types[0] = 5
    if change == "mask":
        mask[0] = 2
    if change == "reversed":
        bounds[0, :2] = [3, 0]
    if change == "infeasible":
        start[0, 0] = 9
    if change == "nan":
        bounds[0, 0] = np.nan
    with pytest.raises((ValueError, RuntimeError, TypeError)):
        gf.fit_constrained(
            data, None, gf.ModelID.LINEAR_1D, start, bounds, types, parameters_to_fit=mask
        )
    assert not native.gpufit_constrained.calls


@pytest.mark.parametrize("abi_version", [1, 2])
def test_reject_incompatible_native_binary(binding, abi_version):
    gf, native = binding
    manifest = json.loads(native.gpufit_get_capabilities.result)
    manifest["abi_version"] = abi_version
    native.gpufit_get_capabilities.result = json.dumps(manifest).encode()
    with pytest.raises(RuntimeError, match="QI2lab|capabilit"):
        gf.get_capabilities()


@pytest.mark.parametrize("model_id, n_parameters, dimensions", [(13, 6, 4), (14, 7, 5)])
def test_upstream_spline_ids_accept_spline_user_info(binding, model_id, n_parameters, dimensions):
    gf, native = binding
    # Upstream splines must bypass the arbitrary-coordinate microscopy layout.
    del native.gpufit_get_capabilities
    info = np.concatenate([np.ones(2 * dimensions), np.zeros(4**dimensions)]).astype(np.float32)
    gf.fit(
        np.ones((2, 16), np.float32),
        None,
        model_id,
        np.ones((2, n_parameters), np.float32),
        user_info=info,
    )
    args = native.gpufit.calls[-1]
    assert args[4].value == model_id
    assert args[10].value == info.nbytes


def test_microscopy_rejects_previous_native_model_ids(binding):
    gf, native = binding
    manifest = json.loads(native.gpufit_get_capabilities.result)
    for model in manifest["models"].values():
        model["id"] -= 2
    native.gpufit_get_capabilities.result = json.dumps(manifest).encode()
    with pytest.raises(RuntimeError, match="Incompatible"):
        gf.get_capabilities()
