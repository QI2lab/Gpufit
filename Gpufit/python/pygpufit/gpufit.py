"""
Python binding for Gpufit, a Levenberg Marquardt curve fitting library written in CUDA
See https://github.com/gpufit/Gpufit, http://gpufit.readthedocs.io/en/latest/bindings.html#python

The binding is based on ctypes.
See https://docs.python.org/3.5/library/ctypes.html, http://www.scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html
"""

import json
import os
import time
from ctypes import POINTER, byref, c_char, c_char_p, c_float, c_int, c_size_t, cdll
from importlib.metadata import PackageNotFoundError, distribution

import numpy as np

# define library loader (actual loading is lazy)
package_dir = os.path.dirname(os.path.realpath(__file__))

if os.name == "nt":
    lib_path = os.path.join(package_dir, "Gpufit.dll")  # library name on Windows
elif os.name == "posix":
    lib_path = os.path.join(package_dir, "libGpufit.so")  # library name on Unix
else:
    raise RuntimeError(f"OS {os.name} not supported by pyGpufit.")

# Native editable installs keep Python in the checkout and binaries in site-packages.
if not os.path.isfile(lib_path):
    try:
        installed_path = distribution("pyGpufit").locate_file(
            "pygpufit/" + os.path.basename(lib_path)
        )
        if installed_path.is_file():
            lib_path = str(installed_path)
            package_dir = str(installed_path.parent)
    except PackageNotFoundError:
        pass

_dll_directory = os.add_dll_directory(package_dir) if os.name == "nt" else None
lib = cdll.LoadLibrary(lib_path)

# gpufit function in the dll
gpufit_func = lib.gpufit
gpufit_func.restype = c_int
gpufit_func.argtypes = [
    c_size_t,
    c_size_t,
    POINTER(c_float),
    POINTER(c_float),
    c_int,
    POINTER(c_float),
    c_float,
    c_int,
    POINTER(c_int),
    c_int,
    c_size_t,
    POINTER(c_char),
    POINTER(c_float),
    POINTER(c_int),
    POINTER(c_float),
    POINTER(c_int),
]

constrained_func = lib.gpufit_constrained
constrained_func.restype = c_int
constrained_func.argtypes = (
    gpufit_func.argtypes[:6] + [POINTER(c_float), POINTER(c_int)] + gpufit_func.argtypes[6:]
)

# gpufit_get_last_error function in the dll
error_func = lib.gpufit_get_last_error
error_func.restype = c_char_p
error_func.argtypes = None

# gpufit_cuda_available function in the dll
cuda_available_func = lib.gpufit_cuda_available
cuda_available_func.restype = c_int
cuda_available_func.argtypes = None

# gpufit_get_cuda_version function in the dll
get_cuda_version_func = lib.gpufit_get_cuda_version
get_cuda_version_func.restype = c_int
get_cuda_version_func.argtypes = [POINTER(c_int), POINTER(c_int)]


class ModelID:
    GAUSS_1D = 0
    GAUSS_2D = 1
    GAUSS_2D_ELLIPTIC = 2
    GAUSS_2D_ROTATED = 3
    CAUCHY_2D_ELLIPTIC = 4
    LINEAR_1D = 5
    FLETCHER_POWELL_HELIX = 6
    BROWN_DENNIS = 7
    SPLINE_1D = 8
    SPLINE_2D = 9
    SPLINE_3D = 10
    SPLINE_3D_MULTICHANNEL = 11
    SPLINE_3D_PHASE_MULTICHANNEL = 12
    GAUSS_3D_ARB = 13
    GAUSS_LOR_3D_ARB = 14
    GAUSS_3D_ROT_ARB = 15
    GAUSS_3D_ASYM_ARB = 16


class EstimatorID:
    LSE = 0
    MLE = 1


class Status:
    Ok = 0
    Error = -1


class ConstraintType:
    """Native constraint enumeration; patterns are shared across a fit batch."""

    FREE = 0
    NONE = 0
    LOWER = 1
    UPPER = 2
    LOWER_UPPER = 3


MODEL_PARAMETER_COUNTS = dict(enumerate((4, 5, 6, 7, 6, 2, 3, 4, 3, 4, 5, 5, 6, 7, 7, 11, 8)))
MICROSCOPY_MODELS = {
    "GAUSS_3D_ARB": {"id": 13, "parameters": 7},
    "GAUSS_LOR_3D_ARB": {"id": 14, "parameters": 7},
    "GAUSS_3D_ROT_ARB": {"id": 15, "parameters": 11},
    "GAUSS_3D_ASYM_ARB": {"id": 16, "parameters": 8},
}


def get_capabilities():
    """Validate the loaded binary's microscopy model/precision ABI and return its manifest."""
    if not hasattr(lib, "gpufit_get_capabilities"):
        raise RuntimeError("QI2lab GPUfit native capability export missing; rebuild the fork")
    function = lib.gpufit_get_capabilities
    function.restype = c_char_p
    function.argtypes = []
    try:
        manifest = json.loads(function().decode("utf-8"))
    except (AttributeError, TypeError, ValueError) as exc:
        raise RuntimeError("Invalid QI2lab GPUfit capability manifest") from exc
    if (
        manifest.get("abi_version") != 2
        or manifest.get("precision") != "float32"
        or manifest.get("models") != MICROSCOPY_MODELS
        or manifest.get("constraints_layout") != "fit-major-interleaved-lower-upper"
    ):
        raise RuntimeError("Incompatible QI2lab GPUfit model/precision ABI")
    manifest["library_path"] = lib_path
    return manifest


def fit(
    data,
    weights,
    model_id,
    initial_parameters,
    tolerance=None,
    max_number_iterations=None,
    parameters_to_fit=None,
    estimator_id=None,
    user_info=None,
):
    """Fit host float32 buffers using the original unconstrained native entry point.

    Returns
    -------
    tuple
        Parameters, states, objectives, iteration counts, and elapsed seconds.
    """
    return _fit(
        data,
        weights,
        model_id,
        initial_parameters,
        tolerance,
        max_number_iterations,
        parameters_to_fit,
        estimator_id,
        user_info,
    )


def fit_constrained(
    data,
    weights,
    model_id,
    initial_parameters,
    constraints=None,
    constraint_types=None,
    tolerance=None,
    max_number_iterations=None,
    parameters_to_fit=None,
    estimator_id=None,
    user_info=None,
):
    """Fit through the supported native constrained interface.

    Parameters
    ----------
    constraints : numpy.ndarray, optional
        Float32 shape (nfits, 2*nparams), interleaved lower/upper endpoints.
        Active endpoints are finite; unused endpoints are ignored by the optimizer.
    constraint_types : numpy.ndarray, optional
        Int32 shape (nparams,), from ConstraintType, shared across the batch.
        Group fits when constraint patterns differ.
    data, weights, model_id, initial_parameters, tolerance, max_number_iterations,
    parameters_to_fit, estimator_id, user_info
        Same host-buffer arguments and controls as fit. Weights are inverse variances
        for LSE. EstimatorID.MLE selects the native Poisson objective.

    Returns
    -------
    tuple
        Parameters, states, objectives, iteration counts, and elapsed seconds.
    """
    return _fit(
        data,
        weights,
        model_id,
        initial_parameters,
        tolerance,
        max_number_iterations,
        parameters_to_fit,
        estimator_id,
        user_info,
        constraints,
        constraint_types,
        True,
    )


def _fit(
    data,
    weights,
    model_id,
    initial_parameters,
    tolerance=None,
    max_number_iterations=None,
    parameters_to_fit=None,
    estimator_id=None,
    user_info=None,
    constraints=None,
    constraint_types=None,
    constrained=False,
):
    """
    Calls the C interface fit function in the library.
    (see also http://gpufit.readthedocs.io/en/latest/bindings.html#python)

    All 2D NumPy arrays must be in row-major order (standard in NumPy), i.e. array.flags.C_CONTIGUOUS must be True
    (see also https://docs.scipy.org/doc/numpy/reference/arrays.ndarray.html#internal-memory-layout-of-an-ndarray)

    :param data: The data - 2D NumPy array of dimension [number_fits, number_points] and data type np.float32
    :param weights: The weights - 2D NumPy array of the same dimension and data type as parameter data or None (no weights available)
    :param model_id: The model ID
    :param initial_parameters: Initial values for parameters - NumPy array of dimension [number_fits, number_parameters] and data type np.float32
    :param tolerance: The fit tolerance or None (will use default value)
    :param max_number_iterations: The maximal number of iterations or None (will use default value)
    :param parameters_to_fit: Which parameters to fit - NumPy array of length number_parameters and type np.int32 or None (will fit all parameters)
    :param estimator_id: The Estimator ID or None (will use default values)
    :param user_info: User info - NumPy array of type np.char or None (no user info available)
    :return: parameters, states, chi_squares, number_iterations, execution_time
    """

    # check all 2D NumPy arrays for row-major memory layout (otherwise interpretation of order of dimensions fails)
    if not data.flags.c_contiguous:
        raise RuntimeError("Memory layout of data array mismatch.")

    if weights is not None and not weights.flags.c_contiguous:
        raise RuntimeError("Memory layout of weights array mismatch.")

    if not initial_parameters.flags.c_contiguous:
        raise RuntimeError("Memory layout of initial_parameters array mismatch.")

    # size check: data is 2D and read number of points and fits
    if data.ndim != 2:
        raise RuntimeError("data is not two-dimensional")
    number_points = data.shape[1]
    number_fits = data.shape[0]

    # size check: consistency with weights (if given)
    if weights is not None and data.shape != weights.shape:
        raise RuntimeError("dimension mismatch between data and weights")
        # the unequal operator checks, type, length and content (https://docs.python.org/3.7/reference/expressions.html#value-comparisons)

    # size check: initial parameters is 2D and read number of parameters
    if initial_parameters.ndim != 2:
        raise RuntimeError("initial_parameters is not two-dimensional")
    number_parameters = initial_parameters.shape[1]
    if initial_parameters.shape[0] != number_fits:
        raise RuntimeError(
            "dimension mismatch in number of fits between data and initial_parameters"
        )

    # size check: consistency with parameters_to_fit (if given)
    if parameters_to_fit is not None and parameters_to_fit.shape[0] != number_parameters:
        raise RuntimeError(
            "dimension mismatch in number of parameters between initial_parameters and parameters_to_fit"
        )

    # default value: tolerance
    if tolerance is None:
        tolerance = 1e-4

    # default value: max_number_iterations
    if max_number_iterations is None:
        max_number_iterations = 25

    # default value: estimator ID
    if estimator_id is None:
        estimator_id = EstimatorID.LSE

    # default value: parameters_to_fit
    if parameters_to_fit is None:
        parameters_to_fit = np.ones(number_parameters, dtype=np.int32)

    # now only weights and user_info could be not given

    # type check: data, weights (if given), initial_parameters are all np.float32
    if data.dtype != np.float32:
        raise RuntimeError("type of data is not np.float32")
    if weights is not None and weights.dtype != np.float32:
        raise RuntimeError("type of weights is not np.float32")
    if initial_parameters.dtype != np.float32:
        raise RuntimeError("type of initial_parameters is not np.float32")

    # type check: parameters_to_fit is np.int32
    if parameters_to_fit.dtype != np.int32:
        raise RuntimeError("type of parameters_to_fit is not np.int32")

    if (
        model_id not in MODEL_PARAMETER_COUNTS
        or number_parameters != MODEL_PARAMETER_COUNTS[model_id]
    ):
        raise ValueError("Invalid model ID or parameter count")
    if estimator_id not in (EstimatorID.LSE, EstimatorID.MLE):
        raise ValueError("Invalid estimator ID")
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be positive and finite")
    if (
        not isinstance(max_number_iterations, (int, np.integer))
        or not 1 <= max_number_iterations <= np.iinfo(np.int32).max
    ):
        raise ValueError("max_number_iterations must be a positive int32")
    if number_points < 1 or number_points > np.iinfo(np.int32).max:
        raise ValueError("number_points must be a positive int32")
    if not np.all(np.isfinite(data)) or not np.all(np.isfinite(initial_parameters)):
        raise ValueError("data and initial_parameters must be finite; zero ROI padding")
    if (
        parameters_to_fit.shape != (number_parameters,)
        or not parameters_to_fit.flags.c_contiguous
        or np.any((parameters_to_fit != 0) & (parameters_to_fit != 1))
    ):
        raise ValueError("parameters_to_fit must be a contiguous binary mask")
    if weights is not None and (np.any(~np.isfinite(weights)) or np.any(weights < 0)):
        raise ValueError("weights must be finite and nonnegative")
    if user_info is not None and (not user_info.flags.c_contiguous or user_info.dtype.hasobject):
        raise ValueError("user_info must be contiguous without Python objects")
    if model_id >= 13:
        get_capabilities()
        extras = 0 if model_id == 14 else (2 if model_id == 13 else 3)
        expected = 3 * number_fits * number_points + number_fits + extras
        if (
            user_info is None
            or user_info.dtype != np.float32
            or user_info.shape != (expected,)
            or not np.all(np.isfinite(user_info))
        ):
            raise ValueError("Invalid microscopy user_info layout/dtype; zero coordinate padding")
        sizes = user_info[
            3 * number_fits * number_points : 3 * number_fits * number_points + number_fits
        ]
        if np.any((sizes < 1) | (sizes > number_points) | (sizes != np.floor(sizes))):
            raise ValueError("user_info ROI sizes must be positive integers within number_points")
    if constrained:
        get_capabilities()
        constraint_types = (
            np.zeros(number_parameters, np.int32) if constraint_types is None else constraint_types
        )
        if (
            constraint_types.dtype != np.int32
            or constraint_types.shape != (number_parameters,)
            or not constraint_types.flags.c_contiguous
            or np.any((constraint_types < 0) | (constraint_types > 3))
        ):
            raise ValueError("Invalid constraint_types shape, dtype or enumeration")
        if np.any(constraint_types) and constraints is None:
            raise ValueError("Active constraint types require constraints")
        if constraints is not None:
            if (
                constraints.dtype != np.float32
                or constraints.shape != (number_fits, 2 * number_parameters)
                or not constraints.flags.c_contiguous
                or np.any(np.isnan(constraints))
            ):
                raise ValueError("Invalid constraints shape, dtype, layout or NaN")
            lower, upper = constraints[:, 0::2], constraints[:, 1::2]
            has_lower = (constraint_types == 1) | (constraint_types == 3)
            has_upper = (constraint_types == 2) | (constraint_types == 3)
            both = has_lower & has_upper
            if np.any(~np.isfinite(lower[:, has_lower])) or np.any(
                ~np.isfinite(upper[:, has_upper])
            ):
                raise ValueError("Active constraint endpoints must be finite")
            if np.any(lower[:, both] > upper[:, both]) or np.any(
                (lower == upper)[:, both & (parameters_to_fit == 1)]
            ):
                raise ValueError("Reversed bounds or equal bounds for free parameters")
            if np.any(initial_parameters[:, has_lower] < lower[:, has_lower]) or np.any(
                initial_parameters[:, has_upper] > upper[:, has_upper]
            ):
                raise ValueError("Infeasible initial parameters, including fixed values")

    # we don't check type of user_info, but we extract the size in bytes of it
    if user_info is not None:
        user_info_size = user_info.nbytes
    else:
        user_info_size = 0

    # pre-allocate output variables
    parameters = np.zeros((number_fits, number_parameters), dtype=np.float32)
    states = np.zeros(number_fits, dtype=np.int32)
    chi_squares = np.zeros(number_fits, dtype=np.float32)
    number_iterations = np.zeros(number_fits, dtype=np.int32)

    if not number_fits:
        return parameters, states, chi_squares, number_iterations, 0.0

    def pointer(array, kind):
        return None if array is None else array.ctypes.data_as(POINTER(kind))

    arguments = [
        c_size_t(number_fits),
        c_size_t(number_points),
        pointer(data, c_float),
        pointer(weights, c_float),
        c_int(model_id),
        pointer(initial_parameters, c_float),
    ]
    if constrained:
        arguments.extend((pointer(constraints, c_float), pointer(constraint_types, c_int)))
    arguments.extend(
        (
            c_float(tolerance),
            c_int(max_number_iterations),
            pointer(parameters_to_fit, c_int),
            c_int(estimator_id),
            c_size_t(user_info_size),
            pointer(user_info, c_char),
            pointer(parameters, c_float),
            pointer(states, c_int),
            pointer(chi_squares, c_float),
            pointer(number_iterations, c_int),
        )
    )
    t0 = time.perf_counter()
    status = (constrained_func if constrained else gpufit_func)(*arguments)
    t1 = time.perf_counter()

    # check status
    if status != Status.Ok:
        # get error from last error and raise runtime error
        error_message = error_func()
        raise RuntimeError(f"status = {status}, message = {error_message}")

    # return output values
    return parameters, states, chi_squares, number_iterations, t1 - t0


def get_last_error():
    """
    :return: Error message of last error.
    """
    return error_func()


def cuda_available():
    """
    :return: True if CUDA is available, False otherwise
    """
    return cuda_available_func() != 0


def get_cuda_version():
    """
    :return: Tuple with runtime and driver version as integers.
    """
    runtime_version = c_int(-1)
    driver_version = c_int(-1)
    status = get_cuda_version_func(byref(runtime_version), byref(driver_version))

    # check status
    if status != Status.Ok:
        # get error from last error and raise runtime error
        error_message = error_func()
        raise RuntimeError(f"status = {status}, message = {error_message}")

    # decode versions
    runtime_version = runtime_version.value
    runtime_version = (runtime_version // 1000, runtime_version % 1000 // 10)
    driver_version = driver_version.value
    driver_version = (driver_version // 1000, driver_version % 1000 // 10)

    return runtime_version, driver_version
