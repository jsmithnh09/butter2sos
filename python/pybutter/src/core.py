import numpy as np
from os import path
from ctypes import *
import atexit
from typing import Tuple


def _initialize_dll():
    """initialize the DLL and set the return types and argtypes."""
    curdir = path.abspath(path.dirname(__file__))
    sosdll = cdll.LoadLibrary(path.join(curdir, "butterlib"))
    sosdll.butter.argtypes = [c_int, c_double, c_double, c_int]
    sosdll.butterband.argtypes = [c_int, c_double, c_double, c_double, c_int]
    sosdll.butter.restype = POINTER(c_double)
    sosdll.butterband.restype = POINTER(c_double)
    return sosdll


# DLL closer-callback on Python instance closure.
def _close_dll(cdll: LibraryLoader):
    del cdll


# instantiate the module-scope variable.
_SOSDLL = _initialize_dll()
_NCOEFFS = 6

# register and pass the library handle.
atexit.register(_close_dll, _SOSDLL)


def stable(matrix: np.ndarray) -> Tuple[bool, int]:
    """indicate if the SOS matrix is stable.

    Instead of constructing an eigensolve companion matrix,
    just use the quadratic formula, (unless np.roots is faster(?))

    Parameters
    ----------
    matrix: np.ndarray
        The Second Order Sections, where each row is a biquad.

    Returns
    -------
    out: Tuple[bool, int]
        A tuple with a flag indicating if stable and an optional
        integer indicating the stage number if the biquad was unstable.
    """
    for stgInd in range(matrix.shape[0]):
        a, b, c = matrix[stgInd, 3:]
        sq = np.sqrt(b**2 - 4 * a * c)
        denom = 2 * a
        roots = [(-b + sq) / denom, (-b - sq) / denom]
        if any(abs(roots)) >= 1:  # beyond unit circle?
            return (False, stgInd)
    return (True, None)


def _numstages(ord, type="lowpass") -> int:
    """compute the number of stages based on the filter type."""
    if type in ["lowpass", "highpass", "allpass"]:
        ncount = (ord - 1) / 2 + 1 if (ord % 2) else ord / 2  # 6 coeffs/stage.
    elif type in ["bandpass", "bandstop"]:
        ncount = ord
    else:
        raise ValueError(f"filter type {type} is not valid.")

    return int(ncount)


def butter(order, fc, fs=44100.0, type="lowpass") -> np.ndarray:
    """Butterworth filter design, (single corner frequency.)

    The C-level API into butterlib.dll. Supports several
    filter types.

    Parameters
    -----------
    order: int
        The order of the filter.
    fc: float
        The corner or "natural" frequency point.
    fs: float
        The sampling rate of the filter.
    type: str
        The type of filter. Supported filters are lowpass,
        highpass, and allpass.

    Returns
    -------
    sosmat: np.ndarray
        The SOS matrix generated the C-level API.

    Raises:
    ------
    In the event the C-level API produces any NaN's, (saw this
    when performing the bilinear xform on floats and not doubles),
    we do a check on the matrix to ensure the poles are stable.

    """
    if fc >= fs / 2:
        raise ValueError("Corner frequency %g cannot exceed Nyquist." % fc)
    elif fs <= 0:
        raise ValueError("Sampling rate %g cannot go below DC." % fs)
    elif type not in ["lowpass", "highpass", "allpass"]:
        raise ValueError(f"Unknown filter type {type}.")

    N = _numstages(order, type)

    # give the integer filter type
    if type[0] == "l":
        e_type = 0
    elif type[0] == "h":
        e_type = 1
    else:
        e_type = 2

    N = _numstages(order, type)
    mat = _SOSDLL.butter(c_int(order), c_double(fc), c_double(fs), c_int(e_type))
    sosmat = np.reshape(
        np.asmatrix([mat[i] for i in range(N * _NCOEFFS)]), (N, _NCOEFFS)
    )
    (status, stg) = stable(sosmat)
    if not status:
        raise ValueError(
            f"SOS matrix from butterlib produces unstable poles at stage ({stg})."
        )
    del mat
    return sosmat


def butterband(order, flow, fhigh, fs=44100.0, type="bandpass") -> np.npdarray:
    """Band-based Butterworth filter design.

    Two corner-frequency based filter design (bandpass/bandstop).

    Parameters
    ----------
    order: int
        The order of the filter. Equal to the number of biquads produced.
    flow: float
        The lower corner frequency of the filter.
    fhigh: float
        The higher corner frequency of the filter.
    fs: float
        The discrete sampling rate of the filter.
    type: str
        The filter-type string, indicating "bandpass" or "bandstop".

    Returns
    -------
    sosmatrix: np.ndarray
        The quantized SOS matrix from the butterlib DLL.

    Raises
    ------
    If the filter has unstable poles due to an error in the C-code,
    we'll raise a value error.
    """
    if type not in ["bandpass", "bandstop"]:
        raise ValueError(f"Unknown filter type {type}")
    elif flow >= fhigh:
        flow, fhigh = fhigh, flow

    if fhigh >= fs / 2:
        raise ValueError("Upper corner frequency exceeds Nyquist.")

    # bandpass = 0, bandstop = 1
    if type[-1] == "s":
        e_type = 0
    else:
        e_type = 1

    # interface and parse the C-function.
    N = _numstages(order, type)
    mat = _SOSDLL.butterband(
        c_int(order), c_double(flow), c_double(fhigh), c_double(fs), c_int(e_type)
    )
    sosmat = np.reshape(
        np.asmatrix([mat[i] for i in range(N * _NCOEFFS)]), (N, _NCOEFFS)
    )
    (status, stg) = stable(sosmat)
    if not status:
        raise ValueError(
            f"SOS matrix from butterlib produces unstable poles at stage ({stg})."
        )
    del mat
    return sosmat
