import numpy as np
from os import path
from ctypes import *
import atexit


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
def close_dll(cdll: LibraryLoader):
    del cdll


# instantiate the module-scope variable.
_SOSDLL = _initialize_dll()
_NCOEFFS = 6

# register and pass the library handle.
atexit.register(close_dll, _SOSDLL)


def _numstages(ord, type="lowpass") -> int:
    """compute the number of stages based on the filter type."""
    if type in ["lowpass", "highpass", "allpass"]:
        ncount = (ord - 1) / 2 + 1 if (ord % 2) else ord / 2  # 6 coeffs/stage.
    elif type in ["bandpass", "bandstop"]:
        ncount = ord
    else:
        raise ValueError(f"filter type {type} is not valid.")

    return int(ncount)


def butter(order, fc, fs=44100.0, type="lowpass") -> np.matrix:
    """single corner-frequency butterworth filter design."""
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
    del mat
    return sosmat


def butterband(order, flow, fhigh, fs=44100.0, type="bandpass") -> np.matrix:
    """two corner frequency butterworth filter design."""
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
    del mat
    return sosmat
