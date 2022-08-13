from src.core import butter as cbutter
from scipy.signal import butter, sosfreqz
import numpy as np

_ERR_TOL = np.sqrt(np.finfo(np.float32).eps)
_ORDERS = range(1, 12, 1)
_FS = 48e3
_NPTS = 1024
_FREQV = np.logspace(np.log10(30), np.log10((_FS / 2) * 0.9), num=10)
_TYPES = {"low": "lowpass", "high": "highpass"}


for tkey in _TYPES.keys():
    for order in _ORDERS:
        for freq in _FREQV:
            # generate matrices and compare resonses.
            pysos = butter(order, freq / (_FS / 2), btype=tkey, fs=_FS, output="sos")
            csos = cbutter(order, freq, fs=_FS, type=_TYPES[tkey])
            (pyh, _) = sosfreqz(pysos, worN=_NPTS, fs=_FS)
            (ch, _) = sosfreqz(csos, worN=_NPTS, fs=_FS)

            # determine if the response is close...
            assert np.allclose(pyh[1:], ch[1:], atol=_ERR_TOL)
