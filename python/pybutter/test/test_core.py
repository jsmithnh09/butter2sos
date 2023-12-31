import pybutter as pb
from scipy.signal import butter, sosfreqz
import numpy as np
import pytest

_ERR_TOL = np.sqrt(np.finfo(np.float32).eps)
_ORDERS = range(1, 12, 1)
_FS = 48e3
_NPTS = 1024
_FREQV = np.logspace(np.log10(30), np.log10((_FS / 2) * 0.9), num=10)
_TYPES = {"low": "lowpass", "high": "highpass"}
_BANDS = np.logspace(np.log10(30), np.log10((_FS / 3) * 0.9), num=3)

def test_singleband():
    for tkey in _TYPES.keys():
        for order in _ORDERS:
            for freq in _FREQV:
                # generate matrices and compare resonses.
                pysos = butter(order, freq, btype=tkey, fs=_FS, output="sos")
                csos = pb.butter(order, freq, fs=_FS, type=_TYPES[tkey])
                (pyh, _) = sosfreqz(pysos, worN=_NPTS, fs=_FS)
                (ch, _) = sosfreqz(csos, worN=_NPTS, fs=_FS)

                # determine if the response is close...
                assert np.allclose(pyh[1:], ch[1:], atol=_ERR_TOL)

def test_multiband():
    for bandtype in ["bandpass", "bandstop"]:
        for order in _ORDERS:
            for freq in _BANDS:
                delta = freq * 0.2
                Bl, Bh = freq - delta, freq + delta
                pysos = butter(order, (Bl, Bh), fs=_FS, btype=bandtype, output="sos")
                csos = pb.butterband(order, Bl, Bh, fs=_FS, type=bandtype)
                (pyh, _) = sosfreqz(pysos, worN=_NPTS, fs=_FS)
                (ch, _) = sosfreqz(csos, worN=_NPTS, fs=_FS)

                # determine the bandpass/bandstop is relatively close as well...
                assert np.allclose(pyh[1:], ch[1:], atol=_ERR_TOL)

def test_sosbin(tmp_path):
    """Use a temporary file to generate the SOSBIN and check the matrix response is OK."""
    fh = tmp_path / "tmp"
    if not fh.parent.is_dir():
        fh.parent.mkdir()
    tpath = str(fh)
    dummy1 = butter(7, 12345, fs=_FS, btype='highpass', output='sos')
    
    # read and write the "dummy" sos matrix.
    pb.writesosbin(tpath, dummy1)
    dummy2 = pb.readsosbin(tpath)

    (H1, _) = sosfreqz(dummy1, worN=_NPTS, fs=_FS)
    (H2, _) = sosfreqz(dummy2, worN=_NPTS, fs=_FS)

    assert np.allclose(H1[1:], H2[1:], atol=_ERR_TOL)
