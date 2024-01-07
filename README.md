# butter2sos
[![CI](https://github.com/jsmithnh09/butter2sos/actions/workflows/CI.yml/badge.svg)](https://github.com/jsmithnh09/butter2sos/actions/workflows/CI.yml)

Z-domain Butterworth Second Order Section (SOS) filter design in C/MATLAB/Julia. Filter design is based on the recursive filter design shown in \[1-3\]. The motivation for this is to learn more about C programming with an interesting DSP project. The supported filter types are lowpass, highpass, and allpass, and soon to include bandpass/bandstop types.

Based on the target filter type, the zero positions are already known, and the biquad coefficients are then easier to evaluate. The gain change from pole/zero movement is based on what's provided in the DSP.jl Julia package \[4\].

## API

For external language interfacing, the `jl_butter` and `jl_butterband` provide the argument inputs where the external language provides a buffer, and the underlying functions will populate the coefficients, (SOS matrix is treated as row-major via C-conventions, column-major orientation would need post-processing.)

To use this within Python, you can run `pip install .` and compare the results against `scipy.signal.butter`. The zero/pole ordering is different, but python tests show that the frequency response is within precision.

## Bulding
To build, a `makefile` is present for generating the different binaries and libraries that are present. In an effort to use more "modern" build tools, [bazel](https://bazel.build/) files are present as well.

## Postmortem Observations
Windows/Microsoft does not [document their complex numbers API well.](https://learn.microsoft.com/en-us/cpp/c-runtime-library/complex-math-support?view=msvc-170) Details of the underlying `_Dcomplex` type definition are not present, unless you dive into code directly. The underlying type uses a `_Val` array member where the first and second indices are respectively the real and imaginary components. This varies from the GNU `complex` type, hence the scattered `_WIN32` directives.

The same lack of documentation can be shared with bazel, the reference on their site is scarse. To get the DLL to properly expose the symbols for Python to interface with via the `ctypes` module, I had to use the `win_def_file` argument. The `.lds` file did not expose the symbols properly unfortunately, only the `.def` format would work on my Windows machine.

## References
\[1\] J.G. Proakis and D.G. Manolakis, Digital Signal Processing, Prentice Hall, 2007, chapter 10, section 3.

\[2\] J.O. Smith III, Digital Filters with Audio Applications, BookSurge Publishing, 2007, appendix I.

\[3\] A.V. Oppenheim and R.W. Schafer, Digital Signal Processing, Prentice Hall, 1975, chapter 5, sections 1 through 3.

\[4\] [JuliaDSP (2020) Filters/design.jl (Version 0.7.3)](https://github.com/JuliaDSP/DSP.jl/blob/master/src/Filters/Filters.jl)
