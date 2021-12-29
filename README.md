# butter2sos
Z-domain Butterworth Second Order Section (SOS) filter design in C/MATLAB/Julia.
Filter design is based on the recursive filter design shown in \[1-3\]. The motivation 
for this is to learn more about C programming with an interesting DSP project.
The supported filter types are lowpass, highpass, and allpass, and soon to include 
bandpass/bandstop types.

To build:
```
gcc -std=c99 butter2sos_design.c butter2sos.c -o butter2sos
```
To build a shared library for interfacing with Julia, (via `jl_butter` and `jl_butterband`):
```
gcc -std=c99 butter2sos_design.c -shared -o butterlib.dll
```

Based on the target filter type, the zero positions are already known,
and the biquad coefficients are then easier to evaluate.

The gain change from pole/zero movement is based on what's provided in
the DSP.jl Julia package \[4\].

Spectral responses between the matrix in this project and the solution
provided by MATLAB is below -100 dB, which is pretty good given 
floating point machine precision is approximately -138 dB.

References:
\[1\] J.G. Proakis and D.G. Manolakis, Digital Signal Processing, Prentice
Hall, 2007, chapter 10, section 3.

\[2\] J.O. Smith III, Digital Filters with Audio Applications, BookSurge
Publishing, 2007, appendix I.

\[3\] A.V. Oppenheim and R.W. Schafer, Digital Signal Processing, Prentice
Hall, 1975, chapter 5, sections 1 through 3.

\[4\] JuliaDSP (2020) Filters/design.jl (Version 0.7.3)
https://github.com/JuliaDSP/DSP.jl/blob/master/src/Filters/Filters.jl
