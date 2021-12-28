"""
    sos = butter2sos(order, Fc, Fs, type="lowpass")

Z-domain IIR Butterworth design of Lowpass/Highpass/Allpass filters.
`Fc` indicates the corner "-3 dB" natural frequency point, whereas `Fs` indicates
the discrete sample rate and `order` specifies the cascade filter order. The returned
`sos` matrix contains the raw second-order section coefficients that can be used to
construct a `SecondOrderSection` type.
"""
function butter2sos(order::Integer, Fc::Float32, Fs::Float32; type="lowpass")
    (Fc > 0 && Fs > 0) || error("Fc and/or Fs must be greater than 0 Hz (DC).")
    (Fc <= Fs/2) || error("Upperband `Fhi` cannot exceed Nyquist.")
    (type ∈ ["lowpass", "highpass", "allpass"]) || error("type must indicate 'lowpass', 'highpass', or 'allpass.'")

    Nstages = ceil(Int, order/2)
    ncoeffs = Nstages * 6
    matptr = Ptr{Float32}(undef)
    
    # based on input, discern filter type.
    if (type == "lowpass")
        ftype = Int(0)
    elseif (type == "highpass")
        ftype = Int(1)
    else
        ftype = Int(2)
    end

    # butterlib.dll has both butter and butterband available.
    matptr = ccall((:butter, "butterlib"), Ptr{Float32}, (Cint, Cfloat, Cfloat, Cint), order, Fc, Fs, ftype)
    
    sosmatrix[1:ncoeffs] = matptr[1:ncoeffs]
    reshape(sosmatrix, (Nstages, ncoeffs))
    sosmatrix
end
        


    



