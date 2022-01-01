using DSP

"""
    sos = butter2sos(order, Fc, Fs, type=:lowpass)

Z-domain IIR Butterworth design of Lowpass/Highpass/Allpass filters.
`Fc` indicates the corner "-3 dB" natural frequency point, whereas `Fs` indicates
the discrete sample rate and `order` specifies the cascade filter order. The returned
`sos` SecondOrderSection type is compatible with DSP.jl.
"""
function butter2sos(order::Integer, Fc::Real, Fs::Real, type::Symbol=:lowpass)
    (Fc > 0 && Fs > 0) || error("Fc and/or Fs must be greater than 0 Hz (DC).")
    (Fc <= Fs/2) || error("Upperband `Fhi` cannot exceed Nyquist.")
    (type ∈ [:lowpass, :highpass, :allpass]) || error("type must indicate ':lowpass', ':highpass', or ':allpass.'")
    Nstages = ceil(Int, order/2)
    matptr = Vector{Float32}(undef, Nstages*6)
    if (type == :lowpass)
        ftype = Int(0)
    elseif (type == :highpass)
        ftype = Int(1)
    else
        ftype = Int(2)
    end
    ccall((:jl_butter, "butterlib"), Cvoid, (Cint, Cfloat, Cfloat, Cint, Ptr{Float32}), order, Fc, Fs, ftype, matptr)
    sosmat = Matrix{Float32}(transpose(reshape(matptr, (6, Nstages))))
    Vb = Vector{Biquad{:z, Float32}}(undef, Nstages)
    K = prod(sosmat[:,1]) # normalize by lead coefficient and extract gain, (whichever ordering.)
    kstg = 1
    for iStg = 1:Nstages
        kstg = sosmat[iStg,1]
        sosmat[iStg,1:3] ./= kstg
        Vb[iStg] = Biquad{:z, Float32}(sosmat[iStg,1], sosmat[iStg,2], sosmat[iStg,3], sosmat[iStg,5], sosmat[iStg,6])
    end
    SecondOrderSections(Vb, K)
end
        
"""
    sos = butterband(order, Fl, Fh, Fs, type=:bandpass)

Z-domain IIR Butterworth design of Bandpass/Bandstop filters.
`Fl` and `Fh` indicate the lower and upper natural frequency edges,
with `Fs` indicating the discrete sample rate, and `order` specifies
the cascaded Biquads filter order.
"""
function butterband(order::Integer, Fl::Real, Fh::Real, Fs::Real, type::Symbol=:bandpass)
    (Fs > 0) || error("Fs sampling rate must be greater than 0 Hz (DC).")
    (order > 0) || error("Filter order must be greater than zero.")
    (type ∈ [:bandpass, :bandstop]) || error("type must indicate ':bandpass' or ':bandstop'.")
    if (Fl > Fh)
        Fl, Fh = Fh, Fl
    end
    ((Fl > 0 && Fl < Fs/2) && (Fh > 0 && Fh < Fs/2)) || error("Upper/lower bounds must be in range (0, Fs/2].")
    matptr = Vector{Float32}(undef, order*6)
    ftype = (type == :bandpass) ? Int(0) : Int(1)
    ccall((:jl_butterband, "butterlib"), Cvoid, (Cint, Cfloat, Cfloat, Cfloat, Cint, Ptr{Float32}), order, Fl, Fh, Fs, ftype, matptr)
    sosmat = Matrix{Float32}(transpose(reshape(matptr, (6, order))))
    Vb = Vector{Biquad{:z, Float32}}(undef, order)
    kstg = 1
    K = prod(sosmat[:,1])
    for iStg = 1:order
        kstg = sosmat[iStg,1]
        sosmat[iStg,1:3] ./= kstg
        Vb[iStg] = Biquad{:z, Float32}(sosmat[iStg,1], sosmat[iStg,2], sosmat[iStg,3], sosmat[iStg,5], sosmat[iStg,6])
    end
    SecondOrderSections(Vb,K)
end
    



    



