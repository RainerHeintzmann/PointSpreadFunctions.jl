module PSFs
using FourierTools: center_pos, FourierJoin
using FourierTools, NDTools, IndexFunArrays, SpecialFunctions, FFTW
using ZernikePolynomials
export PSFParams, sinc_r, jinc_r_2d, pupil_xyz, apsf, psf, k0, kxy, aplanatic_factor
export ModeWidefield, ModeConfocal, Mode4Pi

include("aplanatic.jl")
include("PSF_types.jl")
include("util.jl")
include("pupil_pol.jl")
include("pupil.jl")
include("aPSFs.jl")

"""
    psf(sz::NTuple, pp::PSFParams; sampling=get_sampling(sz, pp))

calculates the point spread function (psf), i.e. the image of a single (very small) emitter. Most of the parameters
(such as refractive index, numerical aperture, vacuum wavelength, aberrations etc.) are hidden in the parameter structure argument `pp`,
which should be generated via the `PSFParams()` constructor. See ``PSFParams()` for details.

See also:
+ apsf():  calculates the underlying aplitude point spread function (apsf)

Example:
```jdoctest
julia> pp = PSFParams(500.0,1.4,1.52);
julia> p = psf((128,128,128),pp; sampling=(50,50,100)); #; 
```
"""
function psf(sz::NTuple, pp::PSFParams; sampling=nothing, use_resampling=true, return_amp=false) # unclear why the resampling seems to be so bad
    if use_resampling == false
        amp = apsf(sz, pp, sampling=sampling)
        if return_amp
            return amp_to_int(amp), amp
        else
            return amp_to_int(amp)
        end
    end
    extra_layers = 2
    small_sz=ceil.(Int,sz./2) .+ extra_layers
    big_sz = small_sz .* 2  # size after upsampling
    if ~isnothing(sampling)
        amp_sampling = sampling .* big_sz ./ small_sz
    else
        amp_sampling = nothing
    end
    amp = apsf(small_sz, pp, sampling=amp_sampling, center_kz=true)
    if pp.FFTPlan != nothing
        P1d = plan_fft(amp,(1,2), flags=pp.FFTPlan)
        P1id = plan_ifft(amp,(1,2), flags=pp.FFTPlan)
    end
    border_in = (0,0,ceil.(Int,sz[3]./2) ./ small_sz[3],0)
    mywin = collect(window_hanning((1,1,small_sz[3],1), border_in=border_in, border_out=1)) # for speed reasons the collect is faster
    amp = upsample2(amp .* mywin,dims=(1,2,3))
    res = sum(abs2.(amp),dims=4)[:,:,:,1]
    if true # any(isodd.(sz))
        if return_amp
            return select_region!(res,new_size=sz), amp
        else
            return select_region!(res,new_size=sz)
        end
    else
        if return_amp
            return res, amp
        else
            return res
        end
    end
end

end # module
