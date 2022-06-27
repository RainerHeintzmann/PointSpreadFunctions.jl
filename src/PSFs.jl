module PSFs
using FourierTools: center_pos, FourierJoin
using FourierTools, NDTools, IndexFunArrays, SpecialFunctions, FFTW
using ZernikePolynomials
export PSFParams, sinc_r, jinc_r_2d, pupil_xyz, apsf, psf, k0, kxy, aplanatic_factor
export ModeWidefield, ModeConfocal, Mode4Pi, get_Abbe_limit, get_Nyquist_limit

include("aplanatic.jl")
include("PSF_types.jl")
include("util.jl")
include("pupil_pol.jl")
include("pupil.jl")
include("aPSFs.jl")

"""
    psf(::Type{ModeWidefield}, sz::NTuple, pp::PSFParams; sampling=get_sampling(sz, pp))

calculates the widefield single-frequency point spread function (psf), i.e. the image of a single (very small) emitter. Most of the parameters
(such as refractive index, numerical aperture, vacuum wavelength, aberrations etc.) are hidden in the parameter structure argument `pp`,
which should be generated via the `PSFParams()` constructor. See ``PSFParams()` for details.

See also:
+ apsf():  calculates the underlying amplitude point spread function (apsf)

Example:
```jdoctest
julia> pp = PSFParams(0.5,1.4,1.52);
julia> p = psf((128,128,128),pp; sampling=(0.50,0.50,0.200));
```
"""
function psf(::Type{ModeWidefield}, sz::NTuple, pp::PSFParams; sampling=nothing, use_resampling=true, return_amp=false) # unclear why the resampling seems to be so bad
    sz, sampling, is2d = let 
        if length(sz)>2
            sz, sampling, false
        else
            (sz[1:2]...,1), (sampling[1:2]..., eps(eltype(sampling))), true
        end
    end

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
    small_sz, big_sz = let 
        if sz[3]==1
            (small_sz[1:2]...,1), (big_sz[1:2]...,1)
        else
            small_sz, big_sz
        end        
    end
    if ~isnothing(sampling)
        amp_sampling = sampling .* big_sz ./ small_sz
    else
        amp_sampling = nothing
    end
    amp = apsf(small_sz, pp, sampling=amp_sampling, center_kz=true)
    if !isnothing(pp.FFTPlan)
        P1d = plan_fft(amp,(1,2), flags=pp.FFTPlan)
        P1id = plan_ifft(amp,(1,2), flags=pp.FFTPlan)
    end
    border_in = (0,0, ceil.(Int,sz[3]./2) ./ small_sz[3],0)
    mywin = collect(window_hanning((1,1,small_sz[3],1), border_in=border_in, border_out=1)) # for speed reasons the collect is faster

    # Note that the upsampling leads to a one-pixel shift of the center for each odd-size dimension
    # This is taken care of in the select_region code below
    amp = let
        if sz[3] == 1
            upsample2(amp .* mywin,dims=(1,2))
        else
            upsample2(amp .* mywin,dims=(1,2,3))
        end
    end
    res = sum(abs2.(amp),dims=4)[:,:,:,1]
    res = let
        if is2d
            dropdims(res, dims=3)
        else
            res
        end
    end
    if true # any(isodd.(sz))
        recenter = size(res).÷2 .- 2 .*(size(res).÷4) 
        if return_amp
            return select_region_view(res, center = 1 .+ size(res).÷2 .- recenter, new_size=sz), amp
        else
            return select_region_view(res, center = 1 .+ size(res).÷2 .- recenter, new_size=sz)
        end
    else
        if return_amp
            return res, amp
        else
            return res
        end
    end
end

"""
    psf(::Type{ModeConfocal}, sz::NTuple, pp_em::PSFParams; pp_ex=nothing, pinhole=nothing, sampling=nothing, use_resampling=true, return_amp=false) # unclear why the resampling seems to be so bad
    Calculates a confocal point spread function. The normalisation is such that a completely open `pinhole` diameter yields the excitation PSF with its normalization. 
    
#Parameters
+ `sz`:         size tuple of the final PSF
+ `pp_em`:      PSF parameters of the emission PSF. This should include the emission wavelength
+ `pp_xe=nothing`:      This is a required named parameter, containing all the settings for the excitation PSF. This should include the excitation wavelength as well as typically `aplanatic=aplanatic_illumination`.
+ `pinhole=nothing`:    The diameter of the pinhole in Airy Units (AU = 1.22 λ/NA). A pinhole size of one AU corresponds to a pinhole border falling onto the first zero of a corresponding paraxial emission PSF.
+ `sampling=nothing`:   The sampling parameters of the resulting PSF.
+ `use_resampling=true`: Exploits a calculation trick, which first calculates the individual PSFs on a twice coarser grid in XY and Z and then upsamples the result. Note that this may be inappropriate due toundersampling due to the Stokes shift which is neglected here. But warnings will then result during the calculations of the subsampled widefield PSFs.
+ `return_amp=false`:    Has to be `false` since confocal amplitude spread functions do not exist for non-zero pinhole sizes. 

```jdoctest
julia> pp_em = PSFParams(0.5,1.4,1.52; mode=ModeConfocal);
julia> pp_ex = PSFParams(pp_em; λ=0.488, aplanatic=aplanatic_illumination);
julia> p = psf((128,128,128),pp_em; pp_ex=pp_ex, pinhole=0.1, sampling=(0.040,0.040,0.100));
```
"""
function psf(::Type{ModeConfocal}, sz::NTuple, pp_em::PSFParams; pp_ex=nothing, pinhole=nothing, sampling=nothing, use_resampling=true, return_amp=nothing) # unclear why the resampling seems to be so bad
    if isnothing(pp_ex) 
        error("The named parameter `pp_ex` is obligatory for confocal calculation. Provide the excitation PSF parameters here.")
    end
    if isnothing(pinhole)
        error("The named parameter `pinhole` is obligatory for confocal calculation. Provide the excitation PSF parameters here.")
    end
    if !isnothing(return_amp) && return_amp == true
        error("A confocal PSF cannot return an amplitude. Please use `return_amp=false`.")
    end

    # creat a pseudo parameter structure with the combined wavelength just to check the individual amplitude samplings of the final result.
    λeff = 1 / (1/pp_ex.λ + 1/pp_em.λ)
    pp_both = PSFParams(pp_em; λ= λeff)
    # the factor of two below, is since the amp psf can be twice undersampled, but the intensity psf not.
    check_amp_sampling(sz, pp_both, sampling .* 2.0)

    psf_ex = let
        if use_resampling
            fct_ex = (sz,my_sampling) -> psf(ModeWidefield, sz, pp_ex; sampling=my_sampling, use_resampling=use_resampling)
            calc_with_resampling(fct_ex, sz, sampling)
        else
            psf(ModeWidefield, sz, pp_ex; sampling=sampling, use_resampling=use_resampling)
        end
    end

    # pp_em = PSFParams(pp, mode = ModeWidefield)
    psf_em = let
        if use_resampling
            fct_em = (sz,my_sampling) -> psf(ModeWidefield, sz,pp_em; sampling=my_sampling, use_resampling=use_resampling)
            calc_with_resampling(fct_em, sz, sampling)
        else
            psf(ModeWidefield, sz, pp_em; sampling=sampling, use_resampling=use_resampling)
        end
    end

    # now we need to modify the sampling such that the pinhole corrsponds to the equivalent of one Airy Unit.
    # The Airy Unit is the diameter of the Airy disc: 1.22 * lamda_em / NA 
    AU = 1.22 * pp_em.λ / pp_em.NA
    AU_pix = AU ./ sampling[1:2]
    # This can be done a lot more efficiently by staying in Fourier space. Ideally even by only calculating half the range of the jinc function:
    # pinhole = real.(ift2d(jinc_r_2d(sz, pinhole .* AU_pix, pp_em.dtype)))
    # pinhole_ft = rfft2d(ifftshift(pinhole))
    pinhole_ft = jinc_r_2d(sz, pinhole .* AU_pix, pp_em.dtype; r_func= PSFs.rr_rfft)
    # pinhole_ft = rfft2d(ifftshift(pinhole))
    # return pinhole, pinhole_ft
    my_em = irfft2d(rfft2d(psf_em) .* pinhole_ft, sz[1])
    return my_em .* psf_ex
end

"""
    psf(sz::NTuple, pp::PSFParams; sampling=get_sampling(sz, pp))

calculates the point spread function (psf), i.e. the image of a single (very small) emitter. Most of the parameters
(such as refractive index, numerical aperture, vacuum wavelength, aberrations etc.) are hidden in the parameter structure argument `pp`,
which should be generated via the `PSFParams()` constructor. See ``PSFParams()` for details.
Note that the field `pp.mode` defines the microscopic mode to simulate. Currently implemented are the default `ModeWidefield` and `ModeConfocal`.

See also:
+ apsf():  calculates the underlying amplitude point spread function (apsf)

Example:
```jdoctest
julia> pp = PSFParams(0.5,1.4,1.52);
julia> p = psf((128,128,128),pp; sampling=(0.050,0.050,0.200));
```
"""
function psf(sz::NTuple, pp::PSFParams; nps...) # unclear why the resampling seems to be so bad
    return psf(pp.mode, sz, pp; nps...)
end

end # module
