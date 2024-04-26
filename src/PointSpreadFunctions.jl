module PointSpreadFunctions
using FourierTools: center_pos, FourierJoin
using FourierTools, NDTools, IndexFunArrays, SpecialFunctions, FFTW
using ZernikePolynomials
export PSFParams, sinc_r, jinc_r_2d, pupil_xyz, apsf, psf, k0, kxy, aplanatic_factor
export get_Abbe_limit, get_Nyquist_limit, get_pupil_aperture
export kz_mid_pos
export modify_STED_exp, modify_STED_hyper, modify_ident, modify_square

export ModeWidefield, ModeConfocal, Mode4Pi, ModeISM, Mode2Photon, ModeSTED, ModeLightsheet
export MethodRichardsWolf, MethodPropagate, MethodPropagateIterative, MethodShell, MethodSincR

include("aplanatic.jl")
include("PSF_types.jl")
include("util.jl")
include("pupil_pol.jl")
include("pupil.jl")
include("aPSFs.jl")

"""
    psf(::Type{ModeWidefield}, sz::NTuple, pp::PSFParams; sampling=get_sampling(sz, pp))

Calculates the widefield single-frequency point spread function (psf), i.e. the image of a single (very small) emitter. 
Most of the parameters (such as refractive index, numerical aperture, vacuum wavelength, aberrations etc.) are hidden in the 
parameter structure argument `pp`,
which should be generated via the `PSFParams()` constructor. See ``PSFParams()` for details.

# Parameters
+ `sz`:         size tuple of the final PSF
+ `pp`:         PSF parameters of the PSF. See `PSFParams()` for details. The argument `pp.aplanatic` defined whether an excitation or emission PSF is calculated.
+ `sampling=nothing`:   The sampling parameters of the resulting PSF. If `nothing` is provided, the PSF will be sampled according to the Abbe limit.
+ `use_resampling=true`: Exploits a calculation trick, which first calculates the amplitudes are calculated on a twice coarser grid and the result is upsampled. This increases the speed but may be less accurate.
+ `return_amp=false`:    Has to be `false` since confocal amplitude spread functions do not exist for non-zero pinhole sizes. 

See also:
+ `apsf`:  calculates the underlying amplitude point spread function (apsf)

Example:
```julia-repl
julia> pp = PSFParams(0.5,1.4,1.52); sz=(128,128,128); sampling=(0.050,0.050,0.200);
# an emission PSF of an isotropic (freely rotating emitter)

julia> p_wf = psf(sz, pp; sampling=sampling);
# an emission PSF of an emitter oriented along the Z direction

julia> pp_dipole = PSFParams(0.5,1.4,1.52; transition_dipole=[0.0,0.0,1.0]);

julia> p_dipole = psf(sz, pp_dipole; sampling=sampling);
```
"""
function psf(::Type{ModeWidefield}, sz::NTuple, pp::PSFParams; sampling=nothing, use_resampling=true, return_amp=false) # unclear why the resampling seems to be so bad
    amp = let
        if use_resampling
            fct_ex = (sz,my_sampling) -> apsf(sz, pp, sampling=my_sampling, center_kz=true)
            calc_with_resampling(fct_ex, sz, sampling, norm_amp=true) # , shift_kz=kz_mid_pos(sz, pp, sampling)
        else
            apsf(sz, pp, sampling=sampling)
        end
    end

    res = amp_to_int(amp, pp)
    if length(sz)<3
        res = dropdims(res, dims=3)
    end
    if return_amp
        return res, amp
    else
        return res
    end
end

"""
    psf(::Type{ModeConfocal}, sz::NTuple, pp_em::PSFParams; pp_ex=nothing, pinhole=nothing, pinhole_ft=disc_pinhole_ft, sampling=nothing, use_resampling=true, return_amp=false, pinhole_positions=[(0.0,0.0)], ex_modifier=modify_ident) # unclear why the resampling seems to be so bad

Calculates a confocal point spread function. The normalisation is such that a completely open `pinhole` diameter yields the excitation PSF with its normalization. 
Returns the PSF or a vector of PSFs.
    
# Parameters
+ `sz`:         size tuple of the final PSF
+ `pp_em`:      PSF parameters of the emission PSF. PSF parameters of the PSF. See `PSFParams()` for details. This should include the emission wavelength
+ `pp_ex=nothing`:      This is a required named parameter, containing all the settings for the excitation PSF. This should include the excitation wavelength as well as typically `aplanatic=aplanatic_illumination`.
+ `pinhole=nothing`:    The diameter of the pinhole in Airy Units (AU = 1.22 λ/NA). A pinhole size of one AU corresponds to a pinhole border falling onto the first zero of a corresponding paraxial emission PSF.
+ `sampling=nothing`:   The sampling parameters of the resulting PSF.
+ `use_resampling=true`: Exploits a calculation trick, which first calculates the individual PSFs on a twice coarser grid in XY and Z and then upsamples the result. Note that this may be inappropriate due to undersampling due to the Stokes shift which is neglected here. But warnings will then result during the calculations of the subsampled widefield PSFs.
+ `return_amp=false`:    Has to be `false` since confocal amplitude spread functions do not exist for non-zero pinhole sizes. 
+ `pinhole_positions=[(0.0,0.0)]`:  A list of pinhole positions in pixels. One PSF will be returned for each pinhole position. If only a single pinhole position is supplied the PSF will directly be returned instead of a vector of PSFs. Specifies the precise position(s) of the pinholes in the detection path. This allows to simulate an offset (misadjusted) pinhole, or (as a vector of tuples) a PSF for a whole set of positions. See the `psf(ModeISM, ...)` for more details. 
+ `pinhole_ft=disc_pinhole_ft`:  Specifies which function is used to calculate the Fourier transform of the pinhole. This allows the user to control the pinhole shape. 
+ `ex_modifier`:   A function that can modify the excitation PSF. By default the identity is used, but other options are `modify_square` to calculate a two-photon confocal PSF. However, you should typically use the `Mode2Photon` to do this. Note that the order of the PSFParams are reversed.
    By default (modify_ident), no modication is performed.

```julia-repl
julia> pp_em = PSFParams(0.5,1.4,1.52; mode=ModeConfocal);
julia> pp_ex = PSFParams(pp_em; λ=0.488, aplanatic=aplanatic_illumination);
julia> p_conf = psf((128,128,128),pp_em; pp_ex=pp_ex, pinhole=0.1, sampling=(0.040,0.040,0.100));
```
"""
function psf(::Type{ModeConfocal}, sz::NTuple, pp_em::PSFParams; sampling=nothing, pp_ex=nothing, use_resampling=true, pinhole=nothing, pinhole_ft=disc_pinhole_ft,  return_amp=nothing, pinhole_positions=[(0.0,0.0)], ex_modifier=modify_ident) # unclear why the resampling seems to be so bad
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
    pp_both = PSFParams(pp_em; λ=λeff)
    # the factor of two below, is since the amp psf can be twice undersampled, but the intensity psf not.
    check_amp_sampling(sz, pp_both, sampling .* 2.0)

    psf_ex = let
        if use_resampling
            fct_ex = (sz,my_sampling) -> psf(ModeWidefield, sz, pp_ex; sampling=my_sampling, use_resampling=use_resampling)
            calc_with_resampling(fct_ex, sz, sampling, norm_amp=false)
        else
            psf(ModeWidefield, sz, pp_ex; sampling=sampling, use_resampling=use_resampling)
        end
    end

    # apply the two-photon effect for excitation if needed.
    psf_ex =  ex_modifier(psf_ex) 

    # pp_em = PSFParams(pp, mode = ModeWidefield)
    psf_em = let
        if use_resampling
            fct_em = (sz,my_sampling) -> psf(ModeWidefield, sz,pp_em; sampling=my_sampling, use_resampling=use_resampling)
            calc_with_resampling(fct_em, sz, sampling, norm_amp=false)
        else
            psf(ModeWidefield, sz, pp_em; sampling=sampling, use_resampling=use_resampling)
        end
    end

    pinhole_pix = pinhole_AU_to_pix(sz, pp_em, sampling, pinhole)
    res = confocal_int(psf_ex, psf_em, pp_em; pinhole_pix=pinhole_pix, pinhole_ft=pinhole_ft, pinhole_positions=pinhole_positions, ex_modifier=ex_modifier)
    if length(sz)<3
        res = dropdims(res, dims=3)
    end
    return res
end

"""
    psf(::Type{ModeLightsheet}, sz::NTuple, pp_em::PSFParams; pp_ex=nothing, sampling=nothing, use_resampling=true, dynamically_scanned=false) # unclear why the resampling seems to be so bad

Calculates a point spread function for a lightsheet microscopy. Not that the coordinate system refers to the detection PSF. The NA of the 
excitation PSF refers to the Z and Y direction. The polarisation of the excitation PSF along X means along Z and along Y means along Y in the final volume.
Note that the excitation is assumed for the best focus position (beam waist). In reality the PSF will be spatially dependent due to the excitation focus.

# Parameters
+ `sz`:         size tuple of the final PSF
+ `pp_em`:      PSF parameters of the emission PSF. PSF parameters of the PSF. See `PSFParams()` for details. This should include the emission wavelength
+ `pp_ex=nothing`:      This is a required named parameter, containing all the settings for the excitation PSF. This should include the excitation wavelength as well as typically `aplanatic=aplanatic_illumination`.
+ `pinhole=nothing`:    The diameter of the pinhole in Airy Units (AU = 1.22 λ/NA). A pinhole size of one AU corresponds to a pinhole border falling onto the first zero of a corresponding paraxial emission PSF.
+ `sampling=nothing`:   The sampling parameters of the resulting PSF.
+ `use_resampling=true`: Exploits a calculation trick, which first calculates the individual PSFs on a twice coarser grid in XY and Z and then upsamples the result. Note that this may be inappropriate due to undersampling due to the Stokes shift which is neglected here. But warnings will then result during the calculations of the subsampled widefield PSFs.
+ 'dynamically_scanned': specifies, whether the lightsheet uses a dynamically scanned beam. This is important for the calculation of the PSF. If the beam is scanned, the PSF will be calculated in the same way as for a confocal PSF. If the beam is not scanned, the excitation PSF will be assuming a focussing via a cylinder lens.
+ `ex_modifier`:   A function that can modify the excitation PSF. By default the identity is used, but other options are `modify_square` to calculate a two-photon confocal PSF. However, you should typically use the `Mode2Photon` to do this. Note that the order of the PSFParams are reversed.
    By default (modify_ident), no modication is performed.
+ sigma_z = nothing: If the user provides a sigma_z value, the excitation calculation is using a Guassian excitation PSF with the specified sigma. pp_ex is ignored.

```julia-repl
julia> pp_em = PSFParams(0.5,0.9,1.33; mode=ModeLightsheet);
julia> pp_ex = PSFParams(0.488, 0.2, 1.33);
julia> p_lightsheet = psf((128,128,128), pp_em; pp_ex=pp_ex, sampling=(0.040,0.040,0.100), dynamically_scanned=true);
julia> p_lightsheet_gauss = psf((128,128,128), pp_em; sampling=(0.040,0.040,0.100), sigma_z=1.0); # Gaussian excitation PSF
```
"""
function psf(::Type{ModeLightsheet}, sz::NTuple, pp_em::PSFParams; sampling=nothing, pp_ex=nothing, use_resampling=true, dynamically_scanned=false, ex_modifier=modify_ident, sigma_z = nothing) # unclear why the resampling seems to be so bad
    if isnothing(pp_ex) && isnothing(sigma_z) 
        error("The named parameter `pp_ex` is obligatory for lightsheet calculation. Provide the excitation PSF parameters here.")
    end

    # check the sampling for the widefield emission PSF
    check_amp_sampling_xy(sz, pp_em, sampling)

    # switch coordinates and sampling for the excitation PSF
    sz_ex = (sz[3], sz[2], sz[3])
    sampling_ex = (sampling[3], sampling[2], sampling[3])

    if isnothing(sigma_z)
        # check th z-sampling.   WRONG!!!:
        max_samp_ex = get_amp_sampling_xy(sz_ex, pp_ex, sampling_ex)
        if (ex_modifier != modify_ident)
            if (ex_modifier == modify_square)
                max_samp_ex = max_samp_ex./2;
            else
                warn("unsupported excitation modifyer for the sampling calculation")
            end
        end
        max_samp_em = get_amp_sampling_z(sz, pp_em, sampling)
        max_samp_total = 1/(1/max_samp_ex[1] + 1/max_samp_em);
        if (sampling[3] > max_samp_total)
            warn("The z-sampling of the lightsheet PSF is undersampled. The z-sampling should be at least $max_samp_total. The current z-sampling is $(sampling[3]).")
        end
    end

    psf_ex = let
        if isnothing(sigma_z)
            psf_ex = let        
                if !(dynamically_scanned)
                    # force it to calculate only the central slice. In this way a coherent cylinder illumination is calculated
                    sz_ex = (sz_ex[1], 1, sz_ex[3])
                end

                if use_resampling
                    fct_ex = (sz,my_sampling) -> psf(ModeWidefield, sz, pp_ex; sampling=my_sampling, use_resampling=use_resampling)
                    calc_with_resampling(fct_ex, sz_ex, sampling_ex, norm_amp=false)
                else
                    psf(ModeWidefield, sz_ex, pp_ex; sampling=sampling, use_resampling=use_resampling)
                end
            end
            psf_ex = permutedims(psf_ex, (3,2,1))
            if (dynamically_scanned)
                # perform the incoherent averaging of intensity.
                psf_ex = sum(psf_ex, dims=2) # because it averages over the Y-direction during scanning
            end
        else
            psf_ex =  gaussian((1,1,sz[3]), sigma=(1,1, sigma_z./sampling[3]));
        end
        psf_ex
    end

    # apply the two-photon effect for excitation if needed.
    psf_ex =  ex_modifier(psf_ex) 

    # pp_em = PSFParams(pp, mode = ModeWidefield)
    psf_em = let
        if use_resampling
            fct_em = (sz,my_sampling) -> psf(ModeWidefield, sz,pp_em; sampling=my_sampling, use_resampling=use_resampling)
            calc_with_resampling(fct_em, sz, sampling, norm_amp=false)
        else
            psf(ModeWidefield, sz, pp_em; sampling=sampling, use_resampling=use_resampling)
        end
    end

    res = psf_ex .* psf_em
    if length(sz)<3
        res = dropdims(res, dims=3)
    end
    return res
end

"""
    modify_ident(h)

A function from a series of modficator functions which serve to modify the excitation intensity.
    This particular version does nothing. It just returns the input `h`.
"""
function modify_ident(h)
    h
end

"""
    modify_square(h)

A function from a series of modficator functions which serve to modify the excitation intensity.
    This version returns the absolute square of the input `h`. It is useful for two-photon microscopy. 
"""
function modify_square(h)
    abs2.(h)
end

"""
    modify_STED_exp(h, h_doughnut)

A function from a series of modficator functions which serve to modify the excitation intensity.
This version applies a second STED PSF `h_doughnut` to deplete the excitation PSF `h`.
The depletion is modeled by an exponential decay as obtained by a relatively short STED pulse following a very short excitation pulse.
Note the `h_doughnut` should be the depletion PSF. Note that this function should be used as a tool to define a new function with fixed second argument.
"""
function modify_STED_exp(h, h_doughnut)
    h .* exp.(.- h_doughnut)
end

"""
    modify_STED_hyper(h, h_doughnut)

A function from a series of modficator functions which serve to modify the excitation intensity.
    This version applies a second STED psf `h_doughnut` to deplete the excitation psf `h`.
    The depletion is modeled by an hyperbolic saturation curve based on anaylsing the steady-state kinetics.
    This is mostly suitable for CW-STED.
    Note the `h_doughnut` should be the depletion PSF. Note that this function should be used as a tool to define a new function with fixed second argument.
    The normalizations are such that `h` as well as `h_doughnut` reache 50% saturation at `h==1.0`. This means that typically your want to use
    `max(h) << 1` and `max(h_doughnut) >> 1`
"""
function modify_STED_hyper(h, h_doughnut)
    h ./  (one(eltype(h)) .+ h .+ h_doughnut)
end

"""
    pinhole_AU_to_pix(sz, pp_em, sampling)

returns the pinhole in AU. Also checks if the pinhole is too big.
"""
function pinhole_AU_to_pix(sz, pp_em, sampling, pinhole)
     # now we need to modify the sampling such that the pinhole corrsponds to the equivalent of one Airy Unit.
    # The Airy Unit is the diameter of the Airy disc: 1.22 * lamda_em / NA 
    # AU = 1.22 * pp_em.λ / pp_em.NA
    if isnothing(pinhole)
        return nothing
    end
    AU_pix = AU_per_pixel(pp_em, sampling) # AU ./ sampling[1:2]
    pinhole_in_AU = pinhole .* AU_pix
    if any(pinhole_in_AU .> sz[1:2]) 
        @warn "Pinhole is larger than image this leads to serious aliasing artefacts. Maximal pinhole size is: $(sz[1:2]./AU_pix) AU. Assuming fully open pinhole."
        return nothing
    end
    return pp_em.dtype.(pinhole_in_AU)
end

"""
    confocal_int(psf_ex, psf_em, pinhole=nothing, pinhole_ft=disc_pinhole_ft, sampling=nothing, use_resampling=true, return_amp=nothing, pinhole_positions=[(0.0,0.0)], ex_modifier=modify_ident)

helper function to calculate the confocal detection, once the excitation and emission PSFs and pinhole parameters are defined. Two photon excitation and ISM detection are also included.
"""
function confocal_int(psf_ex, psf_em, pp_em;  pinhole_pix=nothing, pinhole_ft=disc_pinhole_ft, pinhole_positions=[(0.0,0.0)], ex_modifier=modify_ident)
   
    # This can be done a lot more efficiently by staying in Fourier space. Ideally even by only calculating half the range of the jinc function:
    # pinhole = real.(ift2d(jinc_r_2d(sz, pinhole .* AU_pix, pp_em.dtype)))
    # pinhole_ft = rfft2d(ifftshift(pinhole))

    all_PSFs = Vector{typeof(psf_em)}();
    sz = size(psf_em)

    # apply two-photon or STED effect 
    psf_ex = ex_modifier(psf_ex)

    rfft_psf_em = rfft2d(psf_em)
    for p in pinhole_positions
        if isnothing(pinhole_pix)
            # fully open pinhole assumed. The integral over all emission PSFs should really be one!
            push!(all_PSFs, psf_ex)
        else
            my_pinhole_ft = exp_ikx_rfft(eltype(psf_em), sz, shift_by=p) .* pinhole_ft(sz, pp_em, pinhole_pix)
            # pinhole_ft = rfft2d(ifftshift(pinhole))
            # return pinhole_ft
            my_em =  irfft2d(rfft_psf_em .* my_pinhole_ft, sz[1])
            push!(all_PSFs, my_em .* psf_ex) # 
            # push!(all_PSFs, irfft2d(my_pinhole_ft, sz[1]))  # only for diagnostic purposes
        end
    end

    if length(all_PSFs) == 1
        return all_PSFs[1]
    else
        return all_PSFs
    end
end

"""
    disc_pinhole_ft(sz, pp, pinhole)

calculates the Fourier transform of a disc-shaped pinhole

+ `pp`: PSF parameter structure
+ `sampling`: sampling of the image pixels in object coordinates
+ `pinhole`: diameter of the pinhole in pixels
"""
function disc_pinhole_ft(sz, pp, pinhole)
    jinc_r_2d(sz, pinhole, pp.dtype; r_func= PointSpreadFunctions.rr_rfft)
end

"""
    disc_pinhole_ft(sz, pp, pinhole)

calculates the Fourier transform of a box-shaped pinhole

+ `pp`: PSF parameter structure
+ `sampling`: sampling of the image pixels in object coordinates
+ `pinhole`: tuple of side lengths of the pinhole in pixels
"""
function box_pinhole_ft(sz, pp, pinhole)
    sinc_r_2d(sz, pinhole, pp.dtype)
end

"""
    AU_in_pixels(pp, sampling)

calculates the size of the Airy unites in pixels.
#arguments
+ `pp`: PSF parameter structure
+ `sampling`: sampling of the image pixels in object coordinates
"""
function AU_per_pixel(pp, sampling)
    AU = 1.22 * pp.λ / pp.NA
    AU_pix = AU ./ sampling[1:2]
    return AU_pix
end

"""
    exp_ikx_rfft(dtype, sz; shift_by)

calculates the exp factor for an rfft to result into shifting the corresponding real array by `shift_by`.
"""
function exp_ikx_rfft(dtype, sz; shift_by)
    szrfft = rfft_size(sz)
    # this is needed to account for the exp_ikx function being designed for FFTs rather than rfts
    shift_by = shift_by .* eltype(shift_by).(szrfft[1:2] ./ sz[1:2])
    ifftshift(exp_ikx(Complex{dtype}, szrfft[1:2], offset=CtrRFT, shift_by=shift_by), (2:length(sz)))
end

function ism_positions_rect(pinhole_dist, pinhole_grid)
    [pinhole_dist.*((n,m).-((pinhole_grid.+1)./2 )) for n=1:pinhole_grid[1] for m=1:pinhole_grid[2]], pinhole_dist
end

"""
    psf(::Type{ModeISM}, sz::NTuple, pp_em::PSFParams; pinhole=nothing, sampling=nothing, pinhole_ft=box_pinhole_ft, pinhole_positions=[(0.0,0.0)], pinhole_dist=0.5, pinhole_grid=(5,5), ism_pos=ism_positions_rect, args...)

Calculates a confocal point spread function. The normalisation is such that a completely open `pinhole` diameter yields the excitation PSF with its normalization. 
Returns the PSF or a vector of PSFs.
    
# Parameters
+ `sz`:         size tuple of the final PSF
+ `pp_em`:      PSF parameters of the emission PSF. This should include the emission wavelength
+ `pp_ex=nothing`:      This is a required named parameter, containing all the settings for the excitation PSF. This should include the excitation wavelength as well as typically `aplanatic=aplanatic_illumination`.
+ `pinhole=nothing`:   The diameter of each pinhole in Airy Units (AU = 1.22 λ/NA). By default (pinholes=nothing) the pinhole size is automatically calculated by the `pinhole_dist` parameter below to yield mutually touching pinholes.
+ `sampling=nothing`:   The sampling parameters of the resulting PSF.
+ `use_resampling=true`: Exploits a calculation trick, which first calculates the individual PSFs on a twice coarser grid in XY and Z and then upsamples the result. Note that this may be inappropriate due toundersampling due to the Stokes shift which is neglected here. But warnings will then result during the calculations of the subsampled widefield PSFs.
+ `return_amp=false`:    Has to be `false` since confocal amplitude spread functions do not exist for non-zero pinhole sizes. 
+ `pinhole_positions=nothing`:  A list of pinhole positions in pixel coordinates. One PSF will be returned for each pinhole position. If only a single pinhole position is supplied the PSF will directly be returned instead of a vector of PSFs. Be careful, since this is not the same as the position of the relative shift of the images.
+ `pinhole_dist=0.5`: a value or tuple specified the distances between pinholes, when arranged in a grid. 
+ `pinhole_ft=box_pinhole_ft`:  Specifies which function is used to calculate the Fourier transform of the pinhole. This allows the user to control the pinhole shape. E.g. hexagonal pattern with round pinholes

```julia-repl
julia> sz=(128,128,128); sampling = (0.04,0.04,0.120)
julia> pp_em = PSFParams(0.5,1.4,1.52; mode=ModeISM);
julia> pp_ex = PSFParams(pp_em; λ=0.488, aplanatic=aplanatic_illumination);
julia> p_ism = psf(sz,pp_em; pp_ex=pp_ex, pinhole=0.21, pinhole_dist=0.2, sampling=sampling);
```
"""
function psf(::Type{ModeISM}, sz::NTuple, pp_em::PSFParams; pinhole=nothing, sampling=nothing, pinhole_ft=box_pinhole_ft, pinhole_positions=nothing, pinhole_dist=0.12, pinhole_grid=(5,5), ism_pos=ism_positions_rect, args...)
    if isnothing(pinhole_positions) || isempty(pinhole_positions)
        #create a list of pinhole positions
        pinhole_dist = pinhole_dist .* AU_per_pixel(pp_em, sampling)
        pinhole_positions, ph_diam = ism_pos(pinhole_dist, pinhole_grid)
        # convert to AUs
        ph_diam = ph_diam ./ AU_per_pixel(pp_em, sampling)
        if !isnothing(pinhole)
            if any(pinhole .> ph_diam)
                @warn "The given pinhole size $(pinhole) AU is larger that the mutual pinhole spacing $(ph_diam). This is unphysical."
            end
        else
            pinhole = ph_diam
        end
    end
    psf(ModeConfocal, sz, pp_em; pinhole=pinhole, sampling=sampling, pinhole_ft=pinhole_ft, pinhole_positions=pinhole_positions, args...) # confocal is able to handle this well
end

"""
    psf(::Type{Mode2Photon}, sz::NTuple, pp_ex::PSFParams; pp_em=nothing, pinhole=nothing, sampling=nothing, pinhole_ft=disc_pinhole_ft, args...)

Calculates a 2-photon (potentially confocal) point spread function.  
Returns the PSF or a vector of PSFs.
    
# Parameters
+ `sz`:         size tuple of the final PSF
+ `pp_ex`:      PSF parameters of the excitation PSF. This should include the exission wavelength (typically in the IR region). Please make sure to also set `pp.aplanatic=aplanatic_illumination`.
+ `pp_em`:      PSF parameters of the emission PSF. This only needs to be supplied, if a pinhole is used. Otherwise non-descanned detection (NDD) is assumed.
+ `pinhole=nothing`:   If `nothing`, NDD is assumed and the PSF is only the square of the excitation PSF. The diameter of each pinhole in Airy Units (AU = 1.22 λ/NA). 
+ `sampling=nothing`:   The sampling parameters of the resulting PSF.
+ `use_resampling=true`: Exploits a calculation trick, which first calculates the individual PSFs on a twice coarser grid in XY and Z and then upsamples the result. Note that this may be inappropriate due toundersampling due to the Stokes shift which is neglected here. But warnings will then result during the calculations of the subsampled widefield PSFs.
+ `return_amp=false`:    Has to be `false` since confocal amplitude spread functions do not exist for non-zero pinhole sizes. 
+ `pinhole_positions=[(0.0,0.0)]`:  A list of pinhole positions in pixels. One PSF will be returned for each pinhole position. If only a single pinhole position is supplied the PSF will directly be returned instead of a vector of PSFs.

```julia-repl
julia> sz=(128,128,128); sampling = (0.04,0.04,0.120)
julia> pp_ex = PSFParams(0.8,1.4,1.52; mode=Mode2Photon, aplanatic= aplanatic_illumination);
julia> p_2p = psf(sz,pp_ex; sampling=sampling);
```
"""
function psf(::Type{Mode2Photon}, sz::NTuple, pp::PSFParams; sampling=nothing, pinhole=nothing, pinhole_ft=disc_pinhole_ft, pinhole_positions=[(0.0,0.0)], args...)
    if isnothing(pinhole) 
        return abs2.(psf(ModeWidefield, sz, pp; sampling=sampling))
    else
        psf(ModeConfocal, sz, pp; ex_2p=true, pinhole=pinhole, sampling=sampling, pinhole_ft=pinhole_ft,pinhole_positions=pinhole_positions, args...) # confocal is able to handle this well
    end
end

"""
    psf(::Type{ModeSTED}, sz::NTuple, pp_ex::PSFParams, pp_sted::PSFParams; pp_em=nothing, pinhole=nothing, sampling=nothing, pinhole_ft=disc_pinhole_ft, args...)


Calculates a 2-photon (potentially confocal) point spread function.  
Returns the PSF or a vector of PSFs.
    
# Parameters
+ `sz`:         size tuple of the final PSF
+ `pp_ex`:      PSF parameters of the excitation PSF. This should include the exission wavelength (typically in the IR region). Please make sure to also set `pp.aplanatic=aplanatic_illumination`.
+ `pp_sted`:      PSF parameters of the excitation PSF. This should include the exission wavelength (typically in the IR region). Please make sure to also set `pp.aplanatic=aplanatic_illumination`.
+ `pp_em`:      PSF parameters of the emission PSF. This only needs to be supplied, if a pinhole is used. Otherwise non-descanned detection (NDD) is assumed.
+ `pinhole=nothing`:   If `nothing`, NDD is assumed and the PSF is only the square of the excitation PSF. The diameter of each pinhole in Airy Units (AU = 1.22 λ/NA). 
+ `sampling=nothing`:   The sampling parameters of the resulting PSF.
+ `use_resampling=true`: Exploits a calculation trick, which first calculates the individual PSFs on a twice coarser grid in XY and Z and then upsamples the result. Note that this may be inappropriate due toundersampling due to the Stokes shift which is neglected here. But warnings will then result during the calculations of the subsampled widefield PSFs.
+ `return_amp=false`:    Has to be `false` since confocal amplitude spread functions do not exist for non-zero pinhole sizes. 
+ `pinhole_positions=[(0.0,0.0)]`:  A list of pinhole positions in pixels. One PSF will be returned for each pinhole position. If only a single pinhole position is supplied the PSF will directly be returned instead of a vector of PSFs.

```julia-repl
julia> sz=(128,128,128); sampling = (0.04,0.04,0.120)
julia> pp_em = PSFParams(0.510,1.4,1.52; mode=ModeSTED); # STED suppression beam
julia> pp_ex = PSFParams(0.488,1.4,1.52; aplanatic= aplanatic_illumination);
julia> pp_sted = PSFParams(0.570,1.4,1.52; method=MethodPropagateIterative, aplanatic= aplanatic_illumination, pol=pol_circ_spiral);
julia> p_sted_01 = psf(sz, pp_em; pp_ex=pp_ex, pinhole= 2.0, sampling=sampling, pp_sted=pp_sted, rel_sted_intensity=0.1);
julia> p_sted_1 = psf(sz, pp_em; pp_ex=pp_ex, pinhole= 2.0, sampling=sampling, pp_sted=pp_sted, rel_sted_intensity=1.0);
julia> p_sted_5 = psf(sz, pp_em; pp_ex=pp_ex, pinhole= 2.0, sampling=sampling, pp_sted=pp_sted, rel_sted_intensity=5.0);
julia> pp_sted_b = PSFParams(0.570,1.4,1.52; method=MethodPropagateIterative, aplanatic= aplanatic_illumination, pol=pol_circ_tophat);
julia> p_sted_b5 = psf(sz, pp_em; pp_ex=pp_ex, pinhole= 2.0, sampling=sampling, pp_sted=pp_sted_b, rel_sted_intensity=5.0);
julia> using View5D; @vt p_sted_01 p_sted_1 p_sted_5 p_sted_b5
```
"""
function psf(::Type{ModeSTED}, sz::NTuple, pp::PSFParams; sampling=nothing, pp_sted::PSFParams, rel_sted_intensity = 1.0, pinhole=nothing, pinhole_ft=disc_pinhole_ft, pinhole_positions=[(0.0,0.0)], modify_STED=modify_STED_exp, ex_modifier=modify_ident, args...)

    h_STED = psf(ModeWidefield, sz, pp_sted; sampling=sampling)
    h_STED = rel_sted_intensity .* h_STED ./ maximum(h_STED)
    # return h_STED

    ex_modifier_new = (x) -> modify_STED(ex_modifier(x), h_STED)

    psf(ModeConfocal, sz, pp; ex_modifier=ex_modifier_new, pinhole=pinhole, sampling=sampling, pinhole_ft=pinhole_ft, pinhole_positions=pinhole_positions, args...) # confocal is able to handle this well
end

"""
    psf(::Type{Mode4Pi}, sz::NTuple, pp_ex::PSFParams; sampling=nothing, pp_ex2=pp_ex, pp_em=nothing, pp_em2=nothing, rel_ex_phase = 0.0, rel_em_phase = 0.0, pinhole=nothing, pinhole_ft=disc_pinhole_ft, pinhole_positions=[(0.0,0.0)], ex_modifier=modify_square)

Calculates a 4Pi point spread function. Note that the default is using a two-photon excitation.
Returns the PSF or a vector of PSFs.
    
# Parameters
+ `sz`:         size tuple of the final PSF
+ `pp_ex`:      PSF parameters of the excitation PSF. This should include the exission wavelength (typically in the IR region). Please make sure to also set `pp.aplanatic=aplanatic_illumination`.
+ `pp_ex2`:     PSF parameters of the other side excitation PSF. If `nothing` is provided, single-sided excitation (e.g. 4Pi Type B) is assumed.
+ `ex_modifier=modify_square`:  If `modify_square`, two-photon excitation is assumed. Use `modify_ident` for single-photon excitation.
+ `pp_em`:      PSF parameters of the emission PSF. This only needs to be supplied, if a pinhole is used.
+ `pp_em2`:     PSF parameters of the other side emission PSF. If `nothing is supplied, Type A 4Pi microscopy is assumed with single-sided (or incoherent) detection.`
+ `pinhole=nothing`:   If `nothing`, NDD is assumed and the PSF is only the square of the excitation PSF. The diameter of each pinhole in Airy Units (AU = 1.22 λ/NA). 
+ `sampling=nothing`:   The sampling parameters of the resulting PSF.
+ `pinhole_positions=[(0.0,0.0)]`:  A list of pinhole positions in pixels. One PSF will be returned for each pinhole position. If only a single pinhole position is supplied the PSF will directly be returned instead of a vector of PSFs.

```julia-repl
julia> sampling = (0.04,0.04,0.04)
julia> sz = (128,128,128)
julia> pp_em = PSFParams(0.5,1.3,1.52; mode=Mode4Pi, pol=pol_x);
julia> pp_ex = PSFParams(0.800,1.3,1.52; aplanatic=aplanatic_illumination, mode=Mode4Pi, pol=pol_x);
julia> @time p_4pi = psf(sz,pp_ex; pp_ex2=pp_ex, pp_em=pp_em, pp_em2=pp_em, sampling=sampling, pinhole=2.0, ex_modifier=modify_square);
```
"""
function psf(::Type{Mode4Pi}, sz::NTuple, pp_ex::PSFParams; sampling=nothing, pp_ex2 = pp_ex, pp_em=nothing, pp_em2=nothing, rel_ex_phase = 0.0, rel_em_phase = 0.0, pinhole=nothing, pinhole_ft=disc_pinhole_ft, pinhole_positions=[(0.0,0.0)], ex_modifier=modify_square)
    if isnothing(pp_ex)
        error("For a 4Pi PSF, you have to at least provide one excitation psf using the named argument `pp_ex`.")
    end
    
    if ex_modifier == modify_square # a special case that warrants the resampling trick twice
        Nq_ex = (pp_ex.λ/2) /pp_ex.NA / 2
        Nq_em = pp_em.λ /pp_em.NA / 2
        Nq_XY = 1/(1/Nq_ex + 1/Nq_em)
        Nq_ex = (pp_ex.λ/2) /pp_ex.n / 2
        Nq_em = pp_em.λ /pp_em.n / 2
        Nq_Z = 1/(1/Nq_ex + 1/Nq_em)
        if isnothing(sampling)
            sampling = (Nq_XY,Nq_XY,Nq_Z)
            @warn "Result sampling was calculated to be $(sampling)."
        else
            if sampling[3] > Nq_Z 
                @warn "The 4Pi PSF needs to be sampled at least at $(Nq_Z), but sampling is $(sampling[3])."
            end
        end
    else
        Nq_ex = pp_ex.λ /pp_ex.NA / 2
        Nq_em = pp_em.λ /pp_em.NA / 2
        Nq_XY = 1/(1/Nq_ex + 1/Nq_em)
        Nq_ex = pp_ex.λ /pp_ex.n / 2
        Nq_em = pp_em.λ /pp_em.n / 2
        Nq_Z = 1/(1/Nq_ex + 1/Nq_em)
        if isnothing(sampling)
            sampling = (Nq_XY,Nq_XY,Nq_Z)
            @warn "Result sampling was calculated to be $(sampling)."
        else
            if sampling[3] > Nq_Z
                @warn "The 4Pi PSF needs to be sampled at least at $(Nq_Z), but sampling is $(sampling[3])."
            end
        end
    end
    asf_ex1 = apsf(sz, pp_ex; sampling=sampling)
    psf_ex = let 
        if isnothing(pp_ex2)
            amp_to_int(asf_ex1, pp_ex)
        else
            asf_ex2 = (pp_ex == pp_ex2) ? asf_ex1 : apsf(sz,pp_ex2; sampling=sampling)
            amp_to_int(asf_ex1 .+ exp.(1im.*pp_ex.dtype(rel_ex_phase)) .* conj.(asf_ex2), pp_ex2)
        end        
    end

    asf_em1 = apsf(sz,pp_em; sampling=sampling)
    psf_em = let
        if isnothing(pp_em2) # 4Pi type A
            amp_to_int(asf_em1, pp_em)
        else # 4Pi type C
            asf_em2 = (pp_em == pp_em2) ? asf_em1 : apsf(sz,pp_em2; sampling=sampling)
            amp_to_int(asf_em1 .+  exp.(1im.*pp_em.dtype(rel_em_phase)) .* conj.(asf_em2), pp_em2)
        end
    end

    pinhole_pix = pinhole_AU_to_pix(sz, pp_em, sampling, pinhole)
    return confocal_int(psf_ex, psf_em, pp_em; pinhole_pix=pinhole_pix, pinhole_ft=pinhole_ft, pinhole_positions=pinhole_positions, ex_modifier=ex_modifier)
end


"""
    psf(sz::NTuple, pp::PSFParams; sampling=get_sampling(sz, pp))

Calculates the point spread function (psf), i.e. the image of a single (very small) emitter. Most of the parameters
(such as refractive index, numerical aperture, vacuum wavelength, aberrations etc.) are hidden in the parameter structure argument `pp`,
which should be generated via the `PSFParams()` constructor. See ``PSFParams()` for details.
Note that the field `pp.mode` defines the microscopic mode to simulate. Currently implemented are the default `ModeWidefield` and `ModeConfocal`.

See also:
+ apsf():  calculates the underlying amplitude point spread function (apsf)

Example:
```julia-repl
julia> pp = PSFParams(0.5,1.4,1.52);
julia> p = psf((128,128,128),pp; sampling=(0.050,0.050,0.200));
```
"""
function psf(sz::NTuple, pp::PSFParams; nps...) # unclear why the resampling seems to be so bad
    return psf(pp.mode, sz, pp; nps...)
end

end # module
