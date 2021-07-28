module PSFs
using FourierTools: center_pos, FourierJoin
using FourierTools, NDTools, IndexFunArrays, SpecialFunctions
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
Example:
pp = PSFParams(500.0,1.4,1.52)
p = psf((128,128,128),pp; sampling=(50,50,100)); #; 
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
    small_sz=ceil.(Int,sz./2)
    big_sz = small_sz .* 2
    if ~isnothing(sampling)
        amp_sampling = sampling .* big_sz ./ small_sz
    else
        amp_sampling = nothing
    end
    amp = apsf(small_sz, pp, sampling=amp_sampling, center_kz=true)
    amp = upsample2(amp,dims=(1,2,3))
    res = sum(abs2.(amp),dims=4)[:,:,:,1]
    if any(isodd.(sz))
        if return_amp
            return NDTools.select_region(res,new_size=sz), amp
        else
            return NDTools.select_region(res,new_size=sz)
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
