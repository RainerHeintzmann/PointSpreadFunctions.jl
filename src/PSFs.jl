module PSFs
using FourierTools: center_pos
using FourierTools, NDTools, IndexFunArrays, SpecialFunctions
export PSFParams, sinc_r, jinc_r_2d, pupil_xyz, apsf, psf, k0, kxy, aplanatic_factor
export ModeWidefield, ModeConfocal, Mode4Pi

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
function psf(sz::NTuple, pp::PSFParams; sampling=nothing)
    amp, sampling = apsf(sz, pp, sampling=sampling)
    if isnothing(sampling)
        sum(upsample2_abs2(amp),dims=4)
    else
        amp_to_int(amp)
    end
end

end # module
