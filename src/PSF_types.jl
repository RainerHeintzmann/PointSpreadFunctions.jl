abstract type PSFMode end
struct ModeWidefield <: PSFMode end
struct ModeConfocal <: PSFMode end
struct Mode4Pi <: PSFMode end
struct STED <: PSFMode end

abstract type PSFMethod end
struct MethodSincR <: PSFMethod end
struct MethodPropagate <: PSFMethod end
struct MethodRichardsWolf <: PSFMethod end
struct MethodShell <: PSFMethod end

"""
    PSFParams
    This structure stores all the general parameters needed to calculate the point spread function.
    Only pixel-pitch and image size information is handles seperately.
Members are:
λ:            Vacuum wavelength
NA:           numerical aperture
n.            refractive index of the embedding AND immersion medium
dtype:        real-valued data type to generate PSF for
mode:         microscopy mode to calculate PSF for   ::PSFMode
polarization: a function calculating the polarization from a given Tuple of relative-k pupil vectors
aplanatic:    aplanatic factor. Provided as a function of angle θ
method:       the method of calculation
"""
struct PSFParams
    λ  # Vacuum wavelength
    NA # numerical aperture
    n # refractive index of the embedding AND immersion medium
    dtype # real-valued data type to generate PSF for
    mode # microscopy mode (widefield, confocal, etc.) to calculate PSF for ::PSFMode
    polarization # a function calculating the polarization from a given Tuple of relative-k pupil vectors. Can also be pol_scalar()
    aplanatic # aplanatic factor. Provided as a function of angle θ
    method # the method of calculation (e.g. PSFs.MethodPropagate).
end

function PSFParams(my_λ=500, my_NA=1.2, my_n=1.33; pol=pol_scalar, dtype=Float32, mode=ModeWidefield, aplanatic = (θ) -> sqrt.(max.(0,cos.(θ))), method=MethodPropagate)
    PSFParams(my_λ, my_NA, my_n, dtype, mode, pol, aplanatic, method)
end

