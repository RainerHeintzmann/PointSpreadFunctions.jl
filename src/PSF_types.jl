abstract type PSFMode end
struct ModeWidefield <: PSFMode end
struct ModeConfocal <: PSFMode end
struct Mode4Pi <: PSFMode end
struct STED <: PSFMode end

abstract type PSFMethod end
struct MethodSincR <: PSFMethod end
struct MethodPropagate <: PSFMethod end
struct MethodPropagateIterative<: PSFMethod end
struct MethodRichardsWolf <: PSFMethod end
struct MethodShell <: PSFMethod end
export Aberrations
export Zernike_Piston, Zernike_Tilt, Zernike_Tip, Zernike_ObliqueAstigmatism, Zernike_Defocus, Zernike_VerticalAstigmatism, Zernike_VerticalTrefoil
export Zernike_VerticalComa, Zernike_HorizontalComa, Zernike_ObliqueTrefoil, Zernike_ObliqueQuadrafoil, Zernike_ObliqueSecondaryAstigmatism, Zernike_Spherical
export Zernike_VerticalSecondaryAstigmatism, Zernike_VerticalQuadrafoil

const Zernike_Piston = 0
const Zernike_Tilt = 1
const Zernike_Tip = 2
const Zernike_ObliqueAstigmatism = 3
const Zernike_Defocus = 4
const Zernike_VerticalAstigmatism = 5
const Zernike_VerticalTrefoil = 6
const Zernike_VerticalComa = 7
const Zernike_HorizontalComa= 8
const Zernike_ObliqueTrefoil = 9
const Zernike_ObliqueQuadrafoil = 10
const Zernike_ObliqueSecondaryAstigmatism = 11
const Zernike_Spherical = 12
const Zernike_VerticalSecondaryAstigmatism = 13
const Zernike_VerticalQuadrafoil = 14

"""
    Aberrations
    defining Zernike phase aberrations via a list of indices and coefficients and an indexing style.
"""
struct Aberrations
    indices # zernike indices as a vector
    coefficients # corresponding zernike coefficients
    index_style
end

"""
    Aberrations(indices=[],coefficients=[];index_style = :OSA)

defining Zernike phase aberrations via a list of indices and coefficients and an indexing style.
By default, no Aberrations are defined.

Arguments:
+ indices:  Vector of indices
+ coefficients: Vector of corresponding coefficients
+ index_style: type of indexing used. By defaul :OSA is used (See: https://en.wikipedia.org/wiki/Zernike_polynomials#OSA/ANSI_standard_indices)

Here is a list of constants defining the main indices (OSA style):
`:Zernike_Piston` = 0
`:Zernike_Tilt` = 1
`:Zernike_Tip` = 2
`:Zernike_ObliqueAstigmatism` = 3
`:Zernike_Defocus` = 4
`:Zernike_VerticalAstigmatism` = 5
`:Zernike_VerticalTrefoil` = 6
`:Zernike_VerticalComa` = 7
`:Zernike_HorizontalComa` = 8
`:Zernike_ObliqueTrefoil` = 0
`:Zernike_ObliqueQuadrafoil` = 10
`:Zernike_ObliqueSecondaryAstigmatism` = 11
`:Zernike_Spherical` = 12
`:Zernike_VerticalSecondaryAstigmatism` = 13
`:Zernike_VerticalQuadrafoil` = 14

See also: 
+ package `ZernikePolynomials.jl`, and https://en.wikipedia.org/wiki/Zernike_polynomials
"""
function Aberrations(indices=[],coefficients=[];index_style = :OSA)
    Aberrations(indices,coefficients,index_style)
end

"""
    PSFParams

This structure stores all the general parameters needed to calculate the point spread function.
Only pixel-pitch and image size information is handles seperately.

Structure Members:
+ λ:            Vacuum wavelength
+ NA:           numerical aperture
+ n:            refractive index of the embedding AND immersion medium
+ dtype:        real-valued data type to generate PSF for
+ mode:         microscopy mode to calculate PSF for. See the constructor PSFParams() for more details.
+ polarization: a function calculating the polarization from a given Tuple of relative-k pupil vectors
+ aplanatic:    aplanatic factor. Provided as a function of angle θ
+ method:       the method of calculation
+ FFTPlan:      information on how to calculate the FFTW plan

See also:
+ PSFParams()

Example:
```jdoctest
julia> using FFTW, PSFs

julia> ppm = PSFParams(580.0, 1.4, 1.518;pol=pol_scalar,method=PSFs.MethodSincR, aberrations= aberr, FFTPlan=FFTW.MEASURE)
PSFParams(580.0, 1.4, 1.518, Float32, ModeWidefield, PSFs.pol_scalar, PSFs.var"#42#43"(), PSFs.MethodSincR, 0x00000000, nothing, nothing)

```
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
    FFTPlan # if not nothing this will be the plan of how to perform FFTs
    aberrations 
    pixelshape # here functions can be supplied, which account for the pixelshape influence (itegration over pixel size).
end

"""
    PSFParams(λ=500, NA=1.2, n=1.33; pol=pol_scalar, dtype=Float32, mode=ModeWidefield, 
    aplanatic = aplanatic_detection, method=MethodPropagateIterative, FFTPlan=nothing,
    aberrations=Aberrations(), pixelshape=nothing)

creates the PSFParams structure via this constructor.

Arguments:
+ λ:            Vacuum wavelength
+ NA:           numerical aperture
+ n:            refractive index of the embedding AND immersion medium
+ pol: a function calculating the polarization from a given Tuple of relative-k pupil vectors
+ dtype:        real-valued data type to generate PSF for
+ mode:         microscopy mode to calculate PSF for ::PSFMode. Currently only `ModeWidefield` is implemented
+ method:         microscopy mode to calculate PSF for 
                valid options are currently:
                + MethodPropagate: Angulare spectrum propagation. This version does NOT account for wrap around problems yielding problems at larger out-of-focus distances
                + MethodPropagateIterativ (default): Angulare spectrum propagation accounting from wrap-around problems in each propagation step by applying a perfectly matched layer (PML).
                + MethodShell: Angulare spectrum propagation with a slightly different calculation order. This version does NOT account for wrap around problems yielding problems at larger out-of-focus distances
                + MethodSincR: Based on first calculating a SincR function in real space and applying consecutive filtering steps. It accounts for wrap around problems but requires a quite high sampling along the Z direction.
+ aplanatic:    aplanatic factor. Provided as a function of angle θ. 
+ FFTPlan:      information on how to calculate the FFTW plan. Default: nothing (using FFTW.ESTIMATE)

Example:
```jdoctest
julia> using FFTW, PSFs

julia> ppm = PSFParams(580.0, 1.4, 1.518;pol=pol_scalar,method=PSFs.MethodSincR, aberrations= aberr, FFTPlan=FFTW.MEASURE)
PSFParams(580.0, 1.4, 1.518, Float32, ModeWidefield, PSFs.pol_scalar, PSFs.var"#42#43"(), PSFs.MethodSincR, 0x00000000, nothing, nothing)

```
"""
function PSFParams(my_λ=500, my_NA=1.2, my_n=1.33; pol=pol_scalar, dtype=Float32, mode=ModeWidefield, 
                    aplanatic = aplanatic_detection, method=MethodPropagateIterative, FFTPlan=nothing,
                    aberrations=Aberrations(), pixelshape=nothing)
    PSFParams(my_λ, my_NA, my_n, dtype, mode, pol, aplanatic, method, FFTPlan, aberrations, pixelshape)
end

