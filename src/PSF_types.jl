abstract type PSFMode end
struct ModeWidefield <: PSFMode end
struct ModeConfocal <: PSFMode end
struct ModeLightsheet <: PSFMode end
struct ModeISM <: PSFMode end
struct Mode2Photon <: PSFMode end
struct Mode4Pi <: PSFMode end
struct ModeSTED <: PSFMode end

abstract type PSFMethod end
struct MethodSincR <: PSFMethod end
struct MethodPropagate <: PSFMethod end
struct MethodPropagateIterative<: PSFMethod end
struct MethodRichardsWolf <: PSFMethod end
struct MethodShell <: PSFMethod end
struct MethodParaxial <: PSFMethod end

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
```jdoctest
# spherical aberration with circular detection polarisation
julia> aberr_sp = Aberrations([Zernike_Spherical],[0.1]);
julia> pp_sp = PSFParams(0.488, 1.4, 1.52; method=MethodPropagateIterative, aberrations= aberr_sp);
julia> p_sp = psf(sz, pp_sp; sampling=sampling);

# diagonal astigmatism wiht x-polarised detection
julia> aberr_as = Aberrations([Zernike_VerticalAstigmatism, Zernike_ObliqueAstigmatism],[0.1, 0.1]);
julia> pp_as = PSFParams(0.488, 1.4, 1.52; method=MethodPropagateIterative, aberrations= aberr_as, pol=pol_x);
julia> p_as = psf(sz, pp_as; sampling=sampling);
```
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
+ polarization: a function calculating the polarization from a given Tuple of relative-k pupil vectors. The constructur argument is `pol`.
+ aplanatic:    aplanatic factor. Provided as a function of angle θ
+ method:       the method of calculation
+ FFTPlan:      information on how to calculate the FFTW plan

See also:
+ PSFParams()

Example:
```jdoctest
julia> using FFTW, PointSpreadFunctions

julia> ppm = PSFParams(0.580, 1.4, 1.518;pol=pol_scalar,method=PointSpreadFunctions.MethodRichardsWolf, aberrations= aberr, FFTPlan=FFTW.MEASURE)
PSFParams(0.58, 1.4, 1.518, Float32, ModeWidefield, PointSpreadFunctions.pol_scalar, PointSpreadFunctions.var"#42#43"(), PointSpreadFunctions.MethodSincR, 0x00000000, nothing, nothing)

```
"""
struct PSFParams
    λ  # Vacuum emission wavelength in µm
    NA # numerical aperture
    n # refractive index of the embedding AND immersion medium
    dtype # real-valued data type to generate PSF for
    mode # microscopy mode (widefield, confocal, etc.) to calculate PSF for ::PSFMode
    polarization # a function calculating the polarization from a given Tuple of relative-k pupil vectors. Can also be pol_scalar()
    aplanatic # aplanatic factor. Provided as a function of angle θ
    method # the method of calculation (e.g. PointSpreadFunctions.MethodPropagate).
    FFTPlan # if not nothing this will be the plan of how to perform FFTs
    aberrations 
    pixelshape # here functions can be supplied, which account for the pixelshape influence (itegration over pixel size).
    transition_dipole # a value of `nothing` means that the total power is summed over X, Y and Z. Otherwise a 3D tuple of values needs to be supplied.
    λ_weights # `nothing` means that only a single wavelength is calculated. Otherwise a tuple of lambdas and wavelengths needs to be supplied.
end

"""
    PSFParams(λ=0.5, NA=1.2, n=1.33; pol=pol_circ, dtype=Float32, mode=ModeWidefield, 
    aplanatic = aplanatic_detection, method=MethodRichardsWolf, FFTPlan=nothing,
    aberrations=Aberrations(), pixelshape=nothing, transition_dipole=nothing, λ_weights=nothing)

creates the PSFParams structure via this constructor. You can also call `PSFParams` with the first argument being an existing structure and just specify the parameters to change via one or multiple named arguments.

Arguments:
+ λ:            Vacuum emsision wavelength (same units as sampling, default is 0.5 µm)
+ NA:           numerical aperture
+ n:            refractive index of the embedding AND immersion medium
+ pol:          a function calculating the polarization from a given Tuple of relative-k pupil vectors
+ dtype:        real-valued data type to generate PSF for
+ mode:         microscopy mode to calculate PSF for ::PSFMode. 
    + ModeWidefield (default): Widefield microscopy
    + Mode2Photon: Two-photon microscopy
    + Mode4Pi: 4Pi microscopy
    + ModeConfocal: Confocal microscopy
    + ModeISM: Image scanning microscopy
    + ModeSTED: stimulated emission depletion microscopy
    + ModeLSCylinder: Light-sheet microscopy illuminating with a cylindrical lens
    + ModeDSLM: Dynamically scanned light-sheet microscopy
+ method:         microscopy mode to calculate PSF for 
                valid options are currently:
    + MethodPropagate: Angulare spectrum propagation. This version does NOT account for wrap around problems yielding problems at larger out-of-focus distances
    + MethodPropagateIterativ (default): Angulare spectrum propagation accounting from wrap-around problems in each propagation step by applying a perfectly matched layer (PML).
    + MethodShell: Angulare spectrum propagation with a slightly different calculation order. This version does NOT account for wrap around problems yielding problems at larger out-of-focus distances
    + MethodSincR: Based on first calculating a SincR function in real space and applying consecutive filtering steps. It accounts for wrap around problems but requires a quite high sampling along the Z direction.
    + MethodRichardsWolf: Uses the method described in the paper by B. Richards and E. Wolf, "Electromagnetic diffraction in optical systems. II. structure of the image field in an aplanatic system," Proc. R. Soc. London A, vol. 253, no. 1274, 1959.
                            The terms I0, I1 and I2 are first calculated for an radial Z-dependet profile and then interpolated onto the 3D volume.
+ aplanatic:    aplanatic factor. Provided as a function of angle θ. Choices are `aplanatic_const`, `aplanatic_detection`, `aplanatic_illumination`, `aplanatic_illumination_flux`
+ FFTPlan:      information on how to calculate the FFTW plan. Default: nothing (using FFTW.ESTIMATE)
+ transition_dipole  If supplied a transition-dipole (e.g. sqrt(2) .* (0.0,1.0,1.0)) will be accounted for in the PSF calculation. If not normalized, the strength will be included.
Example:
```jdoctest
julia> using FFTW, PointSpreadFunctions

julia> ppm = PSFParams(0.58, 1.4, 1.518;pol=pol_circ,method=PointSpreadFunctions.MethodSincR, aberrations= Aberrations(), FFTPlan=FFTW.MEASURE)
PSFParams(0.58, 1.4, 1.518, Float32, ModeWidefield, PointSpreadFunctions.pol_scalar, PointSpreadFunctions.var"#42#43"(), PointSpreadFunctions.MethodSincR, 0x00000000, nothing, nothing)

julia> ppem = PSFParams(ppm; λ=0.620)

```
"""
function PSFParams(λ=0.5, NA=1.2, n=1.33; pol=pol_circ, dtype=Float32, mode=ModeWidefield, 
                    aplanatic = aplanatic_detection, method=MethodRichardsWolf, FFTPlan=nothing,
                    aberrations=Aberrations(), pixelshape=nothing, transition_dipole=nothing, λ_weights=nothing)
    PSFParams(λ, NA, n, dtype, mode, pol, aplanatic, method, FFTPlan, aberrations, pixelshape, transition_dipole, λ_weights)
end

"""
    PSFParams(pp; λ=pp.λ, NA=pp.NA, n=pp.n,  pol=pp.polarization, dtype=pp.dtype, mode=pp.mode, 
        aplanatic = pp.aplanatic, method=pp.method, FFTPlan=pp.FFTPlan,
        aberrations=pp.aberrations, pixelshape=pp.pixelshape, transition_dipole=pp.transition_dipole, λ_weights=pp.λ_weights)

    Alternative conveniance version to construct a PSF parameter structure using an existing structure `pp` and overwriting one or multiple parameters via named arguments.
    See the other version of `PSFParams` for documentation.
"""
function PSFParams(pp::PSFParams; λ=pp.λ, NA=pp.NA, n=pp.n,  pol=pp.polarization, dtype=pp.dtype, mode=pp.mode, 
    aplanatic = pp.aplanatic, method=pp.method, FFTPlan=pp.FFTPlan,
    aberrations=pp.aberrations, pixelshape=pp.pixelshape, transition_dipole=pp.transition_dipole, λ_weights=pp.λ_weights)
    PSFParams(λ, NA, n, dtype, mode, pol, aplanatic, method, FFTPlan, aberrations, pixelshape, transition_dipole, λ_weights)
end
