# PSfs.jl Documentation
This toolbox aims at fast and accurate calculations of point spread functions (PSFs) as used in optics. It is particularyl
suited for high numerical aperture PSfs and vectorial aspects of the electric field. 
Aplanatic factors as for aplanatic optical systems fullfilling the Abbe sine condition are accounted for as well as
various pupil properties can be supplied. 
Different modes of calculation are possible and more will be added in the future. Aberrations can be specified in terms of Zernike coefficients.

## PSF calculation

```@docs
PSFParams
psf
apsf
```

## Aberrations

```@docs
Aberrations
```

## Pupils

```@docs
PSFs.pupil_θ(sz, pp::PSFParams, sampling)
PSFs.pupil_ϕ(sz, pp::PSFParams, sampling)
aplanatic_factor(sz, pp::PSFParams, sampling)
pupil_xyz(sz, pp, sampling=nothing)
PSFs.field_xyz(sz, pp, sampling)
PSFs.field_xy_to_xyz(field,pp,sampling)
PSFs.get_propagator(sz,pp,sampling)
PSFs.get_propagator_gradient(prop_phase, scalar, xy_scale)
PSFs.get_zernike_pupil_phase(sz, pp, sampling, J, coefficients; index) 
PSFs.get_zernike_pupil(sz, pp, sampling) 
```

## Polarization
These functions can be conveniently supplied to `PSFParams()` via the named argument `polarization` 
```@docs
pol_scalar
pol_scalar_spiral
pol_x
pol_y
pol_circ
pol_circ_spiral
```

## Aplanatic factors
```@docs
aplanatic_detection
aplanatic_illumination
aplanatic_const
aplanatic_illumination_flux
```
