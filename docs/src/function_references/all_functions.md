# List of Functions

## Parameters

see Workflow

## PSF generation

see Workflow

## Pupils

```@docs
PSFs.pupil_θ(sz, pp::PSFParams, sampling)
PSFs.pupil_ϕ(sz, pp::PSFParams, sampling)
aplanatic_factor(sz, pp::PSFParams, sampling)
pupil_xyz(sz, pp, sampling=nothing)
PSFs.field_xyz(sz, pp, sampling)
PSFs.field_xy_to_xyz(field,pp,sampling)
PSFs.field_pupil
PSFs.get_propagator(sz,pp,sampling)
PSFs.get_propagator_gradient(prop_phase, scalar, xy_scale)
PSFs.apply_propagators(pupil, z_planes, pp::PSFParams; sampling=nothing) 
PSFs.get_zernike_pupil_phase(sz, pp, sampling) 
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

## Utilities

```@docs
PSFs.amp_to_int(field)
PSFs.has_z_symmetry(pp::PSFParams)
PSFs.get_Abbe_limit(pp::PSFParams)
PSFs.get_required_amp_sampling(sz::NTuple, pp::PSFParams)
PSFs.get_Ewald_sampling(sz::NTuple, pp::PSFParams)
PSFs.get_McCutchen_kz_center(sz, pp::PSFParams, sampling)
PSFs.limit_kz(ft_shell, pp::PSFParams, sampling)
PSFs.sinc_r(sz::NTuple, pp::PSFParams; sampling=nothing)
PSFs.jinc_r_2d(sz::NTuple, pp::PSFParams; sampling=nothing)
PSFs.my_disc(sz, pp)
PSFs.iftz(arr)
PSFs.theta_z(sz)
PSFs.k_0(pp::PSFParams)
PSFs.k_pupil(pp::PSFParams)
PSFs.k_dz(pp::PSFParams)
PSFs.k_scale(sz, pp::PSFParams, sampling)
PSFs.k_pupil_pos(sz, pp::PSFParams, sampling)
PSFs.k_0_pos(sz, pp::PSFParams, sampling)
PSFs.k_r(sz, pp::PSFParams, sampling)
PSFs.k_xy(sz,pp,sampling)
PSFs.k_xy_rel_pupil(sz,pp,sampling)
PSFs.check_amp_sampling_xy(sz, pp,sampling)
PSFs.check_amp_sampling_z(sz, pp,sampling)
PSFs.check_amp_sampling(sz, pp,sampling)
PSFs.check_amp_sampling_sincr(sz, pp,sampling)
```
