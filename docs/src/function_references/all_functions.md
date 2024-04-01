# List of Functions

## Parameters

see Workflow

## PointSpreadFunction generation

see Workflow

## Pupils

```@docs
PointSpreadFunctions.pupil_θ(sz, pp::PointSpreadFunctionParams, sampling)
PointSpreadFunctions.pupil_ϕ(sz, pp::PointSpreadFunctionParams, sampling)
aplanatic_factor(sz, pp::PointSpreadFunctionParams, sampling)
pupil_xyz(sz, pp, sampling=nothing)
PointSpreadFunctions.field_xyz(sz, pp, sampling)
PointSpreadFunctions.field_xy_to_xyz(field,pp,sampling)
PointSpreadFunctions.field_pupil
PointSpreadFunctions.get_propagator(sz,pp,sampling)
PointSpreadFunctions.get_propagator_gradient(prop_phase, scalar, xy_scale)
PointSpreadFunctions.apply_propagators(pupil, z_planes, pp::PointSpreadFunctionParams; sampling=nothing) 
PointSpreadFunctions.get_zernike_pupil_phase(sz, pp, sampling) 
PointSpreadFunctions.get_zernike_pupil(sz, pp, sampling) 
```

## Polarization

These functions can be conveniently supplied to `PointSpreadFunctionParams()` via the named argument `polarization` 
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
PointSpreadFunctions.amp_to_int(field)
PointSpreadFunctions.has_z_symmetry(pp::PointSpreadFunctionParams)
PointSpreadFunctions.get_Abbe_limit(pp::PointSpreadFunctionParams)
PointSpreadFunctions.get_required_amp_sampling(sz::NTuple, pp::PointSpreadFunctionParams)
PointSpreadFunctions.get_Ewald_sampling(sz::NTuple, pp::PointSpreadFunctionParams)
PointSpreadFunctions.get_McCutchen_kz_center(sz, pp::PointSpreadFunctionParams, sampling)
PointSpreadFunctions.limit_kz(ft_shell, pp::PointSpreadFunctionParams, sampling)
PointSpreadFunctions.sinc_r(sz::NTuple, pp::PointSpreadFunctionParams; sampling=nothing)
PointSpreadFunctions.jinc_r_2d(sz::NTuple, pp::PointSpreadFunctionParams; sampling=nothing)
PointSpreadFunctions.my_disc(sz, pp)
PointSpreadFunctions.iftz(arr)
PointSpreadFunctions.theta_z(sz)
PointSpreadFunctions.k_0(pp::PointSpreadFunctionParams)
PointSpreadFunctions.k_pupil(pp::PointSpreadFunctionParams)
PointSpreadFunctions.k_dz(pp::PointSpreadFunctionParams)
PointSpreadFunctions.k_scale(sz, pp::PointSpreadFunctionParams, sampling)
PointSpreadFunctions.k_pupil_pos(sz, pp::PointSpreadFunctionParams, sampling)
PointSpreadFunctions.k_0_pos(sz, pp::PointSpreadFunctionParams, sampling)
PointSpreadFunctions.k_r(sz, pp::PointSpreadFunctionParams, sampling)
PointSpreadFunctions.k_xy(sz,pp,sampling)
PointSpreadFunctions.k_xy_rel_pupil(sz,pp,sampling)
PointSpreadFunctions.check_amp_sampling_xy(sz, pp,sampling)
PointSpreadFunctions.check_amp_sampling_z(sz, pp,sampling)
PointSpreadFunctions.check_amp_sampling(sz, pp,sampling)
PointSpreadFunctions.check_amp_sampling_sincr(sz, pp,sampling)
```
