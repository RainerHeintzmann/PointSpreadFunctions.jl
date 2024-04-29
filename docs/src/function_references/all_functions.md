# List of Functions

## Parameters

see Workflow

## PointSpreadFunction generation

see Workflow

## Pupils

```@docs
PointSpreadFunctions.pupil_θ(sz, pp::PSFParams, sampling)
PointSpreadFunctions.pupil_ϕ(sz, pp::PSFParams, sampling)
aplanatic_factor(sz, pp::PSFParams, sampling)
pupil_xyz(sz, pp, sampling=nothing)
PointSpreadFunctions.field_xyz(sz, pp, sampling)
PointSpreadFunctions.field_xy_to_xyz(field,pp,sampling)
PointSpreadFunctions.field_pupil
PointSpreadFunctions.get_propagator(sz,pp,sampling)
PointSpreadFunctions.get_propagator_gradient(prop_phase, scalar, xy_scale)
PointSpreadFunctions.apply_propagators(pupil, z_planes, pp::PSFParams; sampling=nothing) 
PointSpreadFunctions.get_zernike_pupil_phase(sz, pp, sampling) 
PointSpreadFunctions.get_zernike_pupil(sz, pp, sampling) 
PointSpreadFunctions.get_Nyquist_limit
PointSpreadFunctions.xx_rfft
PointSpreadFunctions.yy_rfft
PointSpreadFunctions.rr_rfft
PointSpreadFunctions.to_rfft_pos  
PointSpreadFunctions.pinhole_AU_to_pix  
PointSpreadFunctions.AU_per_pixel   
PointSpreadFunctions.confocal_int
PointSpreadFunctions.jinc_r_2d
PointSpreadFunctions.normalize_amp_to_plane
PointSpreadFunctions.modify_STED_exp
PointSpreadFunctions.get_pupil_aperture
PointSpreadFunctions.kz_mid_pos
PointSpreadFunctions.disc_pinhole_ft
PointSpreadFunctions.modify_ident
PointSpreadFunctions.sinc_r_2d
PointSpreadFunctions.pol_circ_quadrant
PointSpreadFunctions.modify_STED_hyper
PointSpreadFunctions.box_pinhole_ft
PointSpreadFunctions.exp_ikx_rfft
PointSpreadFunctions.limit_θ
PointSpreadFunctions.size_sampling_to3d
PointSpreadFunctions.calc_with_resampling
PointSpreadFunctions.limit_theta
PointSpreadFunctions.pol_circ_tophat
PointSpreadFunctions.modify_square
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
PointSpreadFunctions.amp_to_int(field)
PointSpreadFunctions.has_z_symmetry(pp::PSFParams)
PointSpreadFunctions.get_Abbe_limit(pp::PSFParams)
PointSpreadFunctions.get_required_amp_sampling(sz::NTuple, pp::PSFParams)
PointSpreadFunctions.get_Ewald_sampling(sz::NTuple, pp::PSFParams)
PointSpreadFunctions.get_McCutchen_kz_center(sz, pp::PSFParams, sampling)
PointSpreadFunctions.limit_kz(ft_shell, pp::PSFParams, sampling)
PointSpreadFunctions.sinc_r(sz::NTuple, pp::PSFParams; sampling=nothing)
PointSpreadFunctions.jinc_r_2d(sz::NTuple, pp::PSFParams; sampling=nothing)
PointSpreadFunctions.my_disc(sz, pp)
PointSpreadFunctions.iftz(arr)
PointSpreadFunctions.theta_z(sz)
PointSpreadFunctions.k_0(pp::PSFParams)
PointSpreadFunctions.k_pupil(pp::PSFParams)
PointSpreadFunctions.k_dz(pp::PSFParams)
PointSpreadFunctions.k_scale(sz, pp::PSFParams, sampling)
PointSpreadFunctions.k_pupil_pos(sz, pp::PSFParams, sampling)
PointSpreadFunctions.k_0_pos(sz, pp::PSFParams, sampling)
PointSpreadFunctions.k_r(sz, pp::PSFParams, sampling)
PointSpreadFunctions.k_xy(sz,pp,sampling)
PointSpreadFunctions.k_xy_rel_pupil(sz,pp,sampling)
PointSpreadFunctions.check_amp_sampling_xy(sz, pp,sampling)
PointSpreadFunctions.check_amp_sampling_z(sz, pp,sampling)
PointSpreadFunctions.check_amp_sampling(sz, pp,sampling)
PointSpreadFunctions.check_amp_sampling_sincr(sz, pp,sampling)
```
