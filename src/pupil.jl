
function pupil_θ(sz, pp::PSFParams, sampling)
    asin.(k_r(sz, pp, sampling)./k_0(pp))
end

function pupil_ϕ(sz, pp::PSFParams, sampling)
    ϕ_tuple.(k_xy(sz, pp, sampling))
end

function aplanatic_factor(sz, pp::PSFParams, sampling)
    pp.aplanatic.(pupil_θ(sz, pp, sampling))
end

function field_pupil(sz, pp, sampling)
    idx_to_dim(pp.polarization.(pp.dtype, k_xy_rel_pupil(sz,pp,sampling)), 4)
end

"""
    field_xyz(field, pp, sampling)
    converts a 2D field at the pupil to a 2D field behind the lens, containing E_x, E_Y and E_z.
"""
function field_xy_to_xyz(field,pp,sampling)
    if size(field)[4] == 1
        return field; # just stay scalar and apply nothing.
    end
    sz = size(field)
    ϕ = pupil_ϕ(sz, pp, sampling) # establish angles
    θ = pupil_θ(sz, pp, sampling)
    cosϕ = cos.(ϕ)
    sinϕ = sin.(ϕ)    
    cosθ = cos.(θ)
    sinθ = sin.(θ)
    field_rad = field[:,:,:,1] .* cosϕ + field[:,:,:,2] .* sinϕ # radial field component in the pupil
    field_azi = -field[:,:,:,1] .* sinϕ + field[:,:,:,2] .* cosϕ # azimuthal field component in the pupil
    field_rad_cosθ = cosθ .* field_rad # what remains of the radial field component behind lens (with energy conservation)
    cat(field_rad_cosθ .* cosϕ .- field_azi .* sinϕ, field_rad_cosθ .* sinϕ .+ field_azi .* cosϕ, sinθ .* field_rad, dims=4)
end

"""
    field_xyz(sz, pp, sampling)
    creates a 2D pupil field behind the lens, containing E_x, E_Y and E_z.
"""
function field_xyz(sz, pp, sampling)
    field_xy_to_xyz(field_pupil(sz, pp, sampling), pp, sampling)
end

"""
    pupil_xyz(sz, pp, sampling=nothing)
    creates a pupil with electric field distributions in XYZ. Returns a 4D dataset with the electric field components along the 4th dimension.
"""
function pupil_xyz(sz, pp, sampling=nothing)
    if isnothing(sampling)
        sampling = get_Ewald_sampling(sz, pp)
    end
    field_xyz(sz, pp, sampling) .* aplanatic_factor(sz,pp,sampling) .* ft(jinc_r_2d(sz[1:2], pp, sampling=sampling) .* my_disc(sz[1:2],pp))
end

function get_propagator(sz,pp,sampling)
    if isnothing(sampling)
        sampling = get_Ewald_sampling(sz, pp)
    end
    k_max_rel = sampling[1:2] ./ (pp.λ / pp.n)
    # prop1 = propagator(pp.dtype, sz[1:2], Δz=1, k_max = k_max_rel)
    # prop = cispi.(2 .* phase_kz(pp.dtype, sz[1:2], scale = 1 ./ (k_max_rel .* sz[1:2])))
    scalar = (2π*sampling[3] / (pp.λ / pp.n))
    xy_scale = 1 ./ (k_max_rel .* sz[1:2])
    return scalar .* phase_kz(pp.dtype, sz[1:2], scale = xy_scale), scalar, xy_scale
end

"""

"""
function get_propagator_gradient(prop_phase, scalar, xy_scale)
    #dr_prop_phase = ifelse.(prop_phase .== zero(pp.dtype), zero(pp.dtype), scalar.^2 ./ prop_term .* rr(pp.dtype, sz[1:2], scale = xy_scale))
    dtype = eltype(prop_phase)
    sz = size(prop_phase)[1:2]
    dx_prop_phase = ifelse.(prop_phase .== zero(dtype), zero(dtype), scalar.^2 ./ prop_phase .* xx(dtype, sz, scale = xy_scale.^2))
    dy_prop_phase = ifelse.(prop_phase .== zero(dtype), zero(dtype), scalar.^2 ./ prop_phase .* yy(dtype, sz, scale = xy_scale.^2))
    return dx_prop_phase, dy_prop_phase
end

