"""
    amp_to_int(field) 
    converts a complex-valued amplitude field to intensity via `abs2.` and summing over the 4th dimension.
"""
amp_to_int(field) = sum(abs2.(field), dims=4)

function has_z_symmetry(pp::PSFParams)
    return true; # will be changed later, when assymetric aberrations are allowed
end

function get_Abbe_limit(pp)
    d_xy = pp.λ ./ (2 .* pp.n) 
    d_z = (1 - cos(asin(pp.NA/ pp.n))) * pp.λ / pp.n
    pp.dtype.((d_xy,d_xy,d_z))
end

function get_required_amp_sampling(sz::NTuple, pp::PSFParams)
    abbe = get_Abbe_limit(pp)[1:length(sz)]
    sz2 = sz .÷ 2
    abbe .* (sz2.-2) ./ sz2 # provide a minimum amount of oversampling to avoid problems with the border pixesl.
end

"""
    get_Ewald_sampling(sz::NTuple, pp::PSFParams)
    returns the required minimum sampling for the calculation of a full Ewald sphere.
"""
function get_Ewald_sampling(sz::NTuple, pp::PSFParams)
    s_xyz = pp.λ ./ (2 .* pp.n) 
    sz2 = sz .÷ 2
    s_xyz .* (sz2.-2) ./ sz2 # provide a minimum amount of oversampling to avoid problems with the border pixesl.
end


"""
    get_McCutchen_kz_center(ft_shell, pp, sampling)
    calculates the (rounded) pixels position half way between both, the kz-borders of the McCutchen pupil to extract from the full sized Ewald sphere.
    The pixel z position is returned together with the corresponding kz position.
"""
function get_McCutchen_kz_center(sz, pp, sampling)
    k_z_scale = k_scale(sz, pp, sampling)[3]
    dkz = k_dz(pp) ./ k_z_scale
    pk0 = k_0(pp) ./ k_z_scale
    old_center = sz .÷ 2 .+ 1
    new_center = (old_center[1], old_center[2], round(eltype(old_center), old_center[3] .+ pk0 .- dkz /2))
    # kz_center = (new_center[3] - old_center[3]) * k_z_scale
    return new_center, new_center[3]-old_center[3]
end

"""
    limit_k_z(ft_shell, pp, sampling)
    limits the k_z range of the ewald-sphere.
    returns the extracted region and the new sampling
"""
function limit_kz(ft_shell, pp, sampling)
    sz = size(ft_shell)
    get_kz_center(sz, pp, sampling)
    dkz = k_dz(pp) ./ k_z_scale
    new_size = (sz[1], sz[2], ceil(eltype(sz), dkz).+1) 
    cut_shell = FourierTools.select_region_ft(ft_shell, center = new_center, new_size = new_size)
    sampling = (sampling[1],sampling[2],sampling[2] * sz[3] / new_size[3])
    return cut_shell, sampling 
end

function sinc_r(sz::NTuple, pp::PSFParams; sampling=nothing)
    if isnothing(sampling)
        sampling=get_Ewald_sampling(sz, pp)
    end 
    sinc.(rr(pp.dtype, sz, scale=2 .*sampling ./ (pp.λ./pp.n)))
end

function jinc_r_2d(sz::NTuple, pp::PSFParams; sampling=nothing)
    if isnothing(sampling)
        sampling=get_Ewald_sampling(sz, pp)
    end 
    jinc.(rr(pp.dtype, sz[1:2], scale=2 .*sampling[1:2] ./ (pp.λ./pp.NA)))
end

# my_disc(sz; rel_border=4/3, eps=0.05) = window_hanning(sz, border_in=rel_border.-eps, border_out=rel_border.+eps)
my_disc(sz, pp) = disc(pp.dtype, sz, sz .* (4/6))  # This is the radius where there is no overlap and the diagonals (sqrt(2)) are still covered

iftz(arr) = ift(arr,(3,))

theta_z(sz) = (zz(sz) .> 0) # The direction is important due to the highest frequency position at even-sized FFTs



"""
    k_0(pp::PSFParams)
    k in the medium as n/lambda.   (1/space units
""" 
function k_0(pp::PSFParams)
    pp.dtype(pp.n / pp.λ)
end

"""
    k_pupil(pp::PSFParams)
    maxim radial k-coordinate (1/space units) where the pupil ends
+ `pp`:  PSF parameter structure
"""
function k_pupil(pp::PSFParams)
    pp.dtype(pp.NA / pp.λ)
end

"""
    k_dz(pp::PSFParams)
    relative kz range from pupil boarder to top of Ewald sphere
+ `pp`:  PSF parameter structure
"""
function k_dz(pp::PSFParams)
    pp.dtype((1 - cos(asin(pp.NA/ pp.n))) * k_0(pp))
end

"""
    k_scale(sz,sampling)
    pixelpitch (as NTuple) in k-space
+ `sz`:  size of the real-space array
+ `pp`:  PSF parameter structure
+ `sampling`: pixelpitch in real space as NTuple
"""
function k_scale(sz, pp::PSFParams, sampling)
    pp.dtype.(1 ./ (sz .* sampling))
end

"""
    k_pupil_pos(sz, pp::PSFParams, sampling)
    returns the X and Y position of the pupil border in reciprocal space pixels.
+ `sz`:  size of the real-space array
+ `pp`:  PSF parameter structure
+ `sampling`: pixelpitch in real space as NTuple
"""
function k_pupil_pos(sz, pp::PSFParams, sampling)
    k_pupil(pp) ./ k_scale(sz, pp, sampling)
end

"""
    k_0_pos(sz, pp::PSFParams, sampling)
    returns the X and Y position of the Ewald-sphere border in reciprocal space pixels.
+ `sz`:  size of the real-space array
+ `pp`:  PSF parameter structure
+ `sampling`: pixelpitch in real space as NTuple
"""
function k_0_pos(sz, pp::PSFParams, sampling)
    k_0(pp) ./ k_scale(sz, pp, sampling)
end

"""
    k_r(sz, pp::PSFParams, sampling)
    returns an array of radial k coordinates, |k_xy|
"""
function k_r(sz, pp::PSFParams, sampling)
    min.(k_0(pp), rr(pp.dtype, sz[1:2],scale = k_scale(sz[1:2], pp, sampling[1:2])))
end

"""
    k_xy(sz,pp,sampling)
    yields a 2D array with each entry being a 2D Tuple.
"""
function k_xy(sz,pp,sampling)
    idx(pp.dtype, sz[1:2],scale = k_scale(sz[1:2], pp, sampling[1:2]))
end

"""
    k_xy_rel_pupil(sz,pp,sampling)
    returns an array of relative distance to the pupil border
"""
function k_xy_rel_pupil(sz,pp,sampling)
    idx(pp.dtype, sz[1:2],scale = k_scale(sz[1:2], pp, sampling[1:2]) ./ k_pupil(pp))
end

function check_amp_sampling_xy(sz, pp,sampling)
    sample_factor = k_pupil(pp) ./ ((sz[1:2] .÷2) .* k_scale(sz, pp, sampling)[1:2])
    if any(sample_factor .> 1.0)
        @warn "Your calculation is undersampled along XY by factors of $sample_factor. The PSF calculation will be incorrect.)"
    end
end

function check_amp_sampling_z(sz, pp,sampling)
    sample_factor = k_dz(pp) ./ ((sz[3] .÷2) .* k_scale(sz, pp, sampling)[3])
    if (sample_factor > 1.0)
        @warn "Your calculation is undersampled along Z by factors of $sample_factor. The PSF calculation will be incorrect.)"
    end
end

function check_amp_sampling(sz, pp,sampling)
    check_amp_sampling_xy(sz, pp, sampling)
    check_amp_sampling_z(sz, pp, sampling)
end

function check_amp_sampling_sincr(sz, pp,sampling) # The sinc-r method needs (for now without aliasing) to be sampled extremely high along Z
    check_amp_sampling_xy(sz, pp, sampling)
    @show sample_factor = k_0(pp) ./ ((sz[3] .÷2) .* k_scale(sz, pp, sampling)[3])
    if (sample_factor > 1.0)
        @warn "Your calculation is undersampled along Z by factors of $sample_factor. The PSF calculation will be incorrect.)"
    end
end

