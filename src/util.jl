"""
    amp_to_int(field) 
    converts a complex-valued amplitude field to intensity via `abs2.` and summing over the 4th dimension.
"""
amp_to_int(field) = sum(abs2.(field), dims=4)


function get_sampling(sz::NTuple, pp::PSFParams)
    s_xyz = (pp.λ ./ pp.n) ./ 2
    sz2 = sz .÷ 2
    s_xyz .* (sz2.-2) ./ sz2 
end

function sinc_r(sz::NTuple, pp::PSFParams; sampling=nothing)
    if isnothing(sampling)
        sampling=get_sampling(sz, pp)
    end 
    sinc.(rr(pp.dtype, sz, scale=2 .*sampling ./ (pp.λ./pp.n)))
end

function jinc_r_2d(sz::NTuple, pp::PSFParams; sampling=nothing)
    if isnothing(sampling)
        sampling=get_sampling(sz, pp)
    end 
    jinc.(rr(pp.dtype, sz[1:2], scale=2 .*sampling[1:2] ./ (pp.λ./pp.NA)))
end

# my_disc(sz; rel_border=4/3, eps=0.05) = window_hanning(sz, border_in=rel_border.-eps, border_out=rel_border.+eps)
my_disc(sz, pp) = disc(pp.dtype, sz, sz .* (4/6))  # This is the radius where there is no overlap and the diagonals (sqrt(2)) are still covered

iftz(arr) = ift(arr,(3,))

theta_z(sz) = (zz(sz) .> 0) # The direction is important due to the highest frequency position at even-sized FFTs



"""
    k_0(pp::PSFParams)
    k in the medium as n/lambda
""" 
function k_0(pp::PSFParams)
    pp.dtype(pp.n / pp.λ)
end

"""
    k_pupil(pp::PSFParams)
    maxim radial k-coordinate where the pupil ends
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
+ `sampling`: pixelpitch in real space as NTuple
"""
function k_scale(sz, pp, sampling)
    pp.dtype.(1 ./ (sz .* sampling))
end

"""
    k_r(sz, pp::PSFParams, sampling)
    the radial k, |k_xy|
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

