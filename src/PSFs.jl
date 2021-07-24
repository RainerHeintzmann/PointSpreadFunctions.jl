module PSFs
using FourierTools: center_pos
using FourierTools, NDTools, IndexFunArrays, SpecialFunctions
export PSFParams, sinc_r, jinc_r_2d, pupil_xyz, apsf, psf, k0, kxy, aplanatic_factor
export ModeWidefield, ModeConfocal, Mode4Pi

abstract type PSFMode end
struct ModeWidefield <: PSFMode end
struct ModeConfocal <: PSFMode end
struct Mode4Pi <: PSFMode end

struct PSFParams
    λ
    NA
    n
    dtype # real-valued data type to generate PSF for
    mode # microscopy mode to calculate PSF for   ::PSFMode
    polarization # a function calculating the polarization from a given Tuple of relative-k pupil vectors
    aplanatic # aplanatic factor
end

function PSFParams(my_λ, my_NA, my_n)
    PSFParams(my_λ, my_NA, my_n, Float32, ModeWidefield, pol_x, (θ) -> sqrt.(max.(0,cos.(θ))))
end

# define a bunch of standard pupil polarizations
function pol_x(T, xypos)
    (one(T),zero(T))
end

function pol_y(T, xypos)
    (zero(T), one(T))
end

function pol_circ(T, xypos)
    (one(T)/sqrt(2)+0im, one(T)/sqrt(2)*1im)
end

function pol_circ_spiral(T, xypos) # e.g. for STED microscopy
    pol_circ(T, xypos) .* cis.(atan(xypos...))
end


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

"""
    limit_k_z(ft_shell)
    limits the k_z range of the ewald-sphere.
"""
function limit_kz(ft_shell, pp, sampling)
    sz = size(ft_shell)
    k_z_scale = k_scale(sz, pp, sampling)[3]
    dkz = k_dz(pp) ./ k_z_scale
    pk0 = k_0(pp) ./ k_z_scale
    old_center = sz .÷ 2 .+ 1
    new_center = (old_center[1], old_center[2], round(eltype(old_center), old_center[3] .+ pk0 .- dkz /2))
    new_size = (sz[1], sz[2], ceil(eltype(sz), dkz).+1) 
    cut_shell = NDTools.select_region(ft_shell, center = new_center, new_size = new_size)
    sampling = (sampling[1],sampling[2],sampling[2] * sz[3] / new_size[3])
    return cut_shell, sampling 
end

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

function pupil_xyz(sz, pp, sampling=nothing)
    if isnothing(sampling)
        sampling = get_sampling(sz, pp)
    end
    field_xyz(sz, pp, sampling) .* aplanatic_factor(sz,pp,sampling) .* ft(jinc_r_2d(sz[1:2], pp, sampling=sampling) .* my_disc(sz[1:2],pp))
end

"""
    apsf(sz::NTuple, pp::PSFParams; sampling=get_sampling(sz, pp))
Example:
pp = PSFParams(500.0,1.4,1.52)
p = apsf((128,128,128),pp; sampling=(50,50,100)); #; 
"""
function apsf(sz::NTuple, pp::PSFParams; sampling=nothing)

    big_sz = ((sz[1:2].*2)...,sz[3]+2)  # The 2 extra slices in z seem to be necessary to avoid a problem at the fist positon
    # disc_rad = sz[1:2] .* (4/6) # This is the radius where there is no overlap and the diagonals (sqrt(2)) are still covered
    # pupil = ft(NDTools.select_region(jinc_r_2d(ip,pp) .* disc(ip.sz[1:2],ip.sz[1:2]./2),new_size=ip.sz[1:2].*2))  # disc(ip.sz[1:2],ip.sz[1:2]./2)
    # return cos.(pupil_θ(big_sz,pp,sampling))
    pupil = pupil_xyz(big_sz, pp, sampling) # field_xyz(big_sz,pp, sampling) .* aplanatic_factor(big_sz,pp,sampling) .* ft(jinc_r_2d(big_sz[1:2],pp, sampling=sampling) .* my_disc(big_sz[1:2],pp)) # 
    sinc_r_big = sinc_r(big_sz,pp, sampling=sampling) .* my_disc(big_sz[1:2],pp)  # maybe this should rather already be apodized by angle?
    # shell =  theta_z(big_sz) .* ft(sinc_r_big, (1,2,3)
    if isnothing(sampling)
            sampling = get_sampling(sz, pp)
            shell, sampling =  limit_kz(theta_z(big_sz) .* ft(sinc_r_big, (1,2,3)), pp, sampling)
        else
        shell =  theta_z(big_sz) .* ft(sinc_r_big, (1,2,3))
    end
    print("sampling is $sampling \n")
    pupils = cos.(pupil_θ(big_sz,pp,sampling)) .* pupil .* iftz(shell)
    res=ift2d(pupils) # This should really be a zoomed iFFT
    return NDTools.select_region(res, new_size=sz[1:2]), sampling # extract the central bit, which avoids the wrap-around effects
end

amp_to_int(field) = sum(abs2.(field), dims=4)

"""
    psf(sz::NTuple, pp::PSFParams; sampling=get_sampling(sz, pp))
Example:
pp = PSFParams(500.0,1.4,1.52)
p = psf((128,128,128),pp; sampling=(50,50,100)); #; 
"""
function psf(sz::NTuple, pp::PSFParams; sampling=nothing)
    amp = apsf(sz, pp, sampling=sampling)
    if isnothing(sampling)
        sum(upsample2_abs2(amp),dims=4)
    else
        amp_to_int(amp)
    end
end

end # module
