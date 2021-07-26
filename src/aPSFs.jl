
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


function apsf(::Type{MethodSincR}, sz::NTuple, pp::PSFParams; sampling=nothing) 
    # The big_sz is used to avoid interference effects by the application of the multiplication with the pupil maks in Fourier space
    # The size if 4/3 is such that the diagonals are covered but a wrap-around into the original size is avoided.
    big_sz = ((sz[1:2].*2)...,sz[3]+2)  # The 2 extra slices in z seem to be necessary to avoid a problem at the fist positon
    pupil = pupil_xyz(big_sz, pp, sampling) # field_xyz(big_sz,pp, sampling) .* aplanatic_factor(big_sz,pp,sampling) .* ft(jinc_r_2d(big_sz[1:2],pp, sampling=sampling) .* my_disc(big_sz[1:2],pp)) # 
    sinc_r_big = sinc_r(big_sz,pp, sampling=sampling) .* my_disc(big_sz[1:2],pp)  # maybe this should rather already be apodized by angle?
    if isnothing(sampling)
        sampling = get_sampling(sz, pp)
        shell, sampling =  limit_kz(theta_z(big_sz) .* ft(sinc_r_big, (1,2,3)), pp, sampling)
    else
        shell =  theta_z(big_sz) .* ft(sinc_r_big, (1,2,3))
    end
    check_amp_sampling_sincr(big_sz, pp, sampling)
    print("Amplitude sampling is $sampling \n")
    pupils = cos.(pupil_θ(big_sz,pp,sampling)) .* pupil .* iftz(shell)
    res=ift2d(pupils) # This should really be a zoomed iFFT
    return NDTools.select_region(res, new_size=sz[1:3]) # extract the central bit, which avoids the wrap-around effects
    # , sampling
end


# This function is needed for the propagate method
function apply_propagators(pupil, z_planes, pp::PSFParams; sampling=nothing) 
    sz = (size(pupil)[1:2]...,z_planes,size(pupil)[4])
    # calculate the phase derivatives
    pupils = Array{Complex{pp.dtype}}(undef, sz)
    start_z = -z_planes÷2-1
    prop_phase, scalar, xy_scale = get_propagator(sz, pp, sampling)
    # dx_prop_phase, dy_prop_phase = get_propagator_gradient(prop_phase, scalar, xy_scale)
    for z in 1:sz[3]
        z_pos = z+start_z # reaches zero at center
        # z_pos_in_λ =  z_pos * sampling[3] / (pp.λ / pp.n)
        # pupils[:,:,z:z,:] .= pupil .* (prop.^z_pos_in_λ)
        # pupils[:,:,z:z,:] .= pupil .* cis.(z_pos .* prop_phase)

        # the integral of exp(i z prop(x)) over x is 
        # exp(i z prop(x))   (-2/(z^2) + i 2 prop(x)/z)
        pupils[:,:,z:z,:] .= pupil .* cis.(z_pos .* prop_phase) # .* 
             # (sinc.(z_pos .* dx_prop_phase) .* sinc.(z_pos .* dy_prop_phase))  # This dampling accounts for the aliasing problem
    end
    return pupils
end

function apsf(::Type{MethodPropagate}, sz::NTuple, pp::PSFParams; sampling=nothing) 
    if isnothing(sampling)
        sampling = get_sampling(sz, pp)
        print("Sampling is $sampling \n")
    end
    check_amp_sampling(sz, pp, sampling)
    pupil = pupil_xyz(sz, pp, sampling) # field_xyz(big_sz,pp, sampling) .* aplanatic_factor(big_sz,pp,sampling) .* ft(jinc_r_2d(big_sz[1:2],pp, sampling=sampling) .* my_disc(big_sz[1:2],pp)) # 
    pupils = apply_propagators(pupil, sz[3], pp, sampling=sampling)
    # return pupils
    # pupils .*= window_hanning((1,1,size(pupils,3)),border_in=0.8,border_out=1.0,dims=(3,))
    res=ift2d(pupils) # This should really be a zoomed iFFT
    return res # extract the central bit, which avoids the wrap-around effects
end

# This function is needed for the propagate method, but iterates between real and Fourier space always smoothly deleting "out-of-bound" waves.
# The final result is already in real space.
function apply_propagator_iteratively(pupil, z_planes, pp::PSFParams; sampling=nothing) 
    sz = (size(pupil)[1:2]...,z_planes,size(pupil)[4])
    # calculate the phase derivatives
    slices = Array{Complex{pp.dtype}}(undef, sz)
    start_z = z_planes÷2+1
    prop_phase, scalar, xy_scale = get_propagator(sz, pp, sampling)
    start_pupil = copy(pupil)
    prop_pupil = cis.(prop_phase)
    real_window = window_hanning(size(pupil)[1:2],border_in=0.9,border_out=1.0)
    prop_pupil = conj(prop_pupil) # from now the advancement is in the opposite direction
    for z in start_z:-1:2 # from the middle to the start
        slice = ift2d(pupil)
        slices[:,:,z:z,:] .= slice
        pupil = ft2d(slice .* real_window)
        pupil .*= prop_pupil 
    end
    slices[:,:,1:1,:] .= ift2d(pupil) 
    if has_z_symmetry(pp) # to save some speed
        dz = sz[3] - (start_z+1)
        slices[:,:,start_z+1:start_z+1+dz,:] .= conj(slices[:,:,start_z-1:-1:start_z-1-dz,:]);
    else
        prop_pupil = conj(prop_pupil) # from now the advancement is in the opposite direction
        pupil = start_pupil .* prop_pupil
        for z in start_z+1:sz[3]-1  # from the middle forward
            slice = ift2d(pupil)
            slices[:,:,z:z,:] .= slice
            pupil = ft2d(slice .* real_window)
            pupil .*= prop_pupil 
        end
        slices[:,:,sz[3]:sz[3],:] .= ift2d(pupil) 
    end
    return slices
end

function apsf(::Type{MethodPropagateIterative}, sz::NTuple, pp::PSFParams; sampling=nothing) 
    if isnothing(sampling)
        sampling = get_sampling(sz, pp)
        print("Sampling is $sampling \n")
    end
    check_amp_sampling(sz, pp, sampling)
    pupil = pupil_xyz(sz, pp, sampling) # field_xyz(big_sz,pp, sampling) .* aplanatic_factor(big_sz,pp,sampling) .* ft(jinc_r_2d(big_sz[1:2],pp, sampling=sampling) .* my_disc(big_sz[1:2],pp)) # 
    slices = apply_propagator_iteratively(pupil, sz[3], pp, sampling=sampling)
    return slices # extract the central bit, which avoids the wrap-around effects
end

# just a different way of writing the propagation down
function apply_kz_to_pupil(pupil, z_planes, pp; sampling)
    sz = (size(pupil)[1:2]...,z_planes,size(pupil)[4])
    if isnothing(sampling)
        sampling = get_sampling(sz, pp)
    end
    k_max_rel = sampling[1:2] ./ (pp.λ / pp.n)
    scalar = (2π*sampling[3] / (pp.λ / pp.n))
    xy_scale = 1 ./ (k_max_rel .* sz[1:2])
    sqrt_term = phase_kz(pp.dtype, sz[1:2], scale = xy_scale)
    prop_phase = scalar .* sqrt_term    
    pupils = pupil .* cis.(zz((1,1,z_planes)) .* prop_phase) # if you fourier-transform this along kz you should get the shifted pupil
    return pupils
end


# Basic idea of MethodShell: Displace each kxy by the appropriate kz 
function apsf(::Type{MethodShell}, sz::NTuple, pp::PSFParams; sampling=nothing) 
    if isnothing(sampling)
        sampling = get_sampling(sz, pp)
        print("Sampling is $sampling \n")
    end
    check_amp_sampling(sz, pp, sampling)
    pupil = pupil_xyz(sz, pp, sampling) # field_xyz(big_sz,pp, sampling) .* aplanatic_factor(big_sz,pp,sampling) .* ft(jinc_r_2d(big_sz[1:2],pp, sampling=sampling) .* my_disc(big_sz[1:2],pp)) # 

    pupils = apply_kz_to_pupil(pupil, sz[3], pp, sampling=sampling)    
    res=ift2d(pupils) # This should really be a zoomed iFFT
    return res # extract the central bit, which avoids the wrap-around effects
end

"""
    apsf(sz::NTuple, pp::PSFParams; sampling=get_sampling(sz, pp))
Example:
pp = PSFParams(500.0,1.4,1.52)
p = apsf((128,128,128),pp; sampling=(50,50,100)); #; 
"""
function apsf(sz::NTuple, pp::PSFParams; sampling=nothing) 
    apsf(pp.method, sz, pp, sampling=sampling)
end

