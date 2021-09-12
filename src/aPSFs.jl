
function apsf(::Type{MethodSincR}, sz::NTuple, pp::PSFParams; sampling=nothing, center_kz=false) 
    check_amp_sampling(sz, pp, sampling)

    if isnothing(sampling)
        sampling = get_required_amp_sampling(sz, pp) # this is knowingly too small along kz but fixed down below, if needed.
    end
    Ewald_sampling = get_Ewald_sampling(sz, pp)
    undersampling_factor = sampling[3] ./ Ewald_sampling[3]
    # The big_sz is used to avoid interference effects by the application of the multiplication with the pupil maks in Fourier space
    # The size if 4/3 is such that the diagonals are covered but a wrap-around into the original size is avoided.
    # big_sz = ((sz[1:2].*2)...,sz[3]+2)  # The 2 extra slices in z seem to be necessary to avoid a problem at the fist position
    big_sz = ((sz[1:2].*2)...,sz[3]+4)  # The 2 extra slices in z seem to be necessary to avoid a problem at the fist positions

    if undersampling_factor > 1.0 # the current sampling is insufficient for the SincR method. We need to properly sample the Ewald sphere.
        # the size sz[3] will be cut out in Fourier-spage around the cleaned and processed part of the McCutchen pupil, which should yield the sampling as currently specified.
        nowrap_sz = (big_sz[1:2]..., ceil(Int, undersampling_factor .* big_sz[3]))
        big_sampling = (sampling[1:2]..., sampling[3] * big_sz[3] / nowrap_sz[3]) # calculate the sampling needed to cover the Ewald sphere.
    else
        nowrap_sz = big_sz
        big_sampling = sampling
    end

    pupil = pupil_xyz(nowrap_sz, pp, big_sampling) # field_xyz(big_sz,pp, sampling) .* aplanatic_factor(big_sz,pp,sampling) .* ft(jinc_r_2d(big_sz[1:2],pp, sampling=sampling) .* my_disc(big_sz[1:2],pp)) # 
    sinc_r_big = sinc_r(nowrap_sz,pp, sampling=big_sampling) .* my_disc(nowrap_sz[1:2],pp)  # maybe this should rather already be apodized by angle?
    if pp.FFTPlan != nothing
        P3d = plan_fft(sinc_r_big, flags=pp.FFTPlan)
    end

    if big_sampling[3] != sampling[3] # we need to extract (to reduce the size to the user-defined) and circshift (to get the correct phases)
        kzc, rel_kz = get_McCutchen_kz_center(nowrap_sz,pp,big_sampling)
        shell = select_region(theta_z(nowrap_sz) .* ft(sinc_r_big, (1,2,3)), new_size=(nowrap_sz[1:2]...,big_sz[3]), center = kzc) # centers the Pupil along kz
        if !center_kz
            shell = circshift(shell, (0,0,rel_kz)) # 
        end
    else
        shell =  theta_z(nowrap_sz) .* ft(sinc_r_big, (1,2,3))
        if center_kz
            kzc, rel_kz = get_McCutchen_kz_center(nowrap_sz,pp,big_sampling)
            shell = circshift(shell, (0,0,-rel_kz)) # needed for subsequent upsampling along kz
        end
    end
    # shell, sampling =  limit_kz(theta_z(nowrap_sz) .* ft(sinc_r_big, (1,2,3)), pp, sampling) # remove negative frequencies and limit to useful range

    # check_amp_sampling_sincr(nowrap_sz, pp, sampling)
    # print("Amplitude sampling is $sampling \n")
    pupils = (cos.(pupil_θ(nowrap_sz,pp,sampling)) .* pupil) .* iftz(shell)
    if pp.FFTPlan != nothing
        Pm2d = plan_fft(pupils,(1,2), flags=pp.FFTPlan)
    end

    res=ift2d(pupils) # This should really be a zoomed iFFT
    return select_region(res, new_size=sz[1:3]) # extract the central bit, which avoids the wrap-around effects
    # , sampling
end


# This function is needed for the propagate method
"""
    apply_propagators(pupil, z_planes, pp::PSFParams; sampling=nothing) 

propagates a given `pupil` by a number of `z_planes` (almost) symmetrically in both directions.
The result is a stack of propagated pupils. The slice-to-slice propagator is obtained via the
`get_propagator` method.

See also:
+ get_propagator

"""
function apply_propagators(pupil, z_planes, pp::PSFParams; sampling=nothing) 
    sz = (size(pupil)[1:2]...,z_planes,size(pupil)[4])
    # calculate the phase derivatives
    pupils = Array{Complex{pp.dtype}}(undef, sz)
    start_z = -z_planes÷2-1
    prop_kz , scalar, xy_scale = get_propagator(sz, pp, sampling)
    # dx_prop_phase, dy_prop_phase = get_propagator_gradient(prop_phase, scalar, xy_scale)
    for z in 1:sz[3]
        z_pos = z+start_z # reaches zero at center
        # z_pos_in_λ =  z_pos * sampling[3] / (pp.λ / pp.n)
        # pupils[:,:,z:z,:] .= pupil .* (prop.^z_pos_in_λ)
        # pupils[:,:,z:z,:] .= pupil .* cis.(z_pos .* prop_phase)

        # the integral of exp(i z prop(x)) over x is 
        # exp(i z prop(x))   (-2/(z^2) + i 2 prop(x)/z)
        pupils[:,:,z:z,:] .= pupil .* cis.(z_pos .* prop_kz) # .*  # cis means: exp.(i * z_pos .* prop_phase)
             # (sinc.(z_pos .* dx_prop_phase) .* sinc.(z_pos .* dy_prop_phase))  # This dampling accounts for the aliasing problem
    end
    return pupils
end


function apsf(::Type{MethodPropagate}, sz::NTuple, pp::PSFParams; sampling=nothing, center_kz=false) 
    if isnothing(sampling)
        sampling = get_Ewald_sampling(sz, pp)
        print("Sampling is $sampling \n")
    end
    check_amp_sampling(sz, pp, sampling)
    pupil = pupil_xyz(sz, pp, sampling) # field_xyz(big_sz,pp, sampling) .* aplanatic_factor(big_sz,pp,sampling) .* ft(jinc_r_2d(big_sz[1:2],pp, sampling=sampling) .* my_disc(big_sz[1:2],pp)) # 
    pupils = apply_propagators(pupil, sz[3], pp, sampling=sampling)
    # return pupils
    # pupils .*= window_hanning((1,1,size(pupils,3)),border_in=0.8,border_out=1.0,dims=(3,))
    if center_kz
        _, rel_kz = get_McCutchen_kz_center(sz,pp,sampling)
        pupils .*= cispi.((-2*rel_kz/sz[3]) .* zz((1,1,sz[3]))) # centers the McCutchen pupil to be able to correctly resample it along kz
    end
    if pp.FFTPlan != nothing
        Pm2d = plan_fft(pupils,(1,2), flags=pp.FFTPlan)
    end
    res=ift2d(pupils) # This should really be a zoomed iFFT
    return res # extract the central bit, which avoids the wrap-around effects
end

# This function is needed for the propagate method, but iterates between real and Fourier space always smoothly deleting "out-of-bound" waves.
# The final result is already in real space.
function apply_propagator_iteratively(sz, pp::PSFParams; sampling=nothing, center_kz=false) 
    z_planes = sz[3]
    # calculate the phase derivatives
    start_z = z_planes÷2+1
    max_pix_travel = (tan(asin(pp.NA / pp.n)) * sampling[3]) ./ sampling[1:2] # how much does the maximal anlge travel geometrically
    PMLsz = (ceil.(Int, max_pix_travel .* 8)...,0) # To get the number of perfectly matched layer (PML) pixels to append
    psz = sz .+ 2 .*PMLsz # total size for propagation pupil
    pupil = pupil_xyz(psz[1:2], pp, sampling) # field_xyz(big_sz,pp, sampling) .* aplanatic_factor(big_sz,pp,sampling) .* ft(jinc_r_2d(big_sz[1:2],pp, sampling=sampling) .* my_disc(big_sz[1:2],pp)) # 
    sz = (sz[1:3]...,size(pupil)[4])
    slices = Array{Complex{pp.dtype}}(undef, sz)

    prop_kz, scalar, xy_scale = get_propagator(psz, pp, sampling)
    start_pupil = pupil # copy(pupil)
    prop_pupil = cis.(prop_kz)
    if center_kz
        _, rel_kz = get_McCutchen_kz_center(psz[1:3],pp,sampling)
        prop_pupil .*= cispi(-2*rel_kz/psz[3])
    end
    border_in = 1.0 .- PMLsz[1:2] ./ psz[1:2]
    real_window = exp.(1.75 .* (window_linear(size(start_pupil)[1:2],border_in=border_in,border_out=1) .-1))  # This is maybe not the best PML?
    # real_window = window_gaussian(size(start_pupil)[1:2], border_in=border_in,border_out=border_in.*0.5 .+ 0.5)  # This is maybe not the best PML?
    # real_window = window_hanning(size(start_pupil)[1:2],border_in=border_in,border_out=1)
    prop_pupil = conj(prop_pupil) # from now the advancement is in the opposite direction
    if pp.FFTPlan !== nothing
        P2d = plan_fft(pupil,(1,2), flags=pp.FFTPlan)
    end
    if pp.FFTPlan !== nothing
        Pi2d = plan_ifft(pupil,(1,2), flags=pp.FFTPlan)
    end
    pupil = start_pupil
    for z in start_z:-1:2 # from the middle to the start
        slice = ift2d(pupil) # should be ifft2d for speed reasons
        dst = @view slices[:,:,z:z,:]
        select_region!(slice, dst, new_size=sz[1:2])
        pupil = ft2d(slice .* real_window)
        pupil .*= prop_pupil 
    end
    dst = @view slices[:,:,1:1,:]
    select_region!(ift2d(pupil), dst, new_size=sz[1:2])

    if has_z_symmetry(pp) # to save some speed
        dz = sz[3] - (start_z+1)
        slices[:,:,start_z+1:start_z+1+dz,:] .= conj(slices[:,:,start_z-1:-1:start_z-1-dz,:]);
    else
        prop_pupil = conj(prop_pupil) # from now the advancement is in the opposite direction
        pupil = start_pupil .* prop_pupil
        for z in start_z+1:sz[3]-1  # from the middle forward
            slice = ift2d(pupil)
            dst = @view slices[:,:,z:z,:]
            select_region!(slice, dst, new_size=sz[1:2])
            pupil = ft2d(slice .* real_window)
            pupil .*= prop_pupil 
        end
        dst = @view slices[:,:,sz[3]:sz[3],:]
        select_region!(ift2d(pupil), dst, new_size=sz[1:2])
    end
    return slices
end

function apsf(::Type{MethodPropagateIterative}, sz::NTuple, pp::PSFParams; sampling=nothing, center_kz=false) 
    if isnothing(sampling)
        sampling = get_Ewald_sampling(sz, pp)
        print("Sampling is $sampling \n")
    end
    check_amp_sampling(sz, pp, sampling)
    slices = apply_propagator_iteratively(sz, pp, sampling=sampling, center_kz=center_kz)
    return slices # extract the central bit, which avoids the wrap-around effects
end

# just a different way of writing the propagation down
function apply_kz_to_pupil(pupil, z_planes, pp; sampling)
    sz = (size(pupil)[1:2]...,z_planes,size(pupil)[4])
    if isnothing(sampling)
        sampling = get_Ewald_sampling(sz, pp)
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
function apsf(::Type{MethodShell}, sz::NTuple, pp::PSFParams; sampling=nothing, center_kz=false) 
    if isnothing(sampling)
        sampling = get_Ewald_sampling(sz, pp)
        print("Sampling is $sampling \n")
    end
    check_amp_sampling(sz, pp, sampling)
    pupil = pupil_xyz(sz, pp, sampling) # field_xyz(big_sz,pp, sampling) .* aplanatic_factor(big_sz,pp,sampling) .* ft(jinc_r_2d(big_sz[1:2],pp, sampling=sampling) .* my_disc(big_sz[1:2],pp)) # 

    pupils = apply_kz_to_pupil(pupil, sz[3], pp, sampling=sampling)    
    res=ift2d(pupils) # This should really be a zoomed iFFT
    return res # extract the central bit, which avoids the wrap-around effects
end


# Calculates I0, I1 and I2 according to the Richards and Wolf paper
# The calculation is done on an RZ plane and then interpolated to a 3D volume.
function apsf(::Type{MethodRichardsWolf}, sz::NTuple, pp::PSFParams; sampling=nothing, center_kz=false) 
    @show sampling_r = min(sampling[1],sampling[2])/2.0
    diagonal = sqrt(sum(abs2.(sz[1:2] .* sampling[1:2]))) / 2.0;
    @show sr = ceil(Int64, diagonal/sampling_r) + 1
    r = xx((sr,), scale=sampling_r, offset=CtrCorner) # a generator for radius
    k = 2pi / (pp.λ / pp.n);
    kz = yy((1,sz[3]), scale=k*sampling[3])
    function integrate!(I, w, theta)  # w includes the aplanatic factor
        (sint,cost) = sincos(theta);
        krsint =  (k*sint).*r
        ciskzcost =  cis.(kz.*cost)
        I[:,:,1] .+= (w* sint*(1+cost)).*besselj0.(krsint) .* ciskzcost
        I[:,:,2] .+= (w* sint^2) .*besselj1.(krsint) .* ciskzcost
        I[:,:,3] .+= (w* sint.*(1-cost)).*besselj.(2,krsint) .* ciskzcost
    end
    α = asin(pp.NA / pp.n); # maximal aperture angle
    # Now perform the integration according to Simpsons rule
    N = 80
    I012 = zeros(Complex{pp.dtype}, (sr,sz[3],3)) # I0, I1 and I2
    h = α / N
    for n in 0:N-1 # integrate using Simpsons rule
        theta = h*n
        apl = pp.aplanatic(theta)
        if n==1 || n==N-1
            integrate!(I012, h/6*apl, theta) # first value
        else
            integrate!(I012, 2*h/6*apl, theta) # intermediate values (count twice in the sum)
        end
        theta = h*(n+0.5) # middle value(s)
        integrate!(I012, 4*h/6*apl, theta) # middle values count 4 times
    end
    integrate!(I012, h/6, α*pp.aplanatic(α)) # last value

    # Interpolate the RZ-results onto the 3D grid
    phi = phiphi((sz[1],sz[2]), scale=(sampling[1],sampling[2]))
    rpos = rr((sz[1],sz[2]), scale=(sampling[1],sampling[2]))
    cosphi = cos.(phi)
    cos2phi = cos.(2 .*phi)
    sin2phi = sin.(2 .*phi)
    r_idx = 1 .+ floor.(Int64, rpos ./ sampling_r) # index position
    w = 1.0 .- (rpos/ sampling_r .+ 1 .- r_idx) # index position
    E = zeros(Complex{pp.dtype}, (sz...,3)) # I0, I1 and I2
    for z = 1:sz[3]
        I0 = @view I012[:,z,1]
        I1 = @view I012[:,z,2]
        I2 = @view I012[:,z,3]
        E[:,:,z,1] .= -im.*((w.*I0[r_idx].+(1 .-w).*I0[r_idx.+1]).+cos2phi.*(w.*I2[r_idx].+(1 .-w).*I2[r_idx.+1]))
        E[:,:,z,2] .= -im.*(w.*I2[r_idx].+(1 .-w).*I2[r_idx.+1]).*sin2phi
        E[:,:,z,3] .= -2 .*(w.*I1[r_idx].+(1 .-w).*I1[r_idx.+1]).*cosphi
    end
    return E
end

"""
    apsf(sz::NTuple, pp::PSFParams; sampling=nothing, center_kz=false)

dispatches to various amplitude point spread function calculation routines.
Note that the `method` entry in `pp` defines which calculation method to be used.
Alternatively a different method can be chosen like this: `apsf(::Type{MethodShell}, sz, ..)`.

Arguments:
+ sz: NTuple of size to generate
+ pp: PSF parameter structure, also contains the dtype. See PSFParam() for details
+ sampling: NTuple for pixel pitch information
+ center_kz: if true, the McCutchen pupil will be centered along the kz direction. This is important to be able to apply a consecutive resampling without errors.
             However, the phase values are then not correct, which does not matter for intensity PSFs.

See also:
+ psf():    calculates the intensity point spread function (psf) by taking (sum along field components of) the absolute square of the corresponding apsf. 

Example:
```jdoctest
julia> pp = PSFParams(500.0,1.4,1.52)
julia> p = apsf((128,128,128),pp; sampling=(50,50,100));
```
"""
function apsf(sz::NTuple, pp::PSFParams; sampling=nothing, center_kz=false) 
    apsf(pp.method, sz, pp, sampling=sampling, center_kz=center_kz)
end
