
"""
    apsf(::Type{MethodParaxial}, sz::NTuple, pp::PSFParams; sampling=nothing) 
    Calculates a paraxial amplitude PSF. Typically `pp.polarization` should be `pol_scalar`. However, other polarisation types yield two channels in the 4th dimension.
    One for X and one for Y polarisation. In the paraxial approximation there is no Z polarisation.
"""
function apsf(::Type{MethodParaxial}, sz::NTuple, pp::PSFParams; sampling=nothing, center_kz=false) 
    # yields a scalar psf
    # res = jinc_r_2d(sz, pp;sampling=sampling) .* field_pupil(sz, pp, sampling)
    error("The paraxial approximation for 3D PSF has not yet been implemented. For a 2D psf please use jinc_r_2d(sz, pp;sampling=sampling) .* field_pupil(sz, pp, sampling)")
    res
end

"""
    apsf(::Type{MethodCZT}, sz::NTuple, pp::PSFParams; sampling=nothing, center_kz=false) 
    Calculate amplitude PSF using Chrip-Z transform. the PSF at a given z_depth requires a bigger window to avoid or reduce wrap-around effect of FFT operation.
    Firstly given an input parameter z_depth to calculate the xy_plane zoom factor c then apply on the pupil based on the desired depth. Calculate the inverse CZT of the 
    field at pupil and zoom-out by the factor c.
    + z_depth: Z-Position from the focus to propagate to. If calcualte 3D apsf (sz[3] is larger than one), a stack is generated with z_depth being the position of the middle (floor.(size/2)) slice.
    + c_want: plane zoom factor only vary by λ, NA, sampling and z_depth. c_allow = λ/(NA*sampling)
"""
function apsf(::Type{MethodCZT}, sz::NTuple, pp::PSFParams; sampling=nothing, center_kz=false) 
    if isnothing(sampling)
        sampling = get_Ewald_sampling(sz, pp)
        print("Sampling is $sampling \n")
    end

    sz, sampling = size_sampling_to3d(sz, sampling)
    check_amp_sampling(sz, pp, sampling)

    wz= sz[1:2] # the window size N_xy.    
    z_depth = 0

    z_planes = (length(sz)>2) ? sz[3] : 1
    se = ifelse(pp.polarization == pol_scalar, 1, 3)
    sz = (wz..., z_planes, se) # initialize a new 4D array of xyz complex numbers and signalton. 
    pupils = Array{Complex{pp.dtype}}(undef, sz)

    mid_point= z_planes ÷2 + 1 
    zoomed_pupil_fwd = Array{Complex{pp.dtype}}(undef, (1,1));
    zoomed_pupil_bwd = Array{Complex{pp.dtype}}(undef, (1,1));
    propagator = Array{Complex{pp.dtype}}(undef, (1,1));
    # for n in 1:sz[3] 
    #     my_dist = n - mid_point + z_depth / sampling[3] 
    #     SliceP, zoomed_pupil = zoomed_pupilProp(sz, pp, zoomed_pupil; sampling, my_dist, center_kz=false) 
    #     pupils[:,:,n:n,:] = SliceP 
    # end
    m=0
    for n in mid_point:-1:1
        my_dist = n - mid_point + z_depth / sampling[3] 
        m = mid_point + (mid_point - n)
        if (m <= size(pupils, 3))
            if (n == mid_point)
                pupils[:,:,m:m,:], _, zoomed_pupil_fwd, zoomed_pupil_bwd, propagator = zoomed_pupilProp(sz, pp, zoomed_pupil_fwd, zoomed_pupil_bwd, propagator, sampling, my_dist; center_kz=false) 
            else
                pupils[:,:,m:m,:], pupils[:,:,n:n,:], zoomed_pupil_fwd, zoomed_pupil_bwd, propagator = zoomed_pupilProp(sz, pp, zoomed_pupil_fwd, zoomed_pupil_bwd, propagator, sampling, my_dist; center_kz=false) 
            end
        else
            _, pupils[:,:,n:n,:], zoomed_pupil_fwd, zoomed_pupil_bwd, propagator = zoomed_pupilProp(sz, pp, zoomed_pupil_fwd, zoomed_pupil_bwd, propagator, sampling, my_dist; center_kz=false) 
        end
    end
    return normalize_amp_to_plane(pupils) 
end


"""
    calculate_CZT_target_size(sz, pp::PSFParams, zoom_factor)

Calculate the target size for the CZT method. The target size is calculated based on
the size of the pupil function in the xy-plane and the desired z-depth. A factor (zoom_factor) determines
how much the pupil function is zoomed in or out in each step.

+ `sz`: the size of the pupil function in the xy-plane.
+ `pp`: the PSFParams object containing the parameters of the PSF.
+ `sampling`: the real-space sampling of the PSF.
+ `zoom_factor`: the factor by which the pupil function is zoomed in or out compared to the central slice. This does not depend on the propagation distance

returned is the target size of the pupil function (`tsz`) and the Float64 zoom factor (`c_want`) to apply for the CZT.
"""
function calculate_CZT_target_size(sz, pp::PSFParams, sampling, my_dist, zoom_factor = 1.3)
    wz = sz[1:2] # original size N_xy: the number of pixels in the xy-focal plane where the pupil fully covers the range and CZT zooms in.
    # calculate the zoom factor c and the desired window size N'_xy.
    hwz = floor.(wz./2) # half image size
    PupilRadius = wz .* (pp.NA / pp.λ) .* sampling[1:2]
    c_allow  = hwz ./ PupilRadius  # allowedZoomFactor = HalfImgSize ./ PupilRadius 

    lambda = pp.λ /pp.n # emission wavelength in the immersion medium in real condition
    # calcualted with parameter z_depth as a heiristic margin to avoid wrap-around. 
    k_xymax = pp.NA / pp.λ # pupil radius in reciprocal space units. 
    kz = sqrt(1/(lambda*lambda)-k_xymax*k_xymax) # calculate kz in reciprocal space units.   
    tan_α = k_xymax / kz
    wzn = 2 .* (hwz .+ tan_α * abs(my_dist*sampling[3]) ./ sampling[1:2]) # D = tan(α) * z_depth, wzn = 2 .* (wz ./2 .+ D)  new window size N'_xy, the extension of this beam is therefore given by D.
    c_want = wzn./ wz .* zoom_factor
    c_want = ifelse.(c_want .< c_allow, c_allow, c_want) # get the zoom "for free", zoom in as much as possible. 
    c_want = floor.(c_want) # round the wanted zoom factor.
    
    if any(c_allow .- c_want .> 0)
        tsz = wz # target size is the new size in pixels need to calculate.
    else
        tsz = ceil.(Int, 2 .* c_want .* PupilRadius) # calculate the new target size such that the pupil just about fits inside
        tsz = tsz .+ (tsz .÷2) # force to be even, due to a problem of czt for uneven sizes.
    end  

    PupilRadiusPro = tsz .* (pp.NA / pp.λ) .* sampling[1:2] 
    # the puple has to be sampled different parameters, so that it fills the availabe number of pixels.
    c_want  = floor.(floor.(tsz./ 2.0) ./ PupilRadiusPro) # wantedZoomFactor = New_HalfImgSize ./ PupilRadiusPro;
    return tsz, c_want
end

"""
    zoomed_pupilProp(sz, pp::PSFParams, zoomed_pupil_fwd, zoomed_pupil_bwd, propagator,  sampling, my_dist,  center_kz=false)

Calculate the forward and backward propagation of the pupil function for a given distance `my_dist` in the z-direction.

+ `sz`: the size of the pupil function in the xy-plane.
+ `pp`: the PSFParams object containing the parameters of the PSF.
+ `zoomed_pupil_fwd`: the forward propagated pupil function.
+ `zoomed_pupil_bwd`: the backward propagated pupil function.
+ `propagator`: the propagator phase for the forward propagation.
+ `sampling`: the sampling of the pupil function.
+ `my_dist`: the distance to propagate the pupil function.
+ `center_kz`: whether to center the pupil function along the kz-axis.

"""
function zoomed_pupilProp(sz, pp::PSFParams, zoomed_pupil_fwd, zoomed_pupil_bwd, propagator, sampling, my_dist; center_kz=false)

    tsz, c_apply = calculate_CZT_target_size(sz, pp, sampling, my_dist)

    # re-compute the zoomed pupil function, only if a recalculation is needed. This is judged solely by the size of the pupil.
    if isempty(zoomed_pupil_fwd) || (size(zoomed_pupil_fwd)[1:2] != tsz)
        zoomed_sampling = sampling .* (c_apply...,1) # broadcast
        zoomed_sz = Int.((tsz...,1)) # because idx(tsz, size::NTuple{N, Int}) 
        zoomed_pupil = pupil_xyz(zoomed_sz, pp, zoomed_sampling) # zoomed_pupil =  field_xyz(Zoomedsz, pp, zoomed_sampling) .* aplanatic_factor(Zoomedsz, pp, zoomed_sampling) .* ft(jinc_r_2d(Zoomedsz, pp, sampling=zoomed_sampling))
        szz =  length(sz)>2  ? sz[3] : 1
        if center_kz
            _, rel_kz = get_McCutchen_kz_center(sz,pp,sampling)
            zoomed_pupil .*= cispi.((-2*rel_kz/szz) .* zz((1,1,szz))) # centers the McCutchen pupil to be able to correctly resample it along kz
        end
        wzn = size(zoomed_pupil)[1:2] # the new window size N'_xy.
        # canvas = (wzn...,szz,size(zoomed_pupil)[4]) # initialize a new 4D array of xyz complex numbers and signalton.
        
        prop_phase, scalar, xy_scale = get_propagator(wzn , pp, zoomed_sampling) # retrieves the propagator phase, propagating a single Z-slice.
        propagator = cis.(prop_phase) # calculated ar returned always for the forward direction
        if (my_dist != 0)
            zoomed_pupil_fwd = zoomed_pupil .* cis.(abs(my_dist).*prop_phase)
            zoomed_pupil_bwd = zoomed_pupil .* cis.(-abs(my_dist).*prop_phase)
        else
            zoomed_pupil_fwd = zoomed_pupil
            zoomed_pupil_bwd = zoomed_pupil
        end
    else
        zoomed_pupil_fwd = zoomed_pupil_fwd .* propagator
        zoomed_pupil_bwd = zoomed_pupil_bwd .* conj.(propagator)
    end

    c_4dim = (c_apply..., 1, 1)  # expand the scale size to fit scale argument of iczt(). 
    # 2D inverse CZT for each slice:
    wz = sz[1:2] # original size N_xy: the number of pixels in the xy-focal plane where the pupil fully covers the range and CZT zooms in.
    # has_z-symmetry does NOT work properly for vectorial situations.
    if false # (has_z_symmetry(pp))
        propagated_fwd = iczt(zoomed_pupil_fwd, c_4dim, 1:2, wz) 
        propagated_bwd = conj.(propagated_fwd)
    else
        propagated_fwd = iczt(zoomed_pupil_fwd, c_4dim, 1:2, wz) 
        propagated_bwd = iczt(zoomed_pupil_bwd, c_4dim, 1:2, wz) 
    end

    return propagated_fwd, propagated_bwd, zoomed_pupil_fwd, zoomed_pupil_bwd, propagator
end  

function apsf(::Type{MethodSincR}, sz::NTuple, pp::PSFParams; sampling=nothing, center_kz=false) 
    sz = (length(sz)>2) ? sz : (sz[1:2]..., 1) # to also work for 2D input sizes
    check_amp_sampling(sz, pp, sampling)
    if sz[3] < 2
        error("Method SincR only makes sense for calculating 3D psfs with non-signalton z-size.")
    end

    if isnothing(sampling)
        sampling = get_required_amp_sampling(sz, pp) # this is knowingly too small along kz but fixed down below, if needed.
    end
    Ewald_sampling = get_Ewald_sampling(sz, pp)
    undersampling_factor = sampling[3] ./ Ewald_sampling[3]
    # The big_sz is used to avoid interference effects by the application of the multiplication with the pupil maks in Fourier space
    # The size if 4/3 is such that the diagonals are covered but a wrap-around into the original size is avoided.
    # big_sz = ((sz[1:2].*2)...,sz[3]+2)  # The 2 extra slices in z seem to be necessary to avoid a problem at the fist position
    big_sz = (ceil.(Int, sz[1:2].*2.1)...,sz[3]+4)  # The 2 extra slices in z seem to be necessary to avoid a problem at the fist positions

    if undersampling_factor > 1.0 # the current sampling is insufficient for the SincR method. We need to properly sample the Ewald sphere.
        # the size sz[3] will be cut out in Fourier-spage around the cleaned and processed part of the McCutchen pupil, which should yield the sampling as currently specified.
        nowrap_sz = (big_sz[1:2]..., ceil(Int, undersampling_factor .* big_sz[3]))
        big_sampling = (sampling[1:2]..., sampling[3] * big_sz[3] / nowrap_sz[3]) # calculate the sampling needed to cover the Ewald sphere.
    else
        nowrap_sz = big_sz
        big_sampling = sampling
    end

    # the pupil below does not need the 1/cos(theta) factor, since this is already in the 3D shell.
    pupil = pupil_xyz(nowrap_sz, pp, big_sampling, is_proj=false) # field_xyz(big_sz,pp, sampling) .* aplanatic_factor(big_sz,pp,sampling) .* ft(jinc_r_2d(big_sz[1:2],pp, sampling=sampling) .* my_disc(big_sz[1:2],pp)) # 
    
    # relatively expensive
    sinc_r_big = sinc_r(nowrap_sz,pp, sampling=big_sampling) .* my_disc(nowrap_sz[1:2],pp)  # maybe this should rather already be apodized by angle?

    if !isnothing(pp.FFTPlan) 
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
    # pupils = (cos.(pupil_θ(nowrap_sz,pp,sampling)) .* pupil) .* iftz(shell)
    pupils = pupil .* iftz(shell)
    if !isnothing(pp.FFTPlan) 
        Pm2d = plan_fft(pupils,(1,2), flags=pp.FFTPlan)
    end

    res=ift2d(pupils) # This should really be a zoomed iFFT

    # ToDo: this normalization can probably be avoided if the pupil is normalized correctly.
    return normalize_amp_to_plane(select_region(res, new_size=sz[1:3])) # extract the central bit, which avoids the wrap-around effects
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
    sz, sampling = size_sampling_to3d(sz, sampling)

    check_amp_sampling(sz, pp, sampling)
    # the pupil below is the ft of a jinc in real space and includes a factor of my_disc(sz[1:2],pp) to reduce wrap around
    pupil = pupil_xyz(sz, pp, sampling) # field_xyz(big_sz,pp, sampling) .* aplanatic_factor(big_sz,pp,sampling) .* ft(jinc_r_2d(big_sz[1:2],pp, sampling=sampling) .* my_disc(big_sz[1:2],pp)) 
    szz = (length(sz)>2) ? sz[3] : 1
    pupils = apply_propagators(pupil, szz, pp, sampling=sampling)

    # return pupils
    # pupils .*= window_hanning((1,1,size(pupils,3)),border_in=0.8,border_out=1.0,dims=(3,))
    if center_kz
        _, rel_kz = get_McCutchen_kz_center(sz,pp,sampling)
        pupils .*= cispi.((-2*rel_kz/szz) .* zz((1,1,szz))) # centers the McCutchen pupil to be able to correctly resample it along kz
    end
    if !isnothing(pp.FFTPlan)
        Pm2d = plan_fft(pupils,(1,2), flags=pp.FFTPlan)
    end
    res = ift2d(pupils) # This should really be a zoomed iFFT
     # ToDo: this normalization can probably be avoided if the pupil is normalized correctly.
     return normalize_amp_to_plane(res) # extract the central bit, which avoids the wrap-around effects
end

# This function is needed for the propagate method, but iterates between real and Fourier space always smoothly deleting "out-of-bound" waves.
# The final result is already in real space.
function apply_propagator_iteratively(sz, pp::PSFParams; sampling, center_kz=false) 
    z_planes = (length(sz)>2) ? sz[3] : 1
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
    border_in = pp.dtype.(1.0 .- PMLsz[1:2] ./ psz[1:2])
    real_window = collect(exp.(pp.dtype(1.75) .* (window_linear(pp.dtype, size(start_pupil)[1:2],border_in=border_in,border_out=1) .-1)))  # This is maybe not the best PML?
    # real_window = window_gaussian(size(start_pupil)[1:2], border_in=border_in,border_out=border_in.*0.5 .+ 0.5)  # This is maybe not the best PML?
    # real_window = window_hanning(size(start_pupil)[1:2],border_in=border_in,border_out=1)
    prop_pupil = conj(prop_pupil) # from now the advancement is in the opposite direction
    if !isnothing(pp.FFTPlan) 
        P2d = plan_fft(pupil,(1,2), flags=pp.FFTPlan)
    end
    if !isnothing(pp.FFTPlan)
        Pi2d = plan_ifft(pupil,(1,2), flags=pp.FFTPlan)
    end
    pupil = start_pupil
    for z = start_z:-1:2 # from the middle to the start
        slice = collect(ift2d(pupil)) # should be ifft2d for speed reasons
        # dst = @view slices[:,:,z:z,:]
        # writes the slice into the destination arry in slices
        # select_region!(slice, dst, new_size=sz[1:2])
        slices[:,:,z:z,:] .= select_region_view(slice,  new_size=sz[1:2])
        pupil = ft2d(slice .* real_window)
        pupil .*= prop_pupil 
    end
    slices[:,:,1:1,:] .= select_region_view(collect(ift2d(pupil)),  new_size=sz[1:2])
    # dst = @view slices[:,:,1:1,:]
    # select_region!(ift2d(pupil), dst, new_size=sz[1:2])

    # Does it need also an XY-flip to be correct?
    if false # has_z_symmetry(pp) # to save some speed
        dz = sz[3] - (start_z+1)
        # missing XY-flip:
        slices[:,:,start_z+1:start_z+1+dz,:] .= conj.(slices[:,:,start_z-1:-1:start_z-1-dz,:]);
    else
        prop_pupil = conj(prop_pupil) # from now the advancement is in the opposite direction
        pupil = start_pupil .* prop_pupil
        for z in start_z+1:sz[3]  # from the middle forward
            slice = collect(ift2d(pupil))
            slices[:,:,z:z,:] .= select_region_view(slice,  new_size=sz[1:2])
            if z < sz[3]
                pupil = ft2d(slice .* real_window)
                pupil .*= prop_pupil 
            end
        end
    end        
    return slices
end

function apsf(::Type{MethodPropagateIterative}, sz::NTuple, pp::PSFParams; sampling=nothing, center_kz=false) 
    if isnothing(sampling)
        sampling = get_Ewald_sampling(sz, pp)
        print("Sampling is $sampling \n")
    end

    sz, sampling = size_sampling_to3d(sz, sampling)
    check_amp_sampling(sz, pp, sampling)
     # ToDo: this normalization can probably be avoided if the pupil is normalized correctly.
    return normalize_amp_to_plane(apply_propagator_iteratively(sz, pp, sampling=sampling, center_kz=center_kz))
end

# just a different way of writing the propagation down
function apply_kz_to_pupil(pupil, z_planes, pp; sampling)
    sz = (size(pupil)[1:2]...,z_planes,size(pupil)[4])
    if isnothing(sampling)
        sampling = get_Ewald_sampling(sz, pp)
    end
    k_max_rel = sampling[1:2] ./ (pp.λ / pp.n)
    scalar = pp.dtype((2π*sampling[3] / (pp.λ / pp.n)))
    xy_scale = 1 ./ (k_max_rel .* sz[1:2])
    sqrt_term = phase_kz(pp.dtype, sz[1:2], scale = xy_scale)
    prop_phase = scalar .* sqrt_term    
    pupils = pupil .* cis.(zz(pp.dtype,(1,1,z_planes)) .* prop_phase) # if you fourier-transform this along kz you should get the shifted pupil
    return pupils
end


# Basic idea of MethodShell: Displace each kxy by the appropriate kz 
function apsf(::Type{MethodShell}, sz::NTuple, pp::PSFParams; sampling=nothing, center_kz=false) 
    if isnothing(sampling)
        sampling = get_Ewald_sampling(sz, pp)
        print("Sampling is $sampling \n")
    end
    sz, sampling = size_sampling_to3d(sz, sampling)
    check_amp_sampling(sz, pp, sampling)
    pupil = pupil_xyz(sz, pp, sampling) # field_xyz(big_sz,pp, sampling) .* aplanatic_factor(big_sz,pp,sampling) .* ft(jinc_r_2d(big_sz[1:2],pp, sampling=sampling) .* my_disc(big_sz[1:2],pp)) # 

    pupils = apply_kz_to_pupil(pupil, sz[3], pp, sampling=sampling)
    if center_kz
        _, rel_kz = get_McCutchen_kz_center(sz,pp,sampling)
        pupils .*= cispi.((-2*rel_kz/sz[3]) .* zz((1,1,sz[3]))) # centers the McCutchen pupil to be able to correctly resample it along kz
    end
    if !isnothing(pp.FFTPlan)
        Pm2d = plan_fft(pupils,(1,2), flags=pp.FFTPlan)
    end

    res=ift2d(pupils) # This should really be a zoomed iFFT
     # ToDo: this normalization can probably be avoided if the pupil is normalized correctly.
    return normalize_amp_to_plane(res) # extract the central bit, which avoids the wrap-around effects
end

function simpson!(f,α,N=80)
    h = α / N # N+1 points are really used...
    I012 = h/6*f(0) # middle values count 4 times
    for n in 0:N-1 # integrate using Simpsons rule
        theta = h*(n+0.5) # middle value(s)
        I012 += 4*h/6*f(theta) # middle values count 4 times
        if n<(N-1)
            theta = h*(n+1)
            I012 += 2*h/6*f(theta) # intermediate values (count twice in the sum)
        end
    end
    I012 += h/6*f(α) # last value
    return I012
end

# Calculates I0, I1 and I2 according to the Richards and Wolf paper
# The calculation is done on an RZ plane and then interpolated to a 3D volume.
# This is probably still buggy! At least the Pupil aplanatic factor looks wrong!
function apsf(::Type{MethodRichardsWolf}, sz::NTuple, pp::PSFParams; sampling=nothing, center_kz=false) 
    if length(pp.aberrations.indices) > 0
        error("The Richards & Wolf amplitude spread function calculations does currently not support aberrations. Please choose a different method.")
    end
    sz, sampling = size_sampling_to3d(sz, sampling)

    sampling_r = min(sampling[1],sampling[2])/3
    diagonal = sqrt(sum(abs2.(sz[1:2] .* sampling[1:2]))) / 2.0;
    sr = ceil(Int64, diagonal/sampling_r) + 1
    # @show sr
    r = xx((sr,), scale=sampling_r, offset=CtrCorner) # a generator for radius
    k = 2pi / (pp.λ / pp.n);
    szz = sz[3] 
    kz = yy((1, szz), scale=k*sampling[3])
    function integrate!(I, w, theta)  # w includes the aplanatic factor
        (sint,cost) = @fastmath sincos(theta);
        krsint =  (k*sint).*r
        ciskzcost =  @fastmath cis.(kz.*cost)
        I[:,:,1] .+= (w* sint*(1+cost)).*(@fastmath besselj0.(krsint)) .* ciskzcost
        I[:,:,2] .+= (w* sint^2) .*(@fastmath besselj1.(krsint)) .* ciskzcost
        I[:,:,3] .+= (w* sint.*(1-cost)).*(@fastmath besselj.(2,krsint)) .* ciskzcost
    end
    α = asin(pp.NA / pp.n); # maximal aperture angle
    # Now perform the integration according to Simpsons rule
    N = 150 # 50 does not seem to be good enough
    h = α / N # since these are really N+1 points
    I012 = zeros(Complex{pp.dtype}, (sr, szz, 3)) # I0, I1 and I2

    # I think this code is wrong and according to C. Sheppard this should be simply
    # aplanatic_fct = pp.aplanatic
    # all the other functions need to account for 1/cos(theta) due to the projection of the k-sphere. 
    aplanatic_fct = let
        if pp.polarization == pol_scalar
            # pp.aplanatic
            (θ) -> pp.aplanatic(θ) .* sqrt.(max.(cos.(θ),zero(typeof(θ))))
        else
            # pp.aplanatic
            (θ) -> pp.aplanatic(θ) .* max.(cos.(θ),zero(typeof(θ)))
        end
    end
    integrate!(I012, h/6*aplanatic_fct(0.0), 0.0) # first value
    for n in 0:N-1 # integrate using Simpsons rule
        theta = h*(n+0.5) # middle value(s)
        integrate!(I012, 4*h/6*aplanatic_fct(theta), theta) # middle values count 4 times
        if n<(N-1)
            theta = h*(n+1)
            integrate!(I012, 2*h/6*aplanatic_fct(theta), theta) # intermediate values (count twice in the sum)
        end
    end

    integrate!(I012, h/6*aplanatic_fct(α), α) # last value

    # Interpolate the RZ-results onto the 3D grid
    phi = phiphi((sz[1],sz[2]), scale=(sampling[1],sampling[2]))
    sinphi = @fastmath sin.(phi)
    cosphi = @fastmath cos.(phi)
    cos2phi = @fastmath cos.(2 .*phi)
    sin2phi = @fastmath sin.(2 .*phi)
    rpos = rr((sz[1],sz[2]), scale=(sampling[1],sampling[2]))
    r_idx = 1 .+ floor.(Int64, rpos ./ sampling_r) # index position
    # @show findall(r_idx.<4)
    w = 1.0 .- (rpos./sampling_r .+ 1 .- r_idx) # index position
    numEl = let 
        if pp.polarization == pol_scalar
            1
        else
            3
        end
    end
    E = zeros(Complex{pp.dtype}, (sz..., numEl)) # I0, I1 and I2
    _, rel_kz = get_McCutchen_kz_center((sz[1:2]..., szz),pp,sampling)

    interp_lin(vals) = (w.*(@view vals[r_idx]).+(1 .-w).*(@view vals[r_idx.+1]))

    # Langrange's quadratic interpolation
    # coefficients 
    L0 = (w.-1).*w./2  # (x - x1)(x -x2) / 2  # x0 = -1, x1 = 0, x2 = 1, w = (1 -x)  --> x = (1 - w)
    L1 = (2 .-w).*w./1  # -(x - x0)(x -x2) # 
    L2 = (2 .-w).*(1 .-w)./2 # (x-x0)(x-x1)/2 
    interp_sqr(vals) = L0 .*(@view vals[abs.(r_idx.-2).+1]) .+
                       L1 .*(@view vals[r_idx]) .+
                       L2 .*(@view vals[r_idx.+1])
    # some code to test the quadratic interpolation
    # test_fun(p) = (p + 17.321)^2; # quadratic function to test with
    # test_vals = test_fun.(rpos ./ sampling_r) # true values
    # test_table = test_fun.(0:1000); # a table to interpolate from
    # @show maximum(abs.(test_vals .- interp_sqr(test_table))) # smaller than 10^-9
    # return rpos ./ sampling_r, test_vals, interp_sqr(test_table) 
    # code to test: t1,t2 = apsf(MethodRichardsWolf, (256,256,1), pp, sampling=samp)
    # @vt t1 t2

    interp = interp_sqr;

    for z = 1:szz
        I0 = @view I012[:,z,1] # create views to be able to write by 1D indexing into the appropriate slices.
        I1 = @view I012[:,z,2]
        I2 = @view I012[:,z,3]

        z_pos = z - (szz÷2+1)

        rel_phase = center_kz ? Complex{pp.dtype}(@fastmath cispi.((-2*rel_kz/szz) .* z_pos)) : one(pp.dtype) # centers the McCutchen pupil to be able to correctly resample it along kz
    
        if pp.polarization == pol_x
            E[:,:,z,1] .= rel_phase .* (interp(I0).+cos2phi.*interp(I2))
            E[:,:,z,2] .= rel_phase .* interp(I2).*sin2phi
            E[:,:,z,3] .= rel_phase .* (2im .* interp(I1).*cosphi)
        elseif pp.polarization == pol_y
            E[:,:,z,1] .= rel_phase .* interp(I2).*sin2phi # sin2phi with pi/2 phase shift becomes -sin2phi
            E[:,:,z,2] .= rel_phase .* (interp(I0).-cos2phi.*interp(I2)) # cos2phi becomes -cos2phi
            E[:,:,z,3] .= rel_phase .* 2im .*interp(I1).*sinphi # cosphi becomes -sinphi 
        elseif pp.polarization == pol_circ # pol_x + im* pol_y
            E[:,:,z,1] .= rel_phase .* (1/sqrt(2).*(interp(I0).+cos2phi.*interp(I2)) .+
            im/sqrt(2)*interp(I2).*sin2phi)
            E[:,:,z,2] .= rel_phase .* (1/sqrt(2).*interp(I2).*sin2phi .+
            im/sqrt(2)*(interp(I0).-cos2phi.*interp(I2)))
            E[:,:,z,3] .= rel_phase .* (2im/sqrt(2)*interp(I1).*(cosphi .+ im * sinphi))
        elseif pp.polarization == pol_scalar # scalar approximation
            E[:,:,z,1] .= rel_phase .* interp(I0) # .+ (w.*I2[r_idx].+(1 .-w).*I2[r_idx.+1]) .+ 1 .*(w.*I1[r_idx].+(1 .-w).*I1[r_idx.+1])
        else
            error("unsupported polarization for Richards-Wolf method")
        end
    end
     # ToDo: this normalization can probably be avoided if the pupil is normalized correctly.
     return normalize_amp_to_plane(E) # no idea why this scaling is needed: .* 1.14 
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

