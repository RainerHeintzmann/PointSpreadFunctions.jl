export pol_scalar, pol_scalar_spiral, pol_x, pol_y, pol_circ, pol_circ_spiral, pol_circ_tophat, pol_circ_quadrant, pol_radial,  pol_radial_annulus
export pupil_apodize_hann, pupil_apodize_cos, pupil_annulus
# define a bunch of standard pupil polarizations

"""
    pol_scalar(T, xypos)

ignores polarization aspects in the calculation but calculates only based on (high-NA) scalar theory.
This is a lot faster but not as accurate.

"""
function pol_scalar(T, xypos)  # returns only one "polarization" indicating that the whole calculation is to be performed scalar
    (one(T),)
end

"""
    pol_scalar_spiral(T, xypos)

ignores polarization aspects in the calculation but calculates only based on (high-NA) scalar theory.
This version still includes a (scalar) phase spiral.
This is a lot faster but not as accurate.
"""
function pol_scalar_spiral(T, xypos) # e.g. for STED microscopy
    cis.(atan(T.(xypos...)))
end

"""
    pol_x(T, xypos)

assumes x-polarization in illumination or an x-oriented polarizer in detection. 
In a high-NA objective this is converted into XYZ electric fields at the focus.
"""
function pol_x(T, xypos)
    (one(T), zero(T))
end

"""
    pol_y(T, xypos)

assumes y-polarization in illumination or an x-oriented polarizer in detection. 
In a high-NA objective this is converted into XYZ electric fields at the focus.
"""
function pol_y(T, xypos)
    (zero(T), one(T))
end

"""
    pol_circ(T, xypos)

assumes circular polarization in illumination or an x-oriented polarizer in detection. 
In a high-NA objective this is converted into XYZ electric fields at the focus.
"""
function pol_circ(T, xypos)
    (one(T)/sqrt(T(2))+T(0)im, one(T)/sqrt(T(2))*T(1)im)
end

"""
    pol_circ_spiral(T, xypos)

assumes circular polarization in illumination or an x-oriented polarizer in detection. 
This version includes phase spiral defining the local (xypos-dependent) phase of both x and y polarization.
"""
function pol_circ_spiral(T, xypos) # e.g. for STED microscopy
    pol_circ(T, xypos) .* cis.(atan(T(xypos[2]),T(xypos[1])))
end

"""
    pol_circ_tophat(T, xypos)

assumes circular polarization in illumination. 
This version includes phase 0/π tophat defining the local (xypos-dependent) phase of both x and y polarization.
"""
function pol_circ_tophat(T, xypos) # e.g. for STED microscopy
    pol_circ(T, xypos) .* cispi.((abs2(xypos[1]) + abs2(xypos[2]) > 1/2))
end

"""
    pol_circ_quadrant(T, xypos)

assumes circular polarization in illumination. 
This version implements quadrant phase steps, which is one version of STED phaseramps.
"""
function pol_circ_quadrant(T, xypos) # e.g. for STED microscopy
    pol_circ(T, xypos) .*  cis.(T(pi)/2*floor(atan(T(xypos[2]),T(xypos[1])) / (T(pi)/2)))
end



"""
    pol_radial(T, xypos)

assumes radial polarization in (illumination/detection) of at the pupil.

Example:
```julia
    pp_em = PSFParams(0.532, 1.3, 1.52; mode=ModeWidefield, pol=pol_radial);
    h_p = apsf(MethodPropagate, sz, pp_em, sampling=samp);
    @vv real.(h_p[:,:,1,3])
``
"""
function pol_radial(T, xypos) # e.g. for STED microscopy
    xypos = T.(xypos)
    xypos ./ sqrt.(1e-10+abs2(xypos[1])+abs2(xypos[2]))
end

# """
#     combine_pupils(T, xypos; fct1, fct2)

# combines two pupil functions fct1 and fct2 by multiplication. 
# """
# function combine_pupils(fct1, fct2)
#     return (T, xypos) -> fct1(T, xypos) .* fct2(T, xypos)
# end

"""
    pupil_annulus(T, xypos; r0= 0.8, σ=0.05) 

returns a pupil function with an annulus at relative pupil radius r0 of exponential apodization width σ.
Note that this pupil function should not be used by itself, but needs to be combined with a polarization function.
"""
function pupil_annulus(T, xypos; r0= 0.8, σ=0.05) # e.g. for STED microscopy
    xypos = T.(xypos)
    exp(-abs2(sqrt(abs2(xypos[1])+abs2(xypos[2]))-r0)/(2*abs2(σ))) 
end
    
"""
    pupil_apodize_hann(T, xypos; r0= 0.9) # 

returns a pupil function with a hanning-type apodization starting at radius r0, extending to the pupil limit.
Both, cos-window and Hann window apodization should yield a finite standard deviation of the PSF.
Note that this pupil function should not be used by itself, but needs to be combined with a polarization function.
"""
function pupil_apodize_hann(T, xypos; r0= 0.9) # e.g. for STED microscopy
    xypos = T.(xypos)
    r = max(0, sqrt(abs2(xypos[1])+abs2(xypos[2]))-r0)/(1-r0)
    return (1+cos(pi*r))/2
end
    
"""
    pupil_apodize_cosann(T, xypos; r0= 0.9) # 

returns a pupil function with a cos-type apodization starting at radius r0, extending to the pupil limit.
The cos-type apodization reaches zero at the border of the pupil, but with a slope.
This cos-filter goes back to Gabor. Both, cos-window and Hann window apodization should yield a finite standard deviation of the PSF.
Note that this pupil function should not be used by itself, but needs to be combined with a polarization function.
"""
function pupil_apodize_cos(T, xypos; r0= 0.9) # e.g. for STED microscopy
    xypos = T.(xypos)
    r = max(0, sqrt(abs2(xypos[1])+abs2(xypos[2]))-r0)/(1-r0)
    return cos(pi*r/2)
end

function pol_radial_annulus(T, xypos; r0= 0.8, σ=0.05) # e.g. for STED microscopy
    return pol_radial(T, xypos) .* pupil_annulus(T, xypos; r0=r0, σ=σ)
end

