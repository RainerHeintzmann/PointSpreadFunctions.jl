export pol_scalar, pol_scalar_spiral, pol_x, pol_y, pol_circ, pol_circ_spiral, pol_circ_tophat
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
    pol_circ_tophat(T, xypos, pupil_r)

assumes circular polarization in illumination or an x-oriented polarizer in detection. 
This version includes phase 0/π tophat defining the local (xypos-dependent) phase of both x and y polarization.
"""
function pol_circ_tophat(T, xypos) # e.g. for STED microscopy
    pol_circ(T, xypos) .* cispi.((abs2(xypos[1]) + abs2(xypos[2]) > 1/2))
end



