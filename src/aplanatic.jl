export aplanatic_detection, aplanatic_illumination, aplanatic_const, aplanatic_illumination_flux

"""
    aplanatic_detection = (θ) -> sqrt.(max.(0,cos.(θ)))
This is the aplanatic factor typically used in detection of fluorescence of (randomly oriented) fluorophores. 
"""
aplanatic_detection = (θ) -> sqrt.(max.(0,cos.(θ)))  # returns only one "polarization" indicating that the whole calculation is to be performed scalar

"""
    aplanatic_illumination = (θ) -> sqrt.(max.(0,cos.(θ)))
This is the aplanatic factor typically used in illumination of (randomly oriented) fluorophores. Note that it is identical to detection.
"""
aplanatic_illumination = (θ) -> sqrt.(max.(0,cos.(θ)))  # returns only one "polarization" indicating that the whole calculation is to be performed scalar

"""
    aplanatic_const = (θ) -> one.(eltype(θ))
This is a constant aplanatic factor
"""
aplanatic_const = (θ) -> one.(eltype(θ))  # returns only one "polarization" indicating that the whole calculation is to be performed scalar

"""
    aplanatic_illumination_flux = (θ) -> max.(0,cos.(θ))
This refers the aplanatic factor if interested in the flux through a detector perpendicular to the optical axis.
"""
aplanatic_illumination_flux = (θ) -> max.(0,cos.(θ)) # flux through a voxel

