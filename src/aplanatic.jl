export aplanatic_detection, aplanatic_illumination, aplanatic_const, aplanatic_illumination_flux

aplanatic_detection = (θ) -> sqrt.(max.(0,cos.(θ)))  # returns only one "polarization" indicating that the whole calculation is to be performed scalar
aplanatic_illumination = (θ) -> sqrt.(max.(0,cos.(θ)))  # returns only one "polarization" indicating that the whole calculation is to be performed scalar
aplanatic_const = (θ) -> one.(eltype(θ))  # returns only one "polarization" indicating that the whole calculation is to be performed scalar
aplanatic_illumination_flux = (θ) -> max.(0,cos.(θ)) # flux through a voxel

