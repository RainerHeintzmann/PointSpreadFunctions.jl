using Test
using Random
using NDTools
using BenchmarkTools
using FFTW

using PSFs

pp = PSFParams(580.0, 1.4, 1.518; pol=pol_scalar,method=PSFs.MethodSincR, FFTPlan=nothing)
sz = (256,256,64)
sampling=(80,80,100)

@benchmark p = psf(sz,$pp, sampling=sampling,use_resampling=true)
# (with extra border): 397 ms using standard heuristics
# Old Data (without extra border): 303 ms -> 0.145 µs / voxel  =  0.304 / prod(sz) * 1e6

aberr = PSFs.Aberrations([Zernike_Spherical],[0.1])
ppm = PSFParams(580.0, 1.4, 1.518;pol=pol_scalar,method=PSFs.MethodSincR, aberrations= aberr, FFTPlan=FFTW.MEASURE)
p = psf(sz, ppm, sampling=sampling,use_resampling=true); # just to ensure that the measure is already executed

@benchmark p = psf(sz,$pp, sampling=sampling,use_resampling=true)
# (with extra border): 371 ms  # using MEASURE
# Old Data (without extra border): 267 ms -> 0.127 µs/voxel = 0.267343 / prod(sz) * 1e6

