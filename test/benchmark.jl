using Test
using Random

using PSFs
using NDTools
using BenchmarkTools
using FFTW

pp = PSFParams(580.0, 1.4, 1.518; pol=pol_scalar,method=PSFs.MethodSincR, FFTPlan=nothing)
sz = (256,256,32)
sampling=(80,80,100)

@benchmark p = psf(sz,$pp, sampling=sampling,use_resampling=true)
# 303 ms -> 0.145 µs / voxel  =  0.304 / prod(sz) * 1e6

ppm = PSFParams(580.0, 1.4, 1.518;pol=pol_scalar,method=PSFs.MethodSincR, FFTPlan=FFTW.MEASURE)
p = psf(sz, ppm, sampling=sampling,use_resampling=true);

@benchmark p = psf(sz,$pp, sampling=sampling,use_resampling=true)
# 267 ms -> 0.127 µs/voxel = 0.267343 / prod(sz) * 1e6

