using Test
using Random

using PSFs
using BenchmarkTools

pp = PSFParams()
sz = (128,128,128)
sampling=(100,100,150)

@testset "Compare scalar aPSFs" begin
    a_prop = apsf(PSFs.MethodPropagate, sz,pp, sampling=sampling);
    a_sincR = apsf(PSFs.MethodSincR, sz,pp, sampling=sampling);
    a_shell = apsf(PSFs.MethodShell, sz,pp, sampling=sampling);
    # @vt a_prop a_sincR
end

@testset "vectorial pupil" begin
    ppv = PSFParams(pol=pol_x)
    pupil = pupil_xyz(sz, ppv, sampling)
end

return
