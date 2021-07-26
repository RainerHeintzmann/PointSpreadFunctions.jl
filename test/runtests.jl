using Test
using Random

using PSFs
using NDTools
using BenchmarkTools

pp = PSFParams()
sz = (128,128,128)
sampling=(100,100,150)

@testset "Compare scalar aPSFs" begin
    a_prop = apsf(PSFs.MethodPropagate, sz,pp, sampling=sampling);
    a_sincR = apsf(PSFs.MethodSincR, sz,pp, sampling=sampling);
    a_shell = apsf(PSFs.MethodShell, sz,pp, sampling=sampling);
    a_iter = apsf(PSFs.MethodPropagateIterative, sz,pp, sampling=sampling);
    sz_big = (512,512,128)
    a_prop2 = NDTools.select_region(apsf(PSFs.MethodPropagate, sz_big,pp, sampling=sampling), new_size=sz);
    # @vt a_prop2 a_iter a_sincR  a_prop a_shell
    # @vtp ft2d(a_prop2) ft2d(a_iter) ft2d(a_sincR) ft2d(a_prop) ft2d(a_shell)  
end

pp = PSFParams(pol=pol_x)
sz = (128,128,128)
sampling=(100,100,150)

@testset "Compare scalar aPSFs" begin
    a_prop = apsf(PSFs.MethodPropagate, sz,pp, sampling=sampling);
    a_sincR = apsf(PSFs.MethodSincR, sz,pp, sampling=sampling);
    a_shell = apsf(PSFs.MethodShell, sz,pp, sampling=sampling);
    a_iter = apsf(PSFs.MethodPropagateIterative, sz,pp, sampling=sampling);
    # @vt a_prop a_sincR a_iter
end

@testset "vectorial pupil" begin
    ppv = PSFParams(pol=pol_x)
    pupil = pupil_xyz(sz, ppv, sampling)
end

return
