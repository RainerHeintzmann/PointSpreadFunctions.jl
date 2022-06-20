using Test

using PSFs
using NDTools

sz = (128,128,128)
csz = (32,32,32)
sampling=(0.100,0.100,0.150)

function ctr_test(dat1, dat2, rtol=0.01)
    isapprox(select_region(dat1,new_size=csz), select_region(dat2,new_size=csz), rtol=rtol)
end

function compare_asfs(sz, pp, sampling)
    a_prop = apsf(PSFs.MethodPropagate, sz,pp, sampling=sampling);
    a_sincR = apsf(PSFs.MethodSincR, sz,pp, sampling=sampling);
    a_shell = apsf(PSFs.MethodShell, sz,pp, sampling=sampling);
    a_iter = apsf(PSFs.MethodPropagateIterative, sz,pp, sampling=sampling);
    a_RW = apsf(PSFs.MethodRichardsWolf, sz,pp, sampling=sampling);
    sz_big = (512,512,128)
    a_prop2 = NDTools.select_region(apsf(PSFs.MethodPropagate, sz_big,pp, sampling=sampling), new_size=sz);
    @test ctr_test(a_prop, a_iter, 0.05)
    @test ctr_test(a_iter, a_sincR, 0.15)
    @test ctr_test(a_iter, a_shell, 0.1)
    @test ctr_test(a_iter, a_prop2, 0.1)
    @test ctr_test(a_iter, a_RW, 0.1)
    # @vt a_prop2 a_iter a_sincR  a_prop a_shell a_RW
    # mz = size(a_prop2,3)÷2+1; @vt ft2d(a_prop2[:,:,mz:mz,:]) ft2d(a_iter[:,:,mz:mz,:]) ft2d(a_sincR[:,:,mz:mz,:]) ft2d(a_prop[:,:,mz:mz,:]) ft2d(a_shell[:,:,mz:mz,:]) ft2d(a_RW[:,:,mz:mz,:])
    # mz = size(a_prop2,3)÷2+1; @vt ft2d(a_prop2) ft2d(a_iter) ft2d(a_sincR) ft2d(a_prop) ft2d(a_shell) ft2d(a_RW)
end

@testset "Compare various aPSFs" begin
    pp = PSFParams(pol=pol_scalar)
    compare_asfs(sz, pp, sampling)
    pp = PSFParams(pol=pol_circ)
    compare_asfs(sz, pp, sampling)
    pp = PSFParams(pol=pol_x)
    compare_asfs(sz, pp, sampling)
    pp = PSFParams(pol=pol_y)
    compare_asfs(sz, pp, sampling)
end

ct = sz[1].÷2+1
pb = 86
@testset "vectorial pupil" begin
    ppv = PSFParams(pol=pol_x)
    pupil = pupil_xyz(sz, ppv, sampling)
    @test imag(pupil[ct,ct,1,1]) < 1e-8
    @test real(pupil[ct,ct,1,1]) > 1
    @test abs(pupil[ct,ct,1,2]) < 1e-8
    @test abs(pupil[pb,pb,1,2]) > 1
    @test abs(pupil[ct,pb,1,2]) < 1e-6
    @test abs(pupil[pb,ct,1,2]) < 1e-6
    @test abs(pupil[ct,ct,1,3]) < 1e-6
    @test abs(pupil[ct,pb,1,3]) < 1e-6
    @test abs(pupil[pb,ct,1,3]) > 1
end

return
