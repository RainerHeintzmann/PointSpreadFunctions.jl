using Test

using PointSpreadFunctions
using NDTools

sz = (128,128,128)
csz = (32,32,32)
sampling=(0.100,0.100,0.150)

function ctr_test(dat1, dat2, rtol=0.01)
    isapprox(select_region(dat1,new_size=csz), select_region(dat2,new_size=csz), rtol=rtol)
end

function compare_asfs(sz, pp, sampling)
    @time a_prop = apsf(PointSpreadFunctions.MethodPropagate, sz, pp, sampling=sampling);
    @time a_sincR = apsf(PointSpreadFunctions.MethodSincR, sz ,pp, sampling=sampling);
    @time a_shell = apsf(PointSpreadFunctions.MethodShell, sz, pp, sampling=sampling);  # fast
    @time a_iter = apsf(PointSpreadFunctions.MethodPropagateIterative, sz ,pp, sampling=sampling);
    @time a_RW = apsf(PointSpreadFunctions.MethodRichardsWolf, sz, pp, sampling=sampling);
    sz_big = (512,512,128)
    @time a_prop2 = NDTools.select_region(apsf(PointSpreadFunctions.MethodPropagate, sz_big,pp, sampling=sampling), new_size=sz);
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
    pp = PSFParams(pol=pol_scalar, aplanatic=aplanatic_illumination)
    compare_asfs(sz, pp, sampling)
    pp = PSFParams(pol=pol_scalar, aplanatic=aplanatic_detection)
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

@testset "confocal PSF" begin
    sampling = (0.04,0.04,0.120)
    sz = (128,128,128)
    pp_em = PSFParams(0.5,1.3,1.52; mode=ModeConfocal);
    pp_ex = PSFParams(pp_em; λ=0.488, aplanatic=aplanatic_illumination);
    pinhole = 0.001
    @time pc = psf(sz,pp_em; pp_ex=pp_ex, pinhole=pinhole, sampling=sampling);
    @time pc_open = psf(sz,pp_em; pp_ex=pp_ex, pinhole=5.0, sampling=sampling);
    @time pc_open2 = psf(sz,pp_em; pp_ex=pp_ex, pinhole=5.0, sampling=sampling, use_resampling=false);

    pp_em = PSFParams(0.5,1.3,1.52; mode=ModeWidefield);
    pw_em = psf(sz,pp_em; sampling=sampling);
    pp_ex = PSFParams(pp_em; λ=0.488, aplanatic=aplanatic_illumination);
    pw_ex = psf(sz,pp_ex; sampling=sampling);
    pc2 = pw_ex .* pw_em
    # check if the confocal calculations agree for a very small pinhole
    @test isapprox(pc ./ sum(pc), pc2./sum(pc2); rtol=0.05)
    @test isapprox(pc ./ maximum(pc), pc2./maximum(pc2); rtol=0.01)
    # test if the pinhole normalization is correct. For a very large pinhole the PSF is an excitation WF PSF
    @test ctr_test(pc_open2, pw_ex, 0.15)
    @test ctr_test(pc_open, pw_ex, 0.15)
end

@testset "lightsheet PSF" begin
    sz = (64,64,64)
    pp_em = PSFParams(0.5,0.9,1.33; mode=ModeLightsheet);
    pp_ex = PSFParams(0.488, 0.2, 1.33);
    p_lightsheet1 = psf(sz, pp_em; pp_ex=pp_ex, sampling=(0.040,0.040,0.100), dynamically_scanned=true);
    p_lightsheet2 = psf(sz, pp_em; pp_ex=pp_ex, sampling=(0.040,0.040,0.100), dynamically_scanned=false);
    p_lightsheetg = psf(sz, pp_em; sampling=(0.040,0.040,0.100), sigma_z=1.0);

    @test size(p_lightsheet1) == sz
    @test size(p_lightsheet2) == sz
    @test size(p_lightsheetg) == sz
    @test p_lightsheet1[33,33,33] > p_lightsheet2[33,33,33]
    @test p_lightsheetg[33,33,33] > p_lightsheet1[33,33,33]
    @test p_lightsheetg[33,33,33] > p_lightsheetg[33,33,30]
end

@testset "ISM PSF" begin
    sampling = (0.04,0.04,0.120)
    sz = (128,128,128)
    pp_em = PSFParams(0.5,1.3,1.52; mode=ModeISM);
    pp_ex = PSFParams(pp_em; λ=0.488, aplanatic=aplanatic_illumination);
    pinhole = 0.001
    @time p_ism = psf(sz,pp_em; pp_ex=pp_ex, pinhole=pinhole, sampling=sampling);
    pp_em_conf = PSFParams(0.5,1.3,1.52; mode=ModeConfocal);
    @time p_conf = psf(sz,pp_em_conf; pp_ex=pp_ex, pinhole=pinhole, sampling=sampling);
    # @test ctr_test(p_ism[13], p_conf, 0.0001)
    @test isapprox(p_ism[13], p_conf; rtol=0.000001)
end

@testset "Two-Photon PSF" begin
    sampling = (0.04,0.04,0.120)
    sz = (128,128,128)
    pp_ex = PSFParams(0.800,1.3,1.52; aplanatic=aplanatic_illumination, mode=Mode2Photon);
    @time p_2p = psf(sz,pp_ex; sampling=sampling);
    pp_ex = PSFParams(0.800,1.3,1.52; mode=ModeWidefield, aplanatic=aplanatic_illumination);
    @time p_wf = psf(sz,pp_ex; sampling=sampling);
    # @test ctr_test(p_ism[13], p_conf, 0.0001)
    @test isapprox(p_2p, abs2.(p_wf); rtol=0.000001)
end

@testset "4Pi PSF" begin
    sampling = (0.04,0.04,0.04)
    sz = (128,128,128)
    pp_em = PSFParams(0.5,1.3,1.52; mode=Mode4Pi, pol=pol_x);
    pp_ex = PSFParams(0.800,1.3,1.52; aplanatic=aplanatic_illumination, mode=Mode4Pi, pol=pol_x);
    @time p_4pi = psf(sz,pp_ex; pp_em=pp_em, sampling=sampling, pinhole=2.0, ex_modifier=PointSpreadFunctions.modify_ident);
    @test isapprox(sum(p_4pi[:,:,64]), 2.376, rtol=0.01)
end

@testset "2D vs. 3D" begin
    sampling = (0.04,0.04)
    sz = (128,128)
    pp = PSFParams(0.5,1.3,1.52; mode=ModeWidefield, method=PointSpreadFunctions.MethodRichardsWolf, pol=pol_x);
    @test psf(sz, pp; sampling=sampling) == psf((sz...,1),pp; sampling=(sampling...,1))[:,:,1]
    pp = PSFParams(0.5,1.3,1.52; mode=ModeWidefield, method=PointSpreadFunctions.MethodPropagateIterative, pol=pol_x);
    @test psf(sz, pp; sampling=sampling) == psf((sz...,1),pp; sampling=(sampling...,eps(eltype(sampling))))[:,:,1]
    pp = PSFParams(0.5,1.3,1.52; mode=ModeWidefield, method=PointSpreadFunctions.MethodShell, pol=pol_x);
    @test psf(sz, pp; sampling=sampling) == psf((sz...,1),pp; sampling=(sampling...,eps(eltype(sampling))))[:,:,1]
    pp = PSFParams(0.5,1.3,1.52; mode=ModeWidefield, method=PointSpreadFunctions.MethodPropagate, pol=pol_x);
    @test psf(sz, pp; sampling=sampling) == psf((sz...,1),pp; sampling=(sampling...,eps(eltype(sampling))))[:,:,1]
    pp = PSFParams(0.5,1.3,1.52; mode=ModeWidefield, method=PointSpreadFunctions.MethodSincR, pol=pol_x);

    pp = PSFParams(0.5,1.3,1.52; mode=ModeConfocal, pol=pol_x);
    @test psf(sz, pp; pp_ex=pp, pinhole=1.0, sampling=sampling) == psf((sz...,1),pp; pp_ex=pp, pinhole=1.0, sampling=(sampling...,1))[:,:,1]
    @test psf((16,16,16), PSFParams(0.5,1.4,1.52); sampling=(0.05,0.05,0.05)) != 0 # to test for a bug in NDTools.select_region
end

return
