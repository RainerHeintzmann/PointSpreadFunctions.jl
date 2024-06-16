using PointSpreadFunctions
using Revise


λ_em = 0.5; NA = 1.4; n = 1.52
λ_ex = 0.488 # only needed for some PointSpreadFunctions, such as confocal, ISM or TwoPhoton


sz = (256, 256, 256)
sampling = (0.050,0.050,0.050)
#print(pp)
#println("Is PointSpreadFunctions loaded? ", isdefined(Main, :PointSpreadFunctions))
#println("Is MethodCZT available? ", isdefined(PointSpreadFunctions, :MethodCZT))

pp1 = PSFParams(λ_em, NA, n; method= MethodCZT);@time p_czt = psf(sz, pp1; sampling=sampling)
pp2 = PSFParams(λ_em, NA, n; method= MethodRichardsWolf); @time p_RW = psf(sz, pp2; sampling=sampling)
#Compare the calculation time of APSF using Method CZT and Method Richards-Wolf.
