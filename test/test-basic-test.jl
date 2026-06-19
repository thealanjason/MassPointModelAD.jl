@testset "MassPointModelAD.jl" begin
    @test isdefined(MassPointModelAD, :MPMParams)
    @test isdefined(MassPointModelAD, :DEM)
    @test isdefined(MassPointModelAD, :load_dem)
    @test isdefined(MassPointModelAD, :save_dem)
    @test isdefined(MassPointModelAD, :dem_from_func)
    @test isdefined(MassPointModelAD, :MPM_rhs)
    @test isdefined(MassPointModelAD, :simulate)
end
