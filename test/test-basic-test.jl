@testset "MassPointModelAD.jl" begin
    @test isdefined(MassPointModelAD, :MPMParams)
    @test isdefined(MassPointModelAD, :MPM_rhs)
    @test isdefined(MassPointModelAD, :simulate)
end
