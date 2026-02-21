using Test
using COSMORSAmino

@testset "Boltzmann" begin
    energies = [0.0, 1000.0]
    w = boltzmann_weights(energies)
    @test sum(w) ≈ 1.0
end