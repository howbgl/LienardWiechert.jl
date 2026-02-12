using LienardWiechert
using Test

@testset "LienardWiechert.jl" begin
    # Write your tests here.
    t = [0.0, 1.0, 2.0]
    x = [0.0, 0.5, 1.0]
    y = [0.0, -0.5, -1.0]
    z = [1.0, 1.5, 2.0]

    t_unordered = [0.0, 2.0, 1.0]

    @testset "WorldLine constructors" begin

        @test_throws ArgumentError WorldLine(t_unordered, x, y, z)

        @test WorldLine(t, x, y, z) == WorldLine(hcat(t, x, y, z))
    end
end
