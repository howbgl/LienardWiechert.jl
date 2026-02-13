using LienardWiechert
using Test

pos(t)  = [1.0, 0.05 * sin(2π * t),      0.01 * cos(2π * t)]
dpos(t) = [0.0, 0.05 * 2π * cos(2π * t), -0.01 * 2π * sin(2π * t)]

@testset "LienardWiechert.jl" begin
    # Write your tests here.
    t   = collect(0:0.1:2π)
    dat = vcat(t',hcat(pos.(t)...),hcat(dpos.(t)...))'
    x   = dat[:, 2]
    y   = dat[:, 3]
    z   = dat[:, 4]
    dxdt = dat[:, 5]
    dydt = dat[:, 6]
    dzdt = dat[:, 7]

    t_unordered = [0.0, 2.0, 1.0]

    @testset "WorldLine constructors" begin

        @test_throws ArgumentError WorldLine([0.0], x, y, z)
        @info "The following warning about \"No speed data\" is expected and can be ignored."
        @test_throws ArgumentError WorldLine(
            t_unordered, 
            zero(t_unordered) , 
            zero(t_unordered), 
            zero(t_unordered))
        @test_throws ArgumentError WorldLine(t, x,1e5 .* y, z, dxdt, 1e5 .* dydt, dzdt) # Speed > c
        @test WorldLine(t,x,y,z,dxdt,dydt,dzdt) == WorldLine(hcat(t,x,y,z,dxdt,dydt,dzdt))
        @test WorldLine(t,x,y,z .=1e-15,dxdt,dydt,dzdt) ≈ WorldLine(t,x,y,z,dxdt,dydt,dzdt)
    end
end
