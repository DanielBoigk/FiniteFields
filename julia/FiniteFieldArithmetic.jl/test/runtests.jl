using Test

@testset "FiniteFieldArithmetic" begin
    @testset "nimber fields" include("nimber.jl")
    @testset "printing" include("printing.jl")
    @testset "static checks" begin
        include("aqua.jl")
        include("jet.jl")
    end
end
