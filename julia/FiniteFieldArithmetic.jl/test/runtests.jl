using Test

@testset "FiniteFieldArithmetic" begin
    @testset "nimber fields" include("nimber.jl")
    @testset "nimber fields above GF(256)" include("nimber_large.jl")
    @testset "fast strategies for the large levels (J2a)" include("nimber_fast.jl")
    @testset "bulk and array operations (J8)" include("bulk.jl")
    @testset "constant times array in broadcasts (J8c)" include("constmul.jl")
    @testset "default strategy of inv" include("inv_default.jl")
    @testset "printing" include("printing.jl")
    @testset "static checks" begin
        include("aqua.jl")
        include("jet.jl")
    end
end
