include("$(@__DIR__)/../src/FiniteGroups.jl")
using .FiniteGroups, LinearAlgebra, Test

g = pointgroup("Oh")
rreps = irreps(g, R=true)

rep = oplus(
    rreps[1],   # 1
    rreps[2],   # 2
    rreps[2],   # 3
    rreps[2],   # 4
    rreps[2],   # 5
    rreps[5],   # 6
    rreps[6],   # 7
    rreps[8],   # 8
    rreps[8],   # 9
    rreps[9],   # 10
    rreps[10]   # 11
)
ru = Symmetric(rand(size(rep[1])...)) |> eigvecs
rep = transform_rep(ru', rep, ru)
res = block_decomposition(rep, g)

@testset "Oh Decomp" begin
    @test length(res) == 11
    @test norm(transform_rep(res[1]', rep, res[1]) .- rreps[1]) < 1e-10
    for i = 2:5
        @test norm(transform_rep(res[i]', rep, res[i]) .- rreps[2]) < 1e-10
    end
    for i = 6:6
        @test norm(transform_rep(res[i]', rep, res[i]) .- rreps[5]) < 1e-10
    end
    for i = 7:7
        @test norm(transform_rep(res[i]', rep, res[i]) .- rreps[6]) < 1e-10
    end
    for i = 8:9
        @test norm(transform_rep(res[i]', rep, res[i]) .- rreps[8]) < 1e-10
    end
    for i = 10:10
        @test norm(transform_rep(res[i]', rep, res[i]) .- rreps[9]) < 1e-10
    end
    for i = 11:11
        @test norm(transform_rep(res[i]', rep, res[i]) .- rreps[10]) < 1e-10
    end
end