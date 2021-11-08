include("$(@__DIR__)/../src/FiniteGroups.jl")
using .FiniteGroups, LinearAlgebra, Test

@testset "Character Table" begin
    for i = 1:32
        g = pointgroup(i)
        χ = charactertable(g, method="Burnside")
        # Test Dimension
        @test size(χ) == (length(class(g)), length(class(g)))
        print("\nGroup $(name(g)) character table : ")
        display(χ)
        print("\n")
        # Test orthogonalizty
        for i = 1:length(class(g)), j = 1:length(class(g))
            dt = sum(conj(χ[i,k]) * χ[j,k] * mult(g, k) for k = 1:length(class(g)))/order(g)
            dt -= (i==j)
            @test abs(dt) < 1e-5
        end
    end
end

