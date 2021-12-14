include("$(@__DIR__)/../src/FiniteGroups.jl")
using .FiniteGroups, LinearAlgebra, Test


for gi = 1:32
    g = pointgroup(gi)
    χ = charactertable(g, method="Burnside")
    R = irreps(g, χ)
    @testset verbose = true "$(name(g)) Test" begin
        @testset "Character" begin
            # Test dimensions
            @test size(χ) == (length(class(g)), length(class(g)))
            # Test Orthogonality
            for i = 1:length(class(g)), j = 1:length(class(g))
                dt = sum(conj(χ[i,k]) * χ[j,k] * mult(g, k) for k = 1:length(class(g)))/order(g)
                dt -= (i==j)
                @test abs(dt) < 1e-5
            end
        end
        @testset "Irreps" begin
            for j = 1:size(χ, 1)
                rep = R[j]
                chi = χ[j, :]
                
                for k = 1:order(g)
                    @test chi[inclass(g, k)] ≈ tr(rep[k]) atol=1e-7     # check trace
                    @test norm(rep[k] * rep[k]' - I) < 1e-7             # check unitary
                end
                #test multiplication
                for k = 1:order(g), l = 1:order(g)
                    m = g[k, l]
                    @test rep[m] ≈ rep[k] * rep[l] atol=1e-7
                end
            end
        end
    end
end
