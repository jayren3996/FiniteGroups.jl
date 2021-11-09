include("$(@__DIR__)/../src/FiniteGroups.jl")
using .FiniteGroups, LinearAlgebra, Test


@testset "Irreps Test" begin
    for i = 1:32
        g = pointgroup(i)
        chis = charactertable(g)
        reps = irreps(g, χ=chis)
        for j = 1:size(chis, 1)
            rep = reps[j]
            chi = chis[j, :]
            # check trace
            for k = 1:order(g)
                @test chi[inclass(g, k)] ≈ tr(rep[k]) atol=1e-7
            end
            #test multiplication
            for k = 1:order(g), l = 1:order(g)
                m = g[k, l]
                @test rep[m] ≈ rep[k] * rep[l] atol=1e-7
            end
        end
    end
end
