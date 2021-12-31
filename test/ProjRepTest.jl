include("$(@__DIR__)/../src/FiniteGroups.jl")
using .FiniteGroups, LinearAlgebra, Test
import .FiniteGroups: double_group, check_rep, check_proj_coeff

for gi = 1:32
    g = pointgroup(gi)
    ct = charactertable(g)
    @testset verbose = true "$(name(g)) Test" begin
        for i = 2:size(ct, 1)
            χ = ct[i, :]
            (promote_type(typeof.(χ)...) <: Real && χ[1] == 1) || continue
            @testset "$(repname(g, i))" begin
                bg, coeff = double_group(g, Int.(χ))
                reps = proj_reps(bg, coeff, 2)
                for rep in reps
                    for k = 1:order(bg)
                        @test norm(rep[k] * rep[k]' - I) < 1e-7             # check unitary
                    end
                    @test sum(size(r[1], 1)^2 for r in rep) == order(bg)
                end
                rreps = proj_reps(bg, coeff, 2, real=true)
                for rep in rreps
                    @test promote_type(eltype.(rep)...) == Float64
                    for k = 1:order(bg)
                        @test norm(rep[k] * rep[k]' - I) < 1e-7             # check unitary
                    end
                    @test check_proj_coeff(g, rep, coeff, 2)
                end
            end
        end
    end
end