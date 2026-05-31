using FiniteGroups, LinearAlgebra, Test

# Multiset equality of character rows (irrep order may differ between methods;
# columns are class-aligned because both use the same group's class ordering).
function same_chartable(A::AbstractMatrix, B::AbstractMatrix; tol=1e-6)
    size(A) == size(B) || return false
    n = size(A, 1)
    used = falses(n)
    for i = 1:n
        j = findfirst(k -> !used[k] && norm(ComplexF64.(A[i, :]) - ComplexF64.(B[k, :])) < tol, 1:n)
        j === nothing && return false
        used[j] = true
    end
    true
end

@testset "Dixon vs bundled table (point groups)" begin
    for i = 1:32
        g = pointgroup(i)
        dx  = Matrix(charactertable(g; method=:dixon))
        ref = Matrix(charactertable(g))                 # bundled standard table
        @test same_chartable(dx, ref)
    end
end

@testset "Dixon vs Burnside (symmetric groups)" begin
    for nn = 2:5
        g = permutationgroup(nn)
        dx  = Matrix(charactertable(g; method=:dixon))
        ref = Matrix(charactertable(g; method=:burnside))   # explicit: default is now :dixon
        @test same_chartable(dx, ref)
    end
end

@testset "Dixon vs Burnside (FiniteGroup from multab)" begin
    # cyclic C4 and the Klein four-group, built from raw multiplication tables
    c4 = FiniteGroup([1 2 3 4; 2 3 4 1; 3 4 1 2; 4 1 2 3], "C4")
    v4 = FiniteGroup([1 2 3 4; 2 1 4 3; 3 4 1 2; 4 3 2 1], "V4")
    for g in (c4, v4)
        @test same_chartable(Matrix(charactertable(g; method=:dixon)),
                             Matrix(charactertable(g; method=:burnside)))
    end
end

@testset "Dixon table orthonormality" begin
    for i = 1:32
        g = pointgroup(i)
        ct = charactertable(g; method=:dixon)
        NC = length(class(g))
        for r = 1:NC, s = 1:NC
            val = sum(mult(g, k) * ct[r, k] * conj(ct[s, k]) for k = 1:NC) / order(g)
            @test abs(val - (r == s)) < 1e-6
        end
    end
end
