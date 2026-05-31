using FiniteGroups, LinearAlgebra, Test
import FiniteGroups: Characters, CharacterTable, check_rep, real_irreps,
                     equivalent_transform, check_group, check_unitary, oplus,
                     matrix, double_group, cover_group

# Quaternion group Q8 built from the quaternion algebra; its 2D irrep is
# pseudo-real (Frobenius–Schur −1), which is the interesting case for realify.
function q8_table()
    umul = Dict{Tuple{Int,Int},Tuple{Int,Int}}()
    for u in 0:3; umul[(0,u)] = (1,u); umul[(u,0)] = (1,u); end
    umul[(1,1)] = (-1,0); umul[(2,2)] = (-1,0); umul[(3,3)] = (-1,0)
    umul[(1,2)] = (1,3);  umul[(2,3)] = (1,1);  umul[(3,1)] = (1,2)
    umul[(2,1)] = (-1,3); umul[(3,2)] = (-1,1); umul[(1,3)] = (-1,2)
    enc = [(1,0),(-1,0),(1,1),(-1,1),(1,2),(-1,2),(1,3),(-1,3)]
    idx = Dict(enc[i] => i for i in 1:8)
    tab = Matrix{Int}(undef, 8, 8)
    for a in 1:8, b in 1:8
        sa, ua = enc[a]; sb, ub = enc[b]; s, u = umul[(ua, ub)]
        tab[a, b] = idx[(sa * sb * s, u)]
    end
    tab
end

@testset "Review fixes" begin
    @testset "permutation validates its input" begin
        @test_throws ArgumentError permutation([1 2])        # 1 ↦ 2 only: not a bijection
        @test_throws ArgumentError permutation([1 2; 1 3])   # source 1 mapped twice
        @test_throws ArgumentError permutation([1 2 3])      # not two columns
        # valid permutations still build and (crucially) display without hanging
        @test string(permutation([1 2; 2 1])) == "Cycle: (1, 2)"
        @test string(cycles((1, 2, 3))) == "Cycle: (1, 2, 3)"
        @test string(cycles()) == "Cycle: "                  # identity
    end

    @testset "permutationgroup(n) edge cases" begin
        @test order(permutationgroup(0)) == 1
        @test order(permutationgroup(1)) == 1
        @test order(permutationgroup(2)) == 2
        @test order(permutationgroup(4)) == 24
    end

    @testset "FiniteGroup: square check and matrix views" begin
        @test_throws DimensionMismatch FiniteGroup([1 2 3; 2 1 3])      # 2×3
        g = FiniteGroup(view([1 2; 2 1], :, :))                        # a SubArray
        @test order(g) == 2
        @test g.multab isa Matrix{Int}
    end

    @testset "Burnside handles the one-class group C1" begin
        c1 = FiniteGroup(reshape([1], 1, 1), "C1")
        @test Matrix(charactertable(c1; method = :burnside)) == reshape([1], 1, 1)
        @test Matrix(charactertable(c1; method = :dixon))    == reshape([1], 1, 1)
    end

    @testset "CharacterTable accepts a concrete-eltype vector" begin
        g = FiniteGroup([1 2 3; 2 3 1; 3 1 2], "C3")
        v = [Characters(g, ComplexF64[1, 1, 1])]               # Vector{Characters{ComplexF64,…}}
        @test v isa Vector{<:Characters}
        @test !(v isa Vector{Characters})                      # genuinely the narrow case
        @test length(CharacterTable(v)) == 1
    end

    @testset "real_irreps uses the Frobenius–Schur indicator" begin
        g = FiniteGroup(q8_table(), "Q8")
        # Q8's 2D irrep is pseudo-real; passing its (real-valued) character as a
        # real-typed vector must still produce the 4-dimensional real block.
        rr = real_irreps(g, Float64[2, -2, 0, 0, 0])           # one real irrep = vector of matrices
        @test size(rr[1], 1) == 4
        @test check_rep(g, rr)
        @test all(m -> norm(m'm - I) < 1e-7, rr)               # orthogonal/real
    end

    @testset "block_decomposition of a real rep with a complex irrep" begin
        # C4 has a complex-conjugate irrep pair; a real representation containing
        # it must decompose into *real* orthogonal blocks when R=true.
        g  = FiniteGroup([1 2 3 4; 2 3 4 1; 3 4 1 2; 4 1 2 3], "C4")
        rr = real_irreps(g)                                    # 1D, 1D, 2D(rotation)
        t2 = findfirst(r -> size(r[1], 1) == 2, rr)
        rep0 = oplus(vcat(rr, [rr[t2]])...)                    # 2D (complex) block twice
        D = size(rep0[1], 1)
        Q = Matrix(qr(reshape(collect(1.0:D^2), D, D) + 7I).Q) # deterministic real orthogonal
        rep = transform_rep(Q', rep0, Q)
        res = block_decomposition(rep, g; R = true)

        # The bug: the complex irrep's component came back as a genuinely complex
        # 1D space (non-zero imaginary part). After the fix every block is
        # real-valued and the complex pair is realified to a 2D block.
        @test all(S -> sum(norm, imag.(S)) < 1e-9, res)
        @test any(S -> size(S, 2) == 2, res)
        @test sum(size(S, 2) for S in res) == D
        Q2 = hcat(res...)
        @test norm(Q2' * Q2 - I) < 1e-8
        blocks = [[real(S' * rep[k] * S) for k in 1:order(g)] for S in res]
        for k in 1:order(g)                                    # Q2 block-diagonalizes rep
            recon = oplus([blocks[i][k] for i in 1:length(res)]...)
            @test norm(Q2' * rep[k] * Q2 - recon) < 1e-7
        end
        @test all(B -> check_rep(g, B), blocks)                # each block is a valid irrep
    end

    @testset "repname still labels point-group irreps" begin
        @test repname(pointgroup("C4v"), 1) isa AbstractString
    end
end

@testset "Review fixes (round 2)" begin
    @testset "equivalent_transform on a 1D irrep" begin
        g = FiniteGroup([1 2 3; 2 3 1; 3 1 2], "C3")
        d = irreps(g, charactertable(g))[1]                    # a 1D irrep
        U = equivalent_transform(d, d)                         # used to throw BoundsError
        @test size(U) == (1, 1)
    end

    @testset "check_group gives an informative error, not BoundsError" begin
        @test_throws ErrorException check_group([1 2; 2 3])    # entry 3 out of range
        @test !(try check_group([1 2; 2 3]); catch e; e isa BoundsError end)
        @test check_group([1 2 3; 2 3 1; 3 1 2])               # valid C3
    end

    @testset "permutationgroup rejects negative n" begin
        @test_throws ArgumentError permutationgroup(-1)
        @test order(permutationgroup(0)) == 1
    end

    @testset "oplus rejects mismatched lengths" begin
        @test_throws DimensionMismatch oplus([ones(1,1), ones(1,1)],
                                             [ones(1,1), ones(1,1), ones(1,1)])
    end

    @testset "check_rep / check_unitary validate shape" begin
        g = FiniteGroup([1 2; 2 1], "C2")
        @test check_rep(g, [ones(1,1)]) == false               # too few matrices
        @test check_unitary([[1.0 0; 0 1; 0 0]]) == false      # non-square isometry
        @test check_unitary([[0.0 1; 1 0]]) == true
    end

    @testset "Permutation == / hash work in Set" begin
        @test permutation([1 2; 2 1]) == permutation([2 1; 1 2])
        @test length(Set([permutation([1 2; 2 1]), permutation([2 1; 1 2])])) == 1
    end

    @testset "point-group accessors return copies (no global aliasing)" begin
        g = pointgroup("C2v")
        m0 = matrix(g)[1][1, 1]
        let m = matrix(g); m[1][1, 1] += 100 end
        @test matrix(g)[1][1, 1] == m0
        o0 = operation(g)[1]
        let o = operation(g); o[1] = "ZZZ" end
        @test operation(g)[1] == o0
    end

    @testset "cover_group validates its coefficient matrix" begin
        c3 = FiniteGroup([1 2 3; 2 3 1; 3 1 2], "C3")
        c2 = FiniteGroup([1 2; 2 1], "C2")
        @test_throws ArgumentError cover_group(c3, [0 0 0; 0 0 1; 0 0 0], 2)  # not a 2-cocycle
        @test_throws DimensionMismatch cover_group(c2, zeros(Int, 3, 3), 2)   # wrong shape
        @test_throws ArgumentError cover_group(c2, zeros(Int, 2, 2), 0)       # p ≤ 0
    end

    @testset "proj_reps dedupes equivalent restricted reps" begin
        c2 = FiniteGroup([1 2; 2 1], "C2")
        @test length(proj_reps(c2, zeros(Int, 2, 2), 2)) == 2   # was 4 (Z₂ duplicates)
    end

    @testset "irreps(pointgroup; R=true) has no conjugate-pair duplicates" begin
        for nm in ("C3", "C4", "C6")
            @test length(irreps(pointgroup(nm); R=true)) == length(real_irreps(pointgroup(nm)))
        end
        @test length(irreps(pointgroup("Oh"); R=true)) == 10    # all-real group: unchanged
    end
end
