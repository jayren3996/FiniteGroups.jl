#-------------------------------------------------------------------------------
# Dixon's method: exact character tables over a finite field 𝔽ₚ.
#
# Pick a prime p ≡ 1 (mod n), n = group exponent, p > 2√|G|, and reduce the
# class-algebra structure constants mod p. The class matrices commute, so a
# generic combination of them has NC distinct eigenvalues in 𝔽ₚ; its
# eigenvectors are the central characters ωₖ = |Cₖ|·χ(gₖ)/χ(1). We normalize
# those to genuine characters mod p (mirroring the Burnside path in
# `Character.jl`, which is the same arithmetic over ℝ), then lift each value to
# an exact sum of complex roots of unity via the power map. The result avoids
# the floating-point tolerances of the Burnside method.
#
# Reference: J. D. Dixon, "High speed computation of group characters",
# Numer. Math. 10 (1967) 446-450.
#-------------------------------------------------------------------------------
export dixon

"""
    dixon(g::AbstractFiniteGroup [, h]) -> Matrix{ComplexF64}

Exact character table of `g` (rows = irreps, columns = conjugacy classes) via
Dixon's method over a finite field. `h` defaults to `class_multab(g)`. This is
the engine behind `charactertable(g; method = :dixon)`.
"""
function dixon(g::AbstractFiniteGroup, h::AbstractArray{<:Integer, 3})
    NC = size(h, 1)
    n, p = findprime(g)
    Z = primitive_root(n, p)
    invcls = inverse_class(g)
    table = Matrix{Int}(undef, NC, NC)          # χᵣ(gₖ) mod p
    for (r, ω) in enumerate(central_characters_modp(h, p))
        normalize_modp!(g, ω, invcls, p)
        table[r, :] .= ω
    end
    lift_table(g, table, p, n, Z)
end
dixon(g::AbstractFiniteGroup) = dixon(g, class_multab(g))

"""
Build a `CharacterTable` from Dixon's exact character matrix; the engine behind
`charactertable(g; method = :dixon)`.
"""
function dixon_table(g::AbstractFiniteGroup; tol::Real=1e-7)
    D = dixon(g)
    NC = size(D, 1)
    characters = Vector{Characters}(undef, NC)      # abstract elt: rows may mix real/complex
    for r = 1:NC
        characters[r] = Characters(g, beautify.(D[r, :], tol=tol))
    end
    CharacterTable(sort(characters))
end

"""
Conjugacy class of the inverse of (a representative of) each class.
"""
inverse_class(g::AbstractFiniteGroup) = [inclass(g, inv(g, class(g, i)[1])) for i = 1:length(class(g))]

#-------------------------------------------------------------------------------
# Central characters over 𝔽ₚ
#-------------------------------------------------------------------------------
"""
Common eigenvectors (central characters ωᵣ, mod p) of the class matrices.

The class matrices commute, so they are simultaneously diagonalizable over 𝔽ₚ.
We refine a list of common eigenspaces one class matrix at a time: restrict each
`hᵢ` to the current subspace `S` (solving `hᵢ S = S T` for the `d×d` matrix `T`),
take the eigenvalues of `T` as the roots of its characteristic polynomial over
𝔽ₚ, and split `S` by the nullspace of each `T − λI`. Iterating over all classes
drives every eigenspace down to dimension one — robust even for small `p`, where
no single combination of class matrices has enough distinct eigenvalues.
"""
function central_characters_modp(h::AbstractArray{<:Integer, 3}, p::Integer)
    NC = size(h, 1)
    spaces = [Matrix{Int}(I, NC, NC)]
    for i = 1:NC
        all(s -> isone(size(s, 2)), spaces) && break
        new = Matrix{Int}[]
        for S in spaces
            if isone(size(S, 2))
                push!(new, S)
                continue
            end
            T = restrict_modp(S, view(h, i, :, :), p)
            for λ in eigvals_modp(T, p)
                ns = nullspace_modp(T, λ, p)
                size(ns, 2) > 0 && push!(new, mod.(S * ns, p))
            end
        end
        spaces = new
    end
    @assert all(s -> isone(size(s, 2)), spaces) "Dixon: class algebra not fully split over 𝔽ₚ = $p."
    [s[:, 1] for s in spaces]
end

"""
Matrix `T` of the restriction of `M` to the column space of `S` (over 𝔽ₚ): the
unique `T` with `M·S = S·T`, found by expressing each column of `M·S` in the
basis `S`. Assumes `S` has full column rank and `col(M·S) ⊆ col(S)`.
"""
function restrict_modp(S::AbstractMatrix{<:Integer}, M::AbstractMatrix{<:Integer}, p::Integer)
    d = size(S, 2)
    A = hcat(mod.(S, p), mod.(M * S, p))        # [S | M·S]
    pivots = rref_modp!(A, p)
    T = zeros(Int, d, d)
    for (ri, pc) in enumerate(pivots)
        pc <= d || continue                     # pivots in the S-block index basis vectors
        for c = 1:d
            T[pc, c] = A[ri, d + c]
        end
    end
    T
end

"""
Turn a central-character eigenvector `ω` (mod p, defined up to scale) into the
character row `χ(gₖ)` of one irrep, in place.

`ωₖ = |Cₖ|·χ(gₖ)/χ(1)`, so the identity component `ω₁ = 1`. Scale to that, then
recover the degree from `d² = |G| / Σₖ ωₖ ω_{k̄}/|Cₖ|` (with `k̄` the inverse
class, since `conj(ωₖ) = ω_{k̄}`) and form `χ(gₖ) = d·ωₖ/|Cₖ|`. This mirrors
`normalize_chi!` in `Character.jl`, but over 𝔽ₚ.
"""
function normalize_modp!(
    g::AbstractFiniteGroup,
    ω::AbstractVector{<:Integer},
    invcls::AbstractVector{<:Integer},
    p::Integer
)
    NC = length(ω)
    a = invmod(mod(ω[1], p), p)                 # scale so ω₁ = 1
    for k = 1:NC
        ω[k] = mod(ω[k] * a, p)
    end
    s = 0                                        # s = Σₖ ωₖ ω_{k̄} / |Cₖ|
    for k = 1:NC
        s = mod(s + ω[k] * ω[invcls[k]] * invmod(mod(mult(g, k), p), p), p)
    end
    d2 = mod(mod(order(g), p) * invmod(s, p), p)
    d = sqrt_modp(d2, p)                         # smaller root ≤ p/2 = dimension
    for k = 1:NC                                 # ωₖ → χ(gₖ) = d·ωₖ/|Cₖ|
        ω[k] = mod(d * ω[k] * invmod(mod(mult(g, k), p), p), p)
    end
    ω
end

#-------------------------------------------------------------------------------
# Lift 𝔽ₚ character values to exact complex roots of unity
#-------------------------------------------------------------------------------
"""
Recover the complex character table from its reduction `table` mod p.

For class `k` of element order `m`, `χ(gₖ) = Σₗ aₗ exp(2πi l/m)` where the
eigenvalue multiplicities `aₗ ∈ 0:d` are obtained by an inverse DFT over the
power classes `gₖᵗ`: `aₗ = m⁻¹ Σₜ χ(gₖᵗ) z^{-lt} (mod p)`, with `z = Z^{n/m}` a
primitive m-th root of unity in 𝔽ₚ. Since `aₗ ≤ d < p/2`, the symmetric
representative recovers the exact integer multiplicity.
"""
function lift_table(g::AbstractFiniteGroup, table::AbstractMatrix{<:Integer}, p::Integer, n::Integer, Z::Integer)
    NC = size(table, 1)
    out = Matrix{ComplexF64}(undef, NC, NC)
    for k = 1:NC
        e = class(g, k)[1]
        m = order(g, e)
        z = powermod(Z, n ÷ m, p)               # primitive m-th root in 𝔽ₚ
        zpow = Vector{Int}(undef, m)            # zpow[j+1] = zʲ mod p, j = 0:m-1
        zpow[1] = 1
        for j = 2:m
            zpow[j] = mod(zpow[j-1] * z, p)
        end
        invm = invmod(mod(m, p), p)
        pcls = Vector{Int}(undef, m)            # classes of gₖᵗ, t = 0:m-1
        x = 1
        for t = 1:m
            pcls[t] = inclass(g, x)
            x = g[e, x]
        end
        for r = 1:NC
            val = 0.0 + 0.0im
            for l = 0:m-1
                acc = 0
                for t = 0:m-1
                    acc = mod(acc + table[r, pcls[t+1]] * zpow[mod(-l * t, m) + 1], p)
                end
                a = mod(acc * invm, p)
                a = a > p ÷ 2 ? a - p : a        # symmetric representative
                iszero(a) || (val += a * cis(2π * l / m))
            end
            out[r, k] = val
        end
    end
    out
end

#-------------------------------------------------------------------------------
# Linear algebra over the prime field 𝔽ₚ
#-------------------------------------------------------------------------------
"""
Basis (as columns) of the null space of `M - λI` over 𝔽ₚ.
"""
function nullspace_modp(M::AbstractMatrix{<:Integer}, λ::Integer, p::Integer)
    N = size(M, 1)
    A = Matrix{Int}(undef, N, N)
    for i = 1:N, j = 1:N
        A[i, j] = mod(M[i, j] - (i == j ? λ : 0), p)
    end
    pivots = rref_modp!(A, p)
    pivotset = Set(pivots)
    free = [c for c = 1:N if !(c in pivotset)]
    basis = zeros(Int, N, length(free))
    for (k, f) in enumerate(free)
        basis[f, k] = 1
        for (ri, pc) in enumerate(pivots)
            basis[pc, k] = mod(-A[ri, f], p)
        end
    end
    basis
end

"""
Eigenvalues of `M` in 𝔽ₚ, ascending: the roots of its characteristic polynomial.

We build the char poly with the division-free Samuelson–Berkowitz recurrence (so
it holds even when the size exceeds `p`, e.g. elementary-abelian groups) and read
off its roots by Horner-evaluating it across 𝔽ₚ. This replaces scanning a null
space at every field element with one `O(d³)` polynomial and an `O(p·d)` sweep.
"""
function eigvals_modp(M::AbstractMatrix{<:Integer}, p::Integer)
    c = charpoly_modp(M, p)
    λs = Int[]
    for λ = 0:p-1
        v = 0
        for ck in c
            v = mod(v * λ + ck, p)
        end
        iszero(v) && push!(λs, λ)
    end
    λs
end

"""
Characteristic polynomial of `M` over 𝔽ₚ via the Samuelson–Berkowitz algorithm,
returned high-degree-first (`c[1]·xⁿ + … + c[n+1]`, monic with `c[1] = 1`).

The recurrence is division-free, so it is valid for any prime `p` regardless of
the matrix size `n` — unlike Faddeev–LeVerrier, which divides by `1:n`.
"""
function charpoly_modp(M::AbstractMatrix{<:Integer}, p::Integer)
    n = size(M, 1)
    result = [1]                                    # coeffs so far, high-degree-first
    for i = 1:n
        v = Vector{Int}(undef, i + 1)               # generating column of the iᵗʰ Toeplitz block
        v[1] = 1
        v[2] = mod(-M[i, i], p)
        if i > 1
            R = view(M, i, 1:i-1)
            u = Vector{Int}(M[1:i-1, i])            # Aᵢ₋₁ʲ · Sᵢ, starting at j = 0
            for j = 0:i-2
                s = 0
                for t = 1:i-1
                    s = mod(s + R[t] * u[t], p)
                end
                v[3 + j] = mod(-s, p)               # −Rᵢ Aᵢ₋₁ʲ Sᵢ
                if j < i - 2                         # advance u ← Aᵢ₋₁·u for the next power
                    u2 = zeros(Int, i - 1)
                    for a = 1:i-1, b = 1:i-1
                        u2[a] = mod(u2[a] + M[a, b] * u[b], p)
                    end
                    u = u2
                end
            end
        end
        nr = zeros(Int, i + 1)                       # nr ← Tᵢ · result, Tᵢ lower-triangular Toeplitz of v
        for col = 1:i
            rc = result[col]
            iszero(rc) && continue
            for row = col:i+1
                nr[row] = mod(nr[row] + v[row - col + 1] * rc, p)
            end
        end
        result = nr
    end
    result
end

"""
Reduced row echelon form of `A` over 𝔽ₚ, in place. Returns the pivot columns.
"""
function rref_modp!(A::AbstractMatrix{<:Integer}, p::Integer)
    m, n = size(A)
    pivots = Int[]
    r = 1
    for c = 1:n
        r > m && break
        piv = 0
        for i = r:m
            iszero(A[i, c]) || (piv = i; break)
        end
        iszero(piv) && continue
        if piv != r
            for j = 1:n
                A[r, j], A[piv, j] = A[piv, j], A[r, j]
            end
        end
        ip = invmod(mod(A[r, c], p), p)
        for j = 1:n
            A[r, j] = mod(A[r, j] * ip, p)
        end
        for i = 1:m
            if i != r && !iszero(A[i, c])
                f = A[i, c]
                for j = 1:n
                    A[i, j] = mod(A[i, j] - f * A[r, j], p)
                end
            end
        end
        push!(pivots, c)
        r += 1
    end
    pivots
end

"""
Smaller modular square root of `a` in 𝔽ₚ (the one in `0:p÷2`). For a quadratic
residue `r²` with `2r < p` this returns `r`, so the integer dimension is exact.
"""
function sqrt_modp(a::Integer, p::Integer)
    a = mod(a, p)
    for r = 0:(p ÷ 2)
        mod(r * r, p) == a && return r
    end
    error("$a is not a quadratic residue mod $p (try a larger prime).")
end

"""
A primitive `n`-th root of unity in 𝔽ₚ (requires `n ∣ p-1`).
"""
function primitive_root(n::Integer, modp::Integer)
    isone(n) && return 1
    for i = 2:modp
        x = i
        ok = true
        for _ = 1:n-2
            x = mod(i * x, modp)
            isone(x) && (ok = false; break)
        end
        ok && isone(mod(i * x, modp)) && return i
    end
    error("no primitive $n-th root of unity in 𝔽_$modp.")
end

"""
Smallest prime `p ≡ 1 (mod n)` with `p > 2√|G|`, where `n` is the group
exponent. Returns `(n, p)`.
"""
function findprime(g::AbstractFiniteGroup)
    n = lcm([order(g, c[1]) for c in class(g)])
    i = ceil(Int, 2 * sqrt(order(g)))
    i = i^2 > 4 * order(g) ? i : i + 1
    while !(isprime_(i) && iszero(mod(i - 1, n)))
        i += 1
    end
    (n, i)
end

function isprime_(i::Integer)
    i < 2 && return false
    i == 2 && return true
    iseven(i) && return false
    d = 3
    while d * d <= i
        iszero(mod(i, d)) && return false
        d += 2
    end
    true
end
