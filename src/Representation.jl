export irreps
"""
    irreps(g::FiniteGroup; χ::AbstractVector)

Get the irreducible representations.
"""
function irreps(g::AbstractFiniteGroup, χ::AbstractVector)
    isone(χ[1]) && return [χ[inclass(g, i)] * ones(1,1) for i = 1:order(g)]
    reg = regular_rep(g)
    project_out_rep(g, χ, reg)
end

function irreps(g::AbstractFiniteGroup, χ::AbstractMatrix)
    reg = regular_rep(g)
    [project_out_rep(g, χ[i, :], reg) for i = 1:size(χ, 1)]
end

function irreps(g::AbstractFiniteGroup; R::Bool=false)
    R ? real_irreps(g) : irreps(g, charactertable(g))
end

#-------------------------------------------------------------------------------
# Project out irreps from regular representation
#-------------------------------------------------------------------------------
export regular_rep, proj_operator
"""
    regular_rep(multab::AbstractMatrix{<:Integer})

Return regular representation.
"""
function regular_rep(g::AbstractFiniteGroup)
    n = order(g)
    [sparse([g[i,j] for j=1:n], 1:n, fill(1, n), n, n) for i=1:n]
end
#-------------------------------------------------------------------------------
function proj_operator(
    r::AbstractVector{<:AbstractMatrix}, 
    χ::AbstractVector{<:Number},
    g::AbstractFiniteGroup
)
    dtype = promote_type(eltype(χ), eltype.(r)...)
    m = zeros(dtype, size(r[1]))
    if length(χ) == length(class(g))
        for i = 1:order(g)
            m += conj(χ[inclass(g, i)]) * r[i]
        end
    elseif length(χ) == order(g)
        for i = 1:order(g)
            m += χ[i] * r[i]
        end
    else
        error("Invalid length $(length(χ)) for χ.")
    end
    m
end

function proj_operator(
    r::AbstractVector{<:AbstractMatrix}, 
    χ::AbstractVector{<:Number}
)
    dtype = promote_type(eltype(χ), eltype.(r)...)
    m = zeros(dtype, size(r[1]))
    for i = 1:length(r)
        m += χ[i] * r[i]
    end
    m
end
#-------------------------------------------------------------------------------
function project_out_rep(
    g::AbstractFiniteGroup, 
    χ::AbstractVector{<:Number}, 
    reg::AbstractVector{<:AbstractMatrix};
    tol::Real=1e-7
)
    D = Int(χ[1])
    isone(D) && return [χ[inclass(g, i)] * ones(1,1) for i = 1:order(g)]
    vs = begin
        e, v = proj_operator(reg, χ, g) |> Hermitian |> eigen
        v[:, e .> 1e-2]
    end
    preg = [vs' * reg[cls[1]] * vs for cls in class(g)]
    v0 = find_least_degen(preg, D, tol=tol)
    P = krylov_space(reg[2:order(g)], vs*v0, D)
    rep = transform_rep(P', reg, P)
    isone(check_real_rep(g, χ)) ? real_rep(rep) : rep
end
#-------------------------------------------------------------------------------
function find_least_degen(rep, D::Integer; dim::Integer=D^2, tol::Real=1e-7)
    min_degen = dim
    min_vs = nothing
    for mat in rep
        e, v = eigen(mat)
        e_split = spectrum_split(e, tol=tol)
        for ei in e_split
            if length(ei) == D
                return v[:,ei[1]]
            elseif length(ei) < min_degen
                min_degen = length(ei)
                min_vs = v[:, ei]
            end
        end
    end
    # If single representation matrix is not enough to get rid of ther degeneracy,
    # we then calculate the rep. mats. restriceted to the degenerate space.
    min_vs = svd(min_vs).U
    rep2 = (min_vs' * mat * min_vs for mat in rep)
    min_vs * find_least_degen(rep2, D, dim=min_degen, tol=tol)
end
#-------------------------------------------------------------------------------
"""
    krylov_space(mats::AbstractVector{<:AbstractMatrix}, v0::AbstractVector, n::Integer; tol)

Find the krylov space of the group action from initial vector `v0`.
"""
function krylov_space(
    mats::AbstractVector{<:AbstractMatrix}, 
    v0::AbstractVector, 
    n::Integer; 
    tol::Real=1e-3
)
    isone(n) && return reshape(v0, :, 1)
    r, i, nmat = 1, 0, length(mats)
    vs = v0
    while r < n
        i = mod(i, nmat) + 1
        vs = hcat(vs, mats[i] * vs)
        U, S = svd(vs)
        pos = S .> tol
        vs = U[:, pos]
        r = size(vs, 2)
    end
    @assert r == n "Dimension of Krylov space not match: expect d = $n, get d = $r"
    vs
end

#-------------------------------------------------------------------------------
# Realify
#-------------------------------------------------------------------------------
export check_real_rep, real_rep, real_irreps
"""
    check_real_rep(g::AbstractFiniteGroup, χ)

Check whether a representation is real/complex/pseudo-real:
+1 : Real 
 0 : Complex 
-1 : Pseudo-real
"""
function check_real_rep(g::AbstractFiniteGroup, χ)
    n = sum(χ[inclass(g, g[i,i])] for i=1:order(g)) / order(g)
    if abs(n-1) < 1e-7          # Real 
        1
    elseif abs(n+1) < 1e-7      # Complex 
        -1
    elseif abs(n) < 1e-7        # Pseudo-real
        0
    else
        error("Sum of χ(g²) = $n ≠ ±g,0")
    end
end
#-------------------------------------------------------------------------------
"""
    real_rep(r::AbstractVector{<:AbstractMatrix})

Return real representation for a complex-value real representation.

For a representation D(g), it is a real representation means that
    S⋅D̄(g)⋅S⁺ = D(g) => Sᵀ = S,
Find the sqaure root of S, i.e., 
    S = W⋅W => W⋅D̄(t)⋅W⁺ = W⁺⋅D(g)⋅W
The real representation is W⁺⋅D(g)⋅W. 
"""
function real_rep(r::AbstractVector{<:AbstractMatrix})
    R = sum(kron(m, m) for m in r)
    v = eigen(R).vectors
    U = reshape(v[:, end], size(r[1])) * sqrt(size(r[1], 1))
    W = unitarysqrt(U)
    Wi = W'
    rep = Vector{Matrix{Float64}}(undef, length(r))
    Threads.@threads for i = 1:length(r)
        rep[i] = real.(Wi * r[i] * W)
    end
    rep
end
#-------------------------------------------------------------------------------
"""
    real_irreps(g::AbstractFiniteGroup)

Get all irreducible real representations of g.
"""
function real_irreps(g::AbstractFiniteGroup)
    ct = charactertable(g)
    rows = single_complex_row(ct)
    [real_irreps(g, ct[r, :]) for r in rows]
end

function real_irreps(g::AbstractFiniteGroup, χ::AbstractVector)
    rep = irreps(g, χ)
    if promote_type(typeof.(χ)...) <: Real
        rep
    else
        [[real(m) imag(m); -imag(m) real(m)] for m in rep]
    end
end
#-------------------------------------------------------------------------------
function single_complex_row(ct::AbstractMatrix)
    nrow = size(ct, 1)
    row = fill(true, nrow)
    for i = 1:nrow
        row[i] || continue
        if !(promote_type(typeof.(ct[i, :])...) <: Real)
            χc = conj(ct[i, :])
            for j = i+1:nrow
                (norm(ct[j, :] .- χc) < 1e-7) && (row[j] = false)
            end
        end
    end
    (1:nrow)[row]
end

#-------------------------------------------------------------------------------
# Helper
#-------------------------------------------------------------------------------
export transform_rep
function transform_rep(
    u::AbstractMatrix{Tu}, 
    rep::AbstractVector{<:AbstractMatrix{Tr}}, 
    v::AbstractMatrix{Tv}
) where {Tu, Tr, Tv}
    new_rep = Vector{Matrix{promote_type(Tu, Tr, Tv)}}(undef, length(rep))
    Threads.@threads for i = 1:length(rep)
        new_rep[i] = u * rep[i] * v
    end
    new_rep
end
#-------------------------------------------------------------------------------
"""
Eigen decomposition of unitary matrix.

The eigen vectors are unitary.
"""
function unitaryeigen(U::AbstractMatrix)
    val, vec = eigen(U)
    # Require the eigen vector of a unitary matrix to be unitary.
    # Here we enforce this by orthogonalized the vectors in the degenerate space.
    for s in spectrum_split(val)
        isone(length(s)) && continue
        vec[:, s] .= svd(vec[:, s]).U
    end
    val, vec
end
#-------------------------------------------------------------------------------
function unitarysqrt(U::AbstractMatrix)
    val, vec = unitaryeigen(U)
    sval = @. sqrt(Complex(val))
    vec * Diagonal(sval) * vec'
end

#-------------------------------------------------------------------------------
# Misc
#-------------------------------------------------------------------------------
export oplus, unitary_rep
"""
Direct sum of matrices or representations.
"""
function oplus(mats::AbstractMatrix...)
    dn = size.(mats, 1)
    D = sum(dn)
    dtype = promote_type(eltype.(mats)...)
    m = zeros(dtype, D, D)
    s = 1
    for i in 1:length(mats)
        d = dn[i]
        m[s:s+d-1, s:s+d-1] .= mats[i]
        s += d
    end
    m
end

function oplus(reps::AbstractVector{<:AbstractMatrix{<:Number}}...)
    [oplus((rep[i] for rep in reps)...) for i=1:length(reps[1])]
end
#-------------------------------------------------------------------------------
"""
Convert a representation to unitary representation.

Construct the matrix H = ∑D⁺(g)D(g) = v⋅e⋅v⁺ = (X⋅X⁺)⁻¹ = (X⁻¹)⁺⋅X⁻¹ 
    => X = v⋅(√e)⁻¹⋅v⁺, X⁻¹ = v⋅(√e)⋅v⁺
The new representation D'(g) := X⁻¹⋅D(g)⋅X is unitary:
    D'⁺(g)⋅D'(g) = X⁺⋅D(g)⋅(X⁻¹)⁺⋅X⁻¹⋅D(g)⋅X = X⁺⋅H⋅X = I.
"""
function unitary_rep(rep)
    H = sum(m' * m for m in rep)
    e, v = eigen(Hermitian(H))
    es, vd = sqrt.(e), v'
    X = v * Diagonal(1 ./ es) * vd
    Xi = v * Diagonal(es) * vd
    transform_rep(Xi, rep, X)
end
#-------------------------------------------------------------------------------
"""
Check whether a group is legit
"""
function check_rep(g::AbstractFiniteGroup, r; tol=1e-7)
    Threads.@threads for k = 1:order(g) 
        for l = 1:order(g)
            m = g[k, l]
            norm(r[m] - r[k] * r[l]) > tol && return false
        end
    end
    true
end
#-------------------------------------------------------------------------------
function check_unitary(r; tol=1e-7)
    Threads.@threads for m in r
        norm(m' * m - I) > tol && return false
    end
    true
end