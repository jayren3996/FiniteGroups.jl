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

function irreps(g::AbstractFiniteGroup, χ::Union{AbstractMatrix, CharacterTable})
    reg = regular_rep(g)
    [project_out_rep(g, χ[i, :], reg) for i = 1:size(χ, 1)]
end

function irreps(g::AbstractFiniteGroup, c::Characters)
    irreps(g, c.χ)
end

function irreps(g::AbstractFiniteGroup; R::Bool=false)
    R ? real_irreps(g) : irreps(g, charactertable(g))
end

#-------------------------------------------------------------------------------
# Project out irreps from regular representation
#-------------------------------------------------------------------------------
export regular_rep, proj_operator
"""
    regular_rep(g::AbstractFiniteGroup)

Return regular representation.
"""
function regular_rep(g::AbstractFiniteGroup)
    n = order(g)
    rep = Vector{SparseMatrixCSC{Int64, Int64}}(undef, n)
    all_ones = fill(1, n)
    @threads for i = 1:n
        rep[i] = sparse([g[i,j] for j=1:n], 1:n, all_ones, n, n)
    end
    rep
end
#-------------------------------------------------------------------------------
function proj_operator(r::AbstractVector{<:AbstractMatrix}, χ::AbstractVector{<:Number}, g::AbstractFiniteGroup)
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

function proj_operator(r::AbstractVector{<:AbstractMatrix}, χ::AbstractVector{<:Number})
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
    vsi = vs'
    classes = class(g)
    preg = Vector{Matrix{eltype(vs)}}(undef, length(classes))
    @threads for i = 1:length(classes)
        preg[i] = vsi * reg[classes[i][1]] * vs
    end
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
function single_complex_row(ct::CharacterTable)
    nrow = length(ct)
    row = fill(true, nrow)
    for i = 1:nrow
        row[i] || continue
        if !(eltype(ct[i]) <: Real)
            χc = conj(ct[i])
            for j = i+1:nrow
                ct[j] == χc && (row[j] = false)
            end
        end
    end
    (1:nrow)[row]
end

#-------------------------------------------------------------------------------
# Helper
#-------------------------------------------------------------------------------
export transform_rep
"""
    transform_rep(U, D1, V)

Transform the basis of the representationc `D1`.

The inputs is `U`, `D1`, and `V`, the output is:
    D₂(g) = U⋅D₁(g)⋅V
"""
function transform_rep(U::AbstractMatrix{Tu}, D1::AbstractVector{<:AbstractMatrix{Tr}}, V::AbstractMatrix{Tv}) where {Tu, Tr, Tv}
    D2 = Vector{Matrix{promote_type(Tu, Tr, Tv)}}(undef, length(D1))
    @threads for i = 1:length(D1)
        D2[i] = U * D1[i] * V
    end
    D2
end

function transform_rep(rep::AbstractVector{<:AbstractMatrix}, v::AbstractMatrix)
    if norm(v' * v - I) < 1e-10 
        transform_rep(v', rep, v)
    else
        transform_rep(inv(v), rep, v)
    end
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
export oplus, unitary_rep, equivalent_transform
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
function unitary_rep(rep::AbstractVector{<:AbstractMatrix})
    H = sum(m' * m for m in rep)
    e, v = eigen(Hermitian(H))
    es, vd = sqrt.(e), v'
    X = v * Diagonal(1 ./ es) * vd
    Xi = v * Diagonal(es) * vd
    transform_rep(Xi, rep, X)
end
#-------------------------------------------------------------------------------
"""
Return the transformation matrix between two equivalent irreducible representation
`rep1`, `rep2`, satisfying:
    rep2[i] = U⁻¹ * rep1[i] * U.
"""
function equivalent_transform(rep1::AbstractVector{<:AbstractMatrix}, rep2::AbstractVector{<:AbstractMatrix})
    D = size(rep1[1], 1)
    H = sum(kron(conj(rep2[i]), rep1[i]) for i = 1:length(rep1))
    e, v = eigen(H)
    @assert abs(e[end]) - abs(e[end-1]) > 1e-3 "Encounter possible degeneracy, implying reducibility."
    reshape(v[:, end], D, D) * sqrt(D)
end
#-------------------------------------------------------------------------------
"""
Check whether a group is legit
"""
function check_rep(g::AbstractFiniteGroup, r::AbstractVector{<:AbstractMatrix}; tol::Real=1e-7)
    Q = Matrix{Bool}(undef, order(g), order(g))
    for k = 1:order(g) 
        for l = 1:order(g)
            m = g[k, l]
            Q[k, l] = norm(r[m] - r[k] * r[l]) < tol
        end
    end
    all(Q)
end
#-------------------------------------------------------------------------------
function check_unitary(r::AbstractVector{<:AbstractMatrix}; tol::Real=1e-7)
    Q = Vector{Bool}(undef, length(r))
    for i = 1:length(r)
        Q[i] = norm(r[i]' * r[i] - I) < tol
    end
    all(Q)
end
