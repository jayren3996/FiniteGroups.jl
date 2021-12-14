export irreps
"""
    irreps(g::FiniteGroup; χ::AbstractVector)

Get the irreducible representations.
"""
function irreps(g::AbstractFiniteGroup, χ::AbstractVector)
    if isone(round(Int64, real(χ[1])))
        return [χ[inclass(g, i)] * ones(1,1) for i = 1:order(g)]
    end
    reg = regular_rep(g)
    prep(g, χ, reg)
end

function irreps(g::AbstractFiniteGroup, χ::AbstractMatrix)
    reg = regular_rep(g)
    [prep(g, χ[i, :], reg) for i = 1:size(χ, 1)]
end

irreps(g::AbstractFiniteGroup) = irreps(g, charactertable(g))
#-------------------------------------------------------------------------------
# Project out irreps from regular representation
#-------------------------------------------------------------------------------
function prep(
    g::AbstractFiniteGroup, 
    χ::AbstractVector{<:Number}, 
    reg::AbstractVector{<:AbstractMatrix};
    tol::Real=1e-7
)
    D = round(Int64, real(χ[1]))
    isone(D) && return [χ[inclass(g, i)] * ones(1,1) for i = 1:order(g)]
    vs = begin
        e, v = proj_operator(reg, χ, g) |> Hermitian |> eigen
        v[:, e .> 1e-2]
    end
    preg = [vs' * reg[cls[1]] * vs for cls in class(g)]
    v0 = find_least_degen(preg, D, tol=tol)
    P = krylov_space(reg[2:order(g)], vs*v0, D)
    rep = Vector{Matrix{eltype(P)}}(undef, order(g))
    Threads.@threads for i = 1:order(g)
        rep[i] = P' * reg[i] * P
    end
    isone(check_real_rep(g, χ)) && return real_rep(rep)
    rep
end

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
# Helpers
#-------------------------------------------------------------------------------
export check_real_rep
function check_real_rep(g::AbstractFiniteGroup, χ)
    n = sum(χ[inclass(g, g[i,i])] for i=1:order(g)) / order(g)
    if abs(n-1) < 1e-7
        1
    elseif abs(n+1) < 1e-7
        -1
    elseif abs(n) < 1e-7
        0
    else
        error("Sum of χ(g²) = $n ≠ ±g,0")
    end
end

export real_rep
function real_rep(r::AbstractVector{<:AbstractMatrix})
    R = sum(kron(m, m) for m in r)
    e, v = eigen(R)
    U = reshape(v[:, end], size(r[1])) * sqrt(size(r[1], 1))
    val, vec = eigen(U)
    sval = sqrt.(val)
    W = vec * Diagonal(sval) * inv(vec)
    Wi = inv(W)
    rep = Vector{Matrix{Float64}}(undef, length(r))
    Threads.@threads for i = 1:length(r)
        rep[i] = real.(Wi * r[i] * W)
    end
    rep
end

export regular_rep
"""
    regular_rep(multab::AbstractMatrix{<:Integer})

Return regular representation.
"""
function regular_rep(g::AbstractFiniteGroup)
    n = order(g)
    [sparse([g[i,j] for j=1:n], 1:n, fill(1, n), n, n) for i=1:n]
end

export proj_operator
function proj_operator(
    r::AbstractVector{<:AbstractMatrix}, 
    χ::AbstractVector{<:Number},
    g::AbstractFiniteGroup
)
    dtype = promote_type(eltype(χ), eltype.(r)...)
    m = zeros(dtype, size(r[1]))
    for i = 1:order(g)
        m += conj(χ[inclass(g, i)]) * r[i]
    end
    m
end

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