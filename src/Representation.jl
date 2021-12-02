export irreps
"""
    irreps(g::FiniteGroup; χ::Union{Nothing, AbstractMatrix}=nothing,)

Get the irreducible representations.
"""
function irreps(
    g::AbstractFiniteGroup;
    χ::Union{Nothing, AbstractMatrix}=nothing
)
    isnothing(χ) && ( χ = charactertable(g) )
    reg = regular_rep(g)
    [prep(g, χ[i, :], reg) for i = 1:size(χ, 1)]
end

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
    preg = begin
        vs = begin
            e, v = proj_operator(reg, χ, g) |> eigen
            v[:, e .> 1e-2]
        end
        @assert size(vs, 2) == D^2 "Dimension error: exprect $(D^2), get $(size(vs, 2))"
        [vs' * regmat * vs for regmat in reg]
    end
    isone(D) && return [beautify.(mi) for mi in preg]
    rep_gen = (preg[cls[1]] for cls in class(g))
    v0 = find_least_degen(rep_gen, D, tol=tol)
    vs = krylov_space(preg[2:order(g)], v0, D)
    [vs' * pregmat * vs for pregmat in preg]
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
    m = sum(conj(χ[inclass(g, i)]) * r[i] for i=1:order(g))
    Hermitian(Array(m))
end

function krylov_space(
    mats::AbstractVector{<:AbstractMatrix}, 
    v0::AbstractVector, 
    n::Integer; 
    tol::Real=1e-3
)
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
