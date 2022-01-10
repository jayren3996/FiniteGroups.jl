#-------------------------------------------------------------------------------
# Decomposition
#-------------------------------------------------------------------------------
export proj_to_irrep
"""
Transform matrix from a (reducible) representation to a given irrep.

For a representation D(g) and irrep d(g), find S that
    S⁺⋅D(g)⋅S = d(g)
"""
function proj_to_irrep(
    rep::AbstractVector{<:AbstractMatrix}, 
    irep::AbstractVector{<:AbstractMatrix};
    R::Bool=false
)
    N, D = length(rep), size(irep[1], 1)
    vecs = map(1:D) do i 
        Pii = proj_operator(rep, [conj(m[i,i]) for m in irep])
        e, v = eigen(Hermitian(Pii))
        v[:, e .> 1e-7]
    end
    for i = 2:D
        Pi1 = proj_operator(rep, [conj(m[i,1]) for m in irep])
        U = vecs[i]' * Pi1 * vecs[1] * D / N
        #@assert norm(U*U'-I) < 1e-7 "Not unitary."
        vecs[i] = vecs[i] * U
    end
    spaces = map(1:size(vecs[1], 2)) do i 
        hcat((vecs[j][:, i] for j=1:D)...)
    end
    if R && sum(norm.(imag.(rep))) > 1e-7
        [sqrt(2)*[real(s) imag(s)] for s in spaces] 
    else
        spaces
    end
end
#-------------------------------------------------------------------------------
"""
Find a list of projection Sᵢ that
    Sᵢ⁺⋅D(g)⋅Sᵢ = dᵢ(g),
where {dᵢ} is a list of irreducible representations as an input.
"""
function proj_to_irrep(
    rep::AbstractVector{<:AbstractMatrix}, 
    ireps::AbstractVector{<:AbstractVector{<:AbstractMatrix}};
    R::Bool=false
)
    vcat([proj_to_irrep(rep, irep, R=R) for irep in ireps]...)
end
#-------------------------------------------------------------------------------
function block_decomposition(
    rep::AbstractVector{<:AbstractMatrix},
    group::AbstractFiniteGroup;
    R::Bool=true
)
    ct = charactertable(group)
    row = single_complex_row(ct)
    ireps =  irreps(group, ct[row, :])
    proj_to_irrep(rep, ireps, R=R)
end
#-------------------------------------------------------------------------------
function block_decomposition(rep::AbstractVector{<:AbstractMatrix}; R::Bool=true)
    n = length(rep)
    multab = Matrix{Integer}(undef, n, n)
    for i=1:n, j=1:n
        gij = rep[i]*rep[j]
        gij_index = 0
        for k = 1:n
            if norm(rep[k] - gij) < 1e-10
                @assert iszero(gij_index) "Not faithful representation. An explicit `FiniteGroup` should be provided."
                gij_index = k
            end
        end
        @assert gij_index > 0 "The input `rep` is not a representation."
        multab[i,j] = gij_index
    end
    group = FiniteGroup(multab)
    block_decomposition(rep, group, R=R)
end
