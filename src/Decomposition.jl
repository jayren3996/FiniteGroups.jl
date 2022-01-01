#-------------------------------------------------------------------------------
# Decomposition
#-------------------------------------------------------------------------------
"""
Transform matrix from a (reducible) representation to a given irrep.

For a representation D(g) and irrep d(g), find S that
    S⁺⋅D(g)⋅S = d(g)

Take special note that the reshape function for row-major/column-major language:
    row_reshape(v, i,j) = col_reshape(v, j,i)ᵀ
"""
function proj_to_irreps(
    rep::AbstractVector{<:AbstractMatrix}, 
    irep::AbstractVector{<:AbstractMatrix};
    R::Bool=false
)
    D = size(irep[1], 1)
    P = proj_operator(rep, [conj(m[1,1]) for m in irep])
    e, v = eigen(Hermitian(P))
    vecs = v[:, e .> 1e-7]
    spaces = [krylov_space(rep, vecs[:, i], D) for i = 1:size(vecs, 2)]
    if R && !(promote_type(eltype.(irep)...) <: Real)
        [sqrt(2)*[real(s) imag(s)] for s in spaces] 
    else
        spaces
    end
end

#-------------------------------------------------------------------------------
function block_decomposition(
    rep::AbstractVector{<:AbstractMatrix};
    R::Bool=true
)
    n = length(rep)
    g = begin
        multab = Matrix{Integer}(undef, n, n)
        for i=1:n, j=1:n
            gij = rep[i]*rep[j]
            k = findfirst(x -> norm(x-gij)<1e-7, rep)
            multab[i,j] = k
        end
        FiniteGroup(multab)
    end
    ct = charactertable(g)
    row = single_complex_row(ct)
    reps = irreps(g, ct[row, :])
    vcat([proj_to_irreps(rep, irep, R=R) for irep in reps]...)
end
