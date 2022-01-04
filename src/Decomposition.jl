#-------------------------------------------------------------------------------
# Decomposition
#-------------------------------------------------------------------------------
export proj_to_irrep
"""
Transform matrix from a (reducible) representation to a given irrep.

For a representation D(g) and irrep d(g), find S that
    S⁺⋅D(g)⋅S = d(g)

Take special note that the reshape function for row-major/column-major language:
    row_reshape(v, i,j) = col_reshape(v, j,i)ᵀ
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
    if R && !(promote_type(eltype.(irep)...) <: Real)
        [sqrt(2)*[real(s) imag(s)] for s in spaces] 
    else
        spaces
    end
end

#-------------------------------------------------------------------------------
function block_decomposition(
    rep::AbstractVector{<:AbstractMatrix};
    g::Union{AbstractFiniteGroup, Nothing}=nothing,
    R::Bool=true
)
    group = if isnothing(g)
        n = length(rep)
        multab = Matrix{Integer}(undef, n, n)
        for i=1:n, j=1:n
            gij = rep[i]*rep[j]
            k = findfirst(x -> norm(x-gij)<1e-7, rep)
            multab[i,j] = k
        end
        FiniteGroup(multab)
    else
        g
    end
    reps = begin
        ct = charactertable(group)
        row = single_complex_row(ct)
        irreps(group, ct[row, :])
    end
    vcat([proj_to_irrep(rep, irep, R=R) for irep in reps]...)
end
