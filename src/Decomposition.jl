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
function transform_mat(
    rep::AbstractVector{<:AbstractMatrix}, 
    irep::AbstractVector{<:AbstractMatrix};
    R::Bool=false
)
    D = size(rep[1], 1)
    d = size(irep[1], 1)
    M = sum(kron(conj.(irep[i]), rep[i]) for i=1:length(rep))
    e, v = eigen(M)
    @assert abs(e[end]-length(rep)) < 1e-7 "No such irrep, max value = $(e[end]/length(rep))"
    U = reshape(v[:, end], D, d) * sqrt(d)
    norm(imag(U)) < 1e-7 && return real(U)
    R && norm(imag(U)) > 1e-7 && return [real(U) imag(U)]/sqrt(2)
    U
end

#-------------------------------------------------------------------------------
function block_decomposition(rep::AbstractVector{<:AbstractMatrix}, R::Bool=true)
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
end