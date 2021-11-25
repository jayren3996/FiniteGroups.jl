#---------------------------------------------------------------------------------------------------
# Low-level multiplication
#---------------------------------------------------------------------------------------------------
export tlmul!, trmul!, tcmul!, tcmul
function tlmul!(
    λ::AbstractVector{<:Number}, 
    Γ::AbstractArray
)
    α = size(Γ, 1)
    Γ_r = reshape(Γ, α, :)
    lmul!(Diagonal(λ), Γ_r)
end
#---------------------------------------------------------------------------------------------------
function trmul!(
    Γ::AbstractArray, 
    λ::AbstractVector{<:Number}
)
    β = size(Γ)[end]
    Γ_r = reshape(Γ, :, β)
    rmul!(Γ_r, Diagonal(λ))
end
#---------------------------------------------------------------------------------------------------
function tcmul!(
    Γ_n::AbstractArray{<:Number, 3}, 
    cmat::AbstractMatrix, 
    Γ::AbstractArray{<:Number, 3}
)
    @tensor Γ_n[:] = cmat[-2,1] * Γ[-1,1,-3]
end
#---------------------------------------------------------------------------------------------------
function tcmul(
    cmat::AbstractMatrix, 
    Γ::AbstractArray{<:Number, 3}
)
    @tensor Γ_n[:] := cmat[-2,1] * Γ[-1,1,-3]
    Γ_n
end

#---------------------------------------------------------------------------------------------------
# Grouping
#---------------------------------------------------------------------------------------------------
export tgroup, tgroup!
function tgroup!(
    Γ::AbstractArray{<:Number, 3}, 
    ΓA::AbstractArray{<:Number, 3}, 
    ΓB::AbstractArray{<:Number, 3}
)
    α, β = size(ΓA, 1), size(ΓB, 3)
    d1, d2 = size(ΓA, 2), size(ΓB, 2)
    ΓA_r = reshape(ΓA, α * d1, :)
    ΓB_r = reshape(ΓB, :, d2 * β)
    Γ_r = reshape(Γ, α * d1, d2 * β)
    mul!(Γ_r, ΓA_r, ΓB_r)
end
#---------------------------------------------------------------------------------------------------
function tgroup(
    ΓA::AbstractArray{<:Number, 3}, 
    ΓB::AbstractArray{<:Number, 3}
)
    α, β = size(ΓA, 1), size(ΓB, 3)
    d1, d2 = size(ΓA, 2), size(ΓB, 2)
    ΓA_r = reshape(ΓA, α * d1, :)
    ΓB_r = reshape(ΓB, :, d2 * β)
    reshape(ΓA_r * ΓB_r, α, d1 * d2, β)
end
#---------------------------------------------------------------------------------------------------
function tgroup!(
    Γ::AbstractArray{<:Number, 3}, 
    ΓA::AbstractArray{<:Number, 3}, 
    ΓB::AbstractArray{<:Number, 3},
    ΓC::AbstractArray{<:Number, 3}
)
    α, β = size(ΓA, 1), size(ΓC, 3)
    d1, d2, d3 = size(ΓA, 2), size(ΓB, 2), size(ΓC, 2)
    Γ_r = reshape(Γ, α, d1, d2, d3, β)
    @tensor Γ_r[:] = ΓA[-1,-2,1] * ΓB[1,-3,2] * ΓC[2,-4,-5]
end
#---------------------------------------------------------------------------------------------------
function tgroup(
    ΓA::AbstractArray{<:Number, 3}, 
    ΓB::AbstractArray{<:Number, 3},
    ΓC::AbstractArray{<:Number, 3}
)
    @tensor Γ_r[:] := ΓA[-1,-2,1] * ΓB[1,-3,2] * ΓC[2,-4,-5]
    reshape(Γ_r, size(ΓA, 1), :, size(ΓB, 3))
end

#---------------------------------------------------------------------------------------------------
# Tensor SVD
#---------------------------------------------------------------------------------------------------
export svd_trim, tsvd
function svd_trim(
    mat::AbstractMatrix,
    bound::Integer,
    cutoff::Real,
    renormalize::Bool
)
    res = svd(mat)
    vals = res.S
    len = min(count(i -> (i > cutoff), vals), bound)
    U = res.U[:, 1:len]
    S = vals[1:len]
    V = res.Vt[1:len, :]
    if renormalize
        S ./= norm(S)
    end
    U, S, V
end
#---------------------------------------------------------------------------------------------------
function tsvd(
    T::AbstractArray,
    α::Integer,
    d1::Integer, 
    d2::Integer, 
    β::Integer;
    bound::Integer=1000,
    cutoff::Real=1e-7,
    renormalize::Bool=false
)
    T_r = reshape(T, α * d1, d2 * β)
    U, S, V = svd_trim(T_r, bound, cutoff, renormalize)
    ΓA = reshape(U, α, d1, :)
    ΓB = reshape(V, :, d2, β)
    ΓA, S, ΓB
end