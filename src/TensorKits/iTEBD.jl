#---------------------------------------------------------------------------------------------------
# iTEBD
#---------------------------------------------------------------------------------------------------
# gate
function gate!(
    mat::AbstractMatrix,
    ΓA::AbstractArray{<:Number, 3},
    ΓB::AbstractArray{<:Number, 3},
    λl::AbstractVector{<:Real},
    d1::Integer,
    d2::Integer,
    bound::Integer,
    cutoff::Real
)
    α, β = size(ΓA, 1), size(ΓB, 3)
    tlmul!(λl, ΓA)
    Γ = tcmul(mat, tgroup(ΓA, ΓB))
    ΓA2, λc, ΓB2 = tsvd(Γ, α, d1, d2, β, bound=bound, cutoff=cutoff, renormalize=true)
    trmul!(ΓA2, λc)
    tlmul!(1.0 ./ λl, ΓA2)
    ΓA2, λc, ΓB2
end
#---------------------------------------------------------------------------------------------------
# 2-site iTEBD
export itebd2!
function itebd2!(
    mat::AbstractMatrix,
    ΓA::AbstractArray{<:Number, 3},
    ΓB::AbstractArray{<:Number, 3},
    λl::AbstractVector{<:Real};
    bound::Integer=100,
    cutoff::Real=1e-7
)
    d1, d2 = size(ΓA, 2), size(ΓB, 2)
    ΓA2, λ1, ΓB2 = gate!(mat,  ΓA,  ΓB, λl, d1, d2, bound, cutoff)
    ΓB3, λ2, ΓA3 = gate!(mat, ΓB2, ΓA2, λ1, d2, d1, bound, cutoff)
    ΓA3, ΓB3, λ2
end

#---------------------------------------------------------------------------------------------------
# 3-site iTEBD
export itebd3!
function itebd3!(
    mat::AbstractMatrix,
    ΓAB1::AbstractArray{<:Number, 3},
    ΓC1::AbstractArray{<:Number, 3},
    λl::AbstractVector{<:Real};
    bound::Integer=100,
    cutoff::Real=1e-7
)
    d1, d2, d3 = size(ΓA, 2), size(ΓB, 2),  size(ΓC, 2)
    ΓA2, λ1, ΓBC2 = gate!(mat, ΓAB1, ΓC1, λl, d1, d2 * d3, bound, cutoff)
    ΓB3, λ2, ΓCA3 = gate!(mat, ΓBC2, ΓA2, λ1, d2, d3 * d1, bound, cutoff)
    ΓC4, λ3, ΓAB4 = gate!(mat, ΓBC2, ΓA2, λ2, d3, d1 * d2, bound, cutoff)
    ΓAB4, ΓC4, λ3
end
