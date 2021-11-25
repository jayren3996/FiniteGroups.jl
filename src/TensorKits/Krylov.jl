#---------------------------------------------------------------------------------------------------
# Kraus operator
#---------------------------------------------------------------------------------------------------
import LinearAlgebra: mul!, adjoint
import Base: *, eltype, Array
export kraus
struct Kraus{TL, TU}
    KL::Array{TL, 3}
    KU::Array{TU, 3}
    dir::Symbol
end
eltype(::Kraus{TL, TU}) where TL where TU = promote_type(TL, TU)
adjoint(k::Kraus) = k.dir == :r ? Kraus(k.KL, k.KU, :l) : Kraus(k.KL, k.KU, :r)
kraus(KL, KR; dir::Symbol=:r) = Kraus(Array(KL), Array(KR), dir)
kraus(KL; dir::Symbol=:r) = Kraus(Array(KL), conj.(KL), dir)
function Array(k::Kraus)
    KL, KU, dir = k.KL, k.KU, k.dir
    α, β = size(KL, 1), size(KU, 1)
    dim = α * β
    mat = Matrix{eltype(k)}(undef, dim, dim)
    if dir == :r
        mat_r = reshape(mat, α, β, α, β)
        @tensor mat_r[:] = KL[-1,1,-3] * KU[-2,1,-4]
    else
        mat_r = reshape(mat, β, α, β, α)
        @tensor mat_r[:] = KU[-3,1,-1] * KL[-4,1,-2]
    end
    mat
end
#---------------------------------------------------------------------------------------------------
function mul!(ρ::Matrix, k::Kraus, ρ0::AbstractMatrix)
    KL, KU, dir = k.KL, k.KU, k.dir
    if dir == :r
        @tensor ρ[:] = KL[-1, 3, 1] * ρ0[1, 2] * KU[-2, 3, 2]
    elseif dir == :l
        @tensor ρ[:] = KU[1, 3, -1] * ρ0[1, 2] * KL[2, 3, -2]
    else
        error("Illegal direction: $dir.")
    end
end

function mul!(ρ::Vector, k::Kraus, ρ0::AbstractVector)
    α = size(k.KL, 1)
    ρ_r = reshape(ρ, α, α)
    ρ0_r = reshape(ρ0, α, α)
    mul!(ρ_r, k, ρ0_r)
end
#---------------------------------------------------------------------------------------------------
function *(k::Kraus, ρ0::AbstractVecOrMat)
    ctype = promote_type(eltype.((k, ρ0))...)
    ρ = Array{ctype}(undef, size(ρ0))
    mul!(ρ, k, ρ0)
    ρ
end

(k::Kraus)(v::AbstractVector) = k * v

#---------------------------------------------------------------------------------------------------
# Eigen system using power iteration

# Find dominent eigensystem by iterative multiplication.
# Krylov method ensures Hermicity and semi-positivity.
#---------------------------------------------------------------------------------------------------
function complex_iter!(K, ρ1, ρ2)
    mul!(ρ2, K, ρ1)
    mul!(ρ1, K, ρ2)
    normalize!(ρ1)
end
#---------------------------------------------------------------------------------------------------
function real_iter!(K, ρ1, ρ2)
    mul!(ρ2, K, ρ1)
    mul!(ρ1, K, ρ2)
    ρ1 .+= ρ2
    normalize!(ρ1)
end
#---------------------------------------------------------------------------------------------------
function power_iteration!(K, ρ1::Matrix, itr::Integer; method::Symbol=:r)
    ρ2 = similar(ρ1)
    if method == :r
        for i = 1:itr
            real_iter!(K, ρ1, ρ2)
        end
    elseif method == :c
        for i = 1:itr
            complex_iter!(K, ρ1, ρ2)
        end
    else
        error("Illegal method: $method.")
    end
    ρ1
end

#---------------------------------------------------------------------------------------------------
# Eigen system using Arnoldi algorithm

# Ard/noldi method ensures Hermicity but NOT semi-positivity.
#---------------------------------------------------------------------------------------------------
function krylov_arnoldi(K, ρ::Matrix)
    α = size(ρ, 1)
    v0 = reshape(ρ, :)
    e,v = eigsolve(K, v0, 1, :LR)
    reshape(v[1], α, α)
end
#---------------------------------------------------------------------------------------------------
# steady state from identity mat
function steady_mat(
    K::Array{<:Number, 3},
    itr::Integer=100;
    dir::Symbol=:r
)
    α = size(K, 3)
    Kc = conj(K)
    kraus = Kraus(K, Kc, dir)
    ρ = Array{eltype(K)}(I(α))
    mat = if α < 20
        power_iteration!(kraus, ρ, itr) |> Hermitian
    else
        m = krylov_arnoldi(kraus, ρ)
        if real(tr(m)) < 0.0
            m .= -1
        end
        m
    end
    Hermitian(mat)
end
#---------------------------------------------------------------------------------------------------
# Random fixed-point matrix.
function rand_fixed_mat(
    K::Array{<:Number, 3},
    itr::Integer=100;
    dir::Symbol=:r
)
    α = size(K, 3)
    Kc = conj(K)
    kraus = Kraus(K, Kc, dir)
    ρ = rand(ComplexF64, α, α) |> Hermitian |> Array
    mat = if α < 20
        power_iteration!(kraus, ρ, itr) |> Hermitian
    else
        krylov_arnoldi(kraus, ρ)
    end
    Hermitian(mat)
end
