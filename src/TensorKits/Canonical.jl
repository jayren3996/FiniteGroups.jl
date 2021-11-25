#---------------------------------------------------------------------------------------------------
# Helper Funtions
#
# 1. mat_sqrt computes matrix square root.
# 2. canonical_gauging regauges (sometimes with trucation) the MPS.
# 3. vals_group groups values into different degenerate classes.
#---------------------------------------------------------------------------------------------------
function mat_sqrt(mat::Hermitian; ztol::Real=1e-20)
    vals_all, vecs_all = eigen(mat)
    pos = vals_all .> ztol
    vals_sqrt = sqrt.(vals_all[pos])
    vecs = vecs_all[:, pos]
    L = vecs * Diagonal(vals_sqrt)
    R = Diagonal(1 ./ vals_sqrt) * vecs'
    L, R
end
#---------------------------------------------------------------------------------------------------
function canonical_gauging(Γ, L, R)
    @tensor Γ2[:] := R[-1,1] * Γ[1,-2,2] * L[2,-3]
    Γ2
end
#---------------------------------------------------------------------------------------------------
function vals_group(vals::Vector; sorttol::Real=1e-3)
    pos = Vector{Vector{Int64}}(undef, 0)
    temp = zeros(Int64, 0)
    current_val = vals[1]
    for i=1:length(vals)
        if abs(vals[i]-current_val) .< sorttol
            push!(temp, i)
        else
            push!(pos, temp)
            temp = [i]
            current_val = vals[i]
        end
    end
    push!(pos, temp)
    pos
end

#---------------------------------------------------------------------------------------------------
# Right Canonical Form
#
# 1. The algorithm will tensor that is right-normalized.
# 2. Automatically trim null-space.
# 3. Is NOT guaranteed that the outputs are non-degenerate.
# 4. The algorithm is UNSTABLE if the transfer matrix has huge condition number.
#---------------------------------------------------------------------------------------------------
export right_canonical
function right_canonical(
    Γ::AbstractArray{<:Number, 3};
    itr::Integer=100,
    ztol::Real=1e-20
)
    ρ = steady_mat(Γ, itr)
    L, R = mat_sqrt(ρ, ztol=ztol)
    canonical_gauging(Γ, L, R)
end

#---------------------------------------------------------------------------------------------------
# Block Decomposition

# 1. Further decompose the right-renormalized tensor into right-renormalized tensors.
# 2. Ensure the outputs are non-degenerate.
# 3. block_trim only return one block when degeneracy is encountered.
#---------------------------------------------------------------------------------------------------
function block_decomp(
    Γ::AbstractArray{<:Number, 3};
    itr::Integer=100,
    sorttol::Real=1e-3
)
    vals, vecs = rand_fixed_mat(Γ, itr) |> eigen
    vgroup = vals_group(vals, sorttol=sorttol)
    res = begin
        num = length(vgroup)
        ctype = promote_type(eltype(Γ), eltype(vecs))
        Vector{Array{ctype, 3}}(undef, num)
    end
    for i=1:num
        p = vecs[:, vgroup[i]]
        Γc = canonical_gauging(Γ, p, p')
        res[i] = Γc
    end
    res
end
#---------------------------------------------------------------------------------------------------
function block_trim(
    Γ::AbstractArray{<:Number, 3};
    itr::Integer=100,
    sorttol::Real=1e-3
)
    vals, vecs = rand_fixed_mat(Γ, itr) |> eigen
    p = begin
        vgroup = vals_group(vals, sorttol=sorttol)
        i = length.(vgroup) |> argmin
        vecs[:, vgroup[i]]
    end
    canonical_gauging(Γ, p, p')
end

#---------------------------------------------------------------------------------------------------
# Block Right Canonical Form
#
# 1. All-at-once method.
# 2. Return multiple non-degenerate right-canonical form.
#---------------------------------------------------------------------------------------------------
export block_canonical
function block_canonical(
    Γ::AbstractArray{<:Number, 3};
    itr::Integer=100,
    ztol::Real=1e-20,
    sorttol::Real=1e-3,
    trim::Bool=false
)
    Γ_RC = right_canonical(Γ, itr=itr, ztol=ztol)
    res = if trim
        block_trim(Γ_RC, itr=itr, sorttol=sorttol)
    else
        block_decomp(Γ_RC, itr=itr, sorttol=sorttol)
    end
    res
end

#---------------------------------------------------------------------------------------------------
# Schmidt Canonical Form
#
# 1. Given a right canonical form, return a Schmidt canonical form.
# 2. This algorithm assume there is no degeneracy.
#---------------------------------------------------------------------------------------------------
export schmidt_canonical
function schmidt_canonical(
    Γ::AbstractArray{<:Number,3};
    itr::Integer=100,
    bound::Integer=100,
    cutoff::Real=1e-14,
    renormalize::Bool=true
)
    X, Yt = begin
        R = steady_mat(Γ, itr)
        L = steady_mat(Γ, itr, dir=:l)
        R_res = cholesky(R)
        L_res = cholesky(L)
        R_res.L, L_res.U
    end
    U, S, V = svd_trim(Yt * X, bound, cutoff, renormalize)
    R_mat = inv(Yt) * U * Diagonal(S)
    L_mat = V * inv(X)
    Γ_new = canonical_gauging(Γ, R_mat, L_mat)
    Γ_new, S
end
