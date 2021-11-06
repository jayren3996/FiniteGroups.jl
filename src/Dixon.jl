#-------------------------------------------------------------------------------
# Dixon's Mothod
# In developing
#-------------------------------------------------------------------------------
export dixon
function dixon(
    g::AbstractFiniteGroup, 
    h::AbstractArray{<:Integer, 3}
)
    n, p = findprime(g)
    Z = primitive_root(n, p)
    ξ = exp(1im * 2π/n)
    NC = size(h, 1)
    v = [I(NC)]
    for i = 1:NC
        hi = h[i, :, :]
        v = vcat([vec_split(hi, vi, p) for vi in v]...)
        all(isone(size(vi, 2)) for vi in v) && break
    end
    mat = vcat(transpose.(v)...)
    for i = 1:size(mat, 1)
        normalize_modp!(view(mat, i, :), p)
    end
    maptocomplex.(mat, Z, n, ξ)
end

function eigvecs_modp(m::AbstractMatrix{<:Integer}, modp::Integer)
    vecs = []
    for λ = 0:modp-1
        sol = eigvecs_modp(m, λ, modp)
        if size(sol, 2) > 0
            push!(vecs, sol)
        end
    end
    @assert sum(size(vec, 2) for vec in vecs) == size(m, 1) "Solution imcomplete."
    vecs
end

function eigvecs_modp(m::AbstractMatrix{<:Integer}, λ::Integer, modp::Integer)
    n = size(m, 1)
    herm = hermiteform!([m-λ*I I], modp)
    R = 1
    while R <= n
        iszero(sum(herm[R, 1:n])) && break
        R += 1
    end
    vecs = herm[R:end, n+1:end]
    #for i = 1:size(vecs, 1)
    #    nornmalize_modp!(view(vecs, i, :), modp)
    #end
    Array(vecs')
end

function normalize_modp!(v::AbstractVector{<:Integer}, modp::Integer)
    v0 = copy(v)
    for i = 1:modp
        v .+= v0
        mod_vector!(v, modp)
        isone(v[1]) && return v
    end
    error("Non-renormalizable.")
end

function vec_split(
    h::AbstractMatrix{<:Integer}, 
    vs::AbstractVecOrMat{<:Integer},
    modp::Integer
)
    isone(size(vs, 2)) && return [vs]
    php = vs' * h * vs 
    v = eigvecs_modp(php, modp)
    [vs * vi for vi in v]
end

#-------------------------------------------------------------------------------
# Hermitian normal form
#-------------------------------------------------------------------------------
"""
    hermiteform!(m::AbstractMatrix{<:Integer}, modp::Integer)

Compute the Hermite form of a matrix in the integer field Fₚ 
"""
function hermiteform!(m::AbstractMatrix{<:Integer}, modp::Integer)
    pointer = 1
    x, y = size(m)
    rank = size(m, 1)
    for i = 1:x
        while pointer <= y && set_prime_row!(m, i, col=pointer) 
            pointer += 1
        end
        if pointer > y
            rank = i-1
            break
        end
        for j = i+1:x
            row_reduce!(view(m, i, :), view(m, j, :), modp, col=pointer)
        end
        pointer += 1
    end
    for i = 2:rank
        n = findfirst(x -> x != 0, @view m[i, :])
        if n === nothing
            break
        else
            for j = 1:i-1
                kill_row!(view(m, i, :), view(m, j, :), modp, col=n)
            end
        end
        if n == y
            break
        end
    end
    m[1:rank, :]
end


function row_reduce!(
    v1::AbstractVector{<:Integer}, 
    v2::AbstractVector{<:Integer}, 
    modp::Integer; 
    col::Integer=1
)
    @assert length(v1) == length(v2) "v1, v2 should be of the same size."
    if abs(v1[col]) < abs(v2[col])
        if v1[col] == 0
            (v2[col] > 0) || (v2 .*= -1)
            mod_vector!(v2, modp)
            swap_vector!(v1, v2)
        else
            kill_row!(v1, v2, modp, col=col)
            row_reduce!(v1, v2, modp, col=col)
        end
    else
        if v2[col] == 0
            return nothing
        else
            kill_row!(v2, v1, modp, col=col)
            row_reduce!(v1, v2, modp, col=col)
        end
    end
end

function kill_row!(
    v1::AbstractVector{<:Integer}, 
    v2::AbstractVector{<:Integer}, 
    modp::Integer; 
    col::Integer=1
)
    n, r = divrem(v2[col], v1[col])
    v2 .-= n * v1
    mod_vector!(v2, modp)
end

function set_prime_row!(
    m::AbstractMatrix{<:Integer}, 
    row::Integer; 
    col::Integer=1
)
    if m[row, col] < 0
        m[row, :] .*= -1
    elseif m[row, col] == 0
        n = findfirst(x -> x != 0, @view m[row+1:end, col])
        if n === nothing
            return true
        else
            swap_vector!(view(m, row, :), view(m, row+n, :), col=col)
        end
    end
    return false
end

function swap_vector!(
    v1::AbstractVector, 
    v2::AbstractVector; 
    col::Integer=1
)
    @assert length(v1) == length(v2) "v1, v2 should be of the same size."
    if v2[col] < 0
        for i = 1:length(v1)
            v1[i], v2[i] = -v2[i], v1[i]
        end
    else
        for i = 1:length(v1)
            v1[i], v2[i] = v2[i], v1[i]
        end
    end
end

function mod_vector!(v::AbstractVector{<:Integer}, modp::Integer)
    for i = 1:length(v)
        v[i] = mod(v[i], modp)
    end
end

# Integer field calculation Fₚ
function primitive_root(n::Integer, modp::Integer)
    isone(n) && (return n)
    for i = 2:modp 
        p = i
        flag = true
        for j = 1:n-2
            p = mod(i * p, modp)
            isone(p) && (flag = false; break)
        end
        flag && isone(mod(i * p, modp)) && (return i)
    end
end

function findprime(g::AbstractFiniteGroup)
    n = lcm([order(g, c[1]) for c in class(g)])
    primes = [2, 3, 5, 7, 11, 13, 15, 17, 19, 23, 29, 31]
    i = ceil(Int64, 2 * sqrt(order(g)))
    i = i^2 > 4 * order(g) ? min(32, i) : min(32, i+1)
    while true
        Q = true
        for p in primes
            iszero(mod(i, p)) && (Q=false; break)
            (p^2 > i) && break
        end
        if Q
            iszero(mod(i-1, n)) ? (return (n, i)) : push!(primes, i)
        end
        i += 1
    end
end

function maptocomplex(N::Integer, Z::Integer, n::Integer, ξ::Number)
    list = zeros(Int, n)
    for i = 1:n-1
        N, list[i] = divrem(N, Z)
    end
    list[n] = N
    reverse!(list)
    C = list[1]
    for i = 2:n
        C = ξ * C + list[i]
    end
    C
end