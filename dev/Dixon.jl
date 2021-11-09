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
        normalize_modp!(g, view(mat, i, :), p)
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
    herm = hermiteform!([m-λ*I; I], modp)
    C = 1
    while C <= n
        iszero(sum(herm[1:n, C])) && break
        C += 1
    end
    herm[n+1:end, C:end]
end

function vec_split(
    h::AbstractMatrix{<:Integer}, 
    vs::AbstractVecOrMat{<:Integer},
    modp::Integer
)
    isone(size(vs, 2)) && return [vs]
    php = vs' * h * vs 
    v = eigvecs_modp(php, modp)
    [mod!(vs * vi, modp) for vi in v]
end

function normalize_modp!(g::AbstractFiniteGroup, v::AbstractVector{<:Integer}, modp::Integer)
    v0 = copy(v)
    for i = 1:modp
        v .+= v0
        mod!(v, modp)
        isone(v[1]) && break
    end
    d2 = order(g) // sum(v[i] * inv_modp(v[i], modp) // mult(g, i) for i =1:length(v))
    d2 = mod(Int(d2), modp)
    d = round(Int, sqrt(d2))
    @assert d^2 == d2
    v .*= d
    mod!(v, modp)
end

#-------------------------------------------------------------------------------
# Calculation on prime integer field Fₚ
#-------------------------------------------------------------------------------
function inv_modp(i::Integer, p::Integer)
    for j = 2:p-1
        isone(mod(i*j, p)) && (return j)
    end
    error("Did not find inverse.")
end

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
            d, r = divrem(i, p)
            iszero(r) && (Q=false; break)
            (d <= p) && break
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

#-------------------------------------------------------------------------------
# Hermitian normal form on prime integer field Fₚ
#-------------------------------------------------------------------------------
"""
    hermiteform!(m::AbstractMatrix{<:Integer}, modp::Integer)

Compute the Hermite form of a matrix in the prime integer field Fₚ 
"""
function hermiteform!(m::AbstractMatrix{<:Integer}, modp::Integer)
    mod!(m, modp)
    pointer = 1
    x, y = size(m)
    for i = 1:y-1
        while pointer <= x && set_prime_column!(m, i, row=pointer) 
            pointer += 1
        end
        pointer > x && break
        for j = i+1:y
            column_reduce!(view(m, :, i), view(m, :, j), modp, row=pointer)
        end
        pointer += 1
    end
    m
end

"""
    column_reduce!(v1, v2, modp; col)

Reduce column `v2` using column `v1`, with prime row set to `row`.
"""
function column_reduce!(
    v1::AbstractVector{<:Integer}, 
    v2::AbstractVector{<:Integer}, 
    modp::Integer; 
    row::Integer=1
)
    iszero(v2[row]) && (return nothing)
    v1[row] > v2[row] && swap!(v1, v2)
    kill_column!(v1, v2, modp, row=row)
    column_reduce!(v1, v2, modp, row=row)
end

"""
Basic column manipulation.
"""
function kill_column!(
    v1::AbstractVector{<:Integer}, 
    v2::AbstractVector{<:Integer}, 
    modp::Integer; 
    row::Integer=1
)
    v2 .-= (v2[row] ÷ v1[row]) * v1
    mod!(v2, modp)
end

"""
    set_prime_column!(m, row; col)

Swap with the remaining column so that m[row, col] > 0.
Return true if its impossible, and false otherwise.
"""
function set_prime_column!(
    m::AbstractMatrix{<:Integer}, 
    col::Integer; 
    row::Integer=1
)
    iszero(m[row, col]) || (return false)
    for j = col+1:size(m, 2)
        if !iszero(m[row, j])
            swap!(view(m, :, col), view(m, :, j))
            return false
        end
    end
    true
end

"""
In place swap of `v1` and `v2`.
"""
@inline function swap!(v1::AbstractArray, v2::AbstractArray)
    for i = 1:length(v1)
        v1[i], v2[i] = v2[i], v1[i]
    end
end

"""
In place mod p on `v`.
"""
@inline function mod!(v::AbstractArray, modp::Integer)
    for i = 1:length(v)
        v[i] = mod(v[i], modp)
    end
    v
end