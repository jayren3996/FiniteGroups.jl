include("$(@__DIR__)/../src/FiniteGroups.jl")
using .FiniteGroups, LinearAlgebra
import Base: *
struct Permutation{T <: Integer}
    M::Matrix{T}
    function Permutation(M::AbstractMatrix{<:Integer})
        row = [i for i = 1:size(M, 1) if M[i,1] ≠ M[i,2]]
        dup = [i for i = 2:length(row) if M[row[i-1], 1] == M[row[i], 1]]
        deleteat!(row, dup)
        new{eltype(M)}(Array(M[row, :]))
    end
end

function Base.string(p::Permutation)
    cyc = tocycles(p)
    str = "Cycle: "
    for c in cyc
        n = length(c)
        str *= "("
        for i = 1:n-1
            str *= "$(c[i]), "
        end
        str *= "$(c[n]))"
    end
    str
end
Base.repr(p::Permutation) = string(p)
Base.display(p::Permutation) = println(string(p))

Base.isequal(p1::Permutation, p2::Permutation) = isequal(p1.M, p2.M)
Base.isless(p1::Permutation, p2::Permutation) = isless(reshape(p1.M, :), reshape(p2.M, :))

function permutation(M::AbstractMatrix{<:Integer})
    p = sortperm(view(M, :, 1))
    M = M[p, 1:2]
    row = [i for i = 1:size(M, 1) if M[i,1] ≠ M[i,2]]
    dup = [i for i = 2:length(row) if M[row[i-1], 1] == M[row[i], 1]]
    deleteat!(row, dup)
    Permutation(Array(M[row, :]))
end

function cycles(cyc::Union{AbstractVector{<:Integer}, Tuple}...)
    n = sum(length(c) for c in cyc)
    iszero(n) && (return Permutation(zeros(Int64, 0, 2)))
    M = Matrix{Int64}(undef, n, 2)
    i = 1
    for c in cyc
        j = i + length(c) - 1
        M[i:j, 1] .= c
        M[j, 2] = c[1]
        M[i:j-1, 2] .= c[2:end]
        i = j + 1
    end
    permutation(M)
end

function tocycles(p::Permutation)
    cyc = []
    ind = p.M[:, 1]
    n = length(ind)
    used = fill(false, n)
    i = 1
    while i <= n
        used[i] ? (i += 1; continue) : (used[i] = true)
        j = ind[i]
        c = [j]
        k = p(j)
        while k ≠ j
            used[binary_search(ind, k)[2]] = true
            push!(c, k)
            k = p(k)
        end
        isone(length(c)) || push!(cyc, c)
    end
    cyc
end

function (p::Permutation)(i::Integer)
    q, j = binary_search(view(p.M, :, 1), i)
    q ? p.M[j, 2] : i
end

function *(p1::Permutation, p2::Permutation)
    i = sort([view(p1.M, :, 1); view(p2.M, :, 1)]) |> delete_duplicate!
    j = p1.(p2.(i))
    Permutation([i j])
end

function _get_position!(eles, i, j)
    ek = eles[i] * eles[j]
    k = findfirst(x -> isequal(x, ek), eles)
    if isnothing(k)
        push!(eles, ek)
        length(eles)
    else
        k
    end
end

function generate_group(gens::AbstractVector)
    ng = length(gens)
    eles = [gen for gen in gens]
    tb = [_get_position!(eles, i, j) for i=1:ng, j=1:ng]
    N1, N2 = ng, length(eles)
    while N1 < N2
        tb2 = [_get_position!(eles, i, j) for i=1:N1, j=N1+1:N2]
        tb3 = [_get_position!(eles, i, j) for i=N1+1:N2, j=1:N2]
        tb = [tb tb2; tb3]
        N1 = N2
        N2 = length(eles)
    end
    tb, eles
end

"""
    binary_search
"""
function binary_search(list::AbstractVector, i)
    l, r = 1, length(list)
    if iszero(r) || i < list[l]
        return false, 1
    elseif i > list[r]
        return false, r+1
    elseif i == list[l]
        return true, l
    elseif i == list[r]
        return true, r
    end
    while true
        c = (l + r) ÷ 2
        t = list[c]
        if i < t
            r = c
        elseif i > t
            l = c
        else
            return true, c
        end
        isone(r - l) && (return false, r)
    end
end

function delete_duplicate!(v)
    N = length(v)
    i = 2
    while i<=N
        if v[i] == v[i-1]
            deleteat!(v, i)
            N -= 1
        else
            i += 1
        end
    end
    v
end


e = cycles(())
a = cycles((1,2))
b = cycles((1,2,3,4,5,6))
res = generate_group([e,a,b]);
g = FiniteGroup(res[1])
ct = charactertable(g)
reg = regular_rep(g)
res = FiniteGroups.prep(g, ct[11,:], reg);
