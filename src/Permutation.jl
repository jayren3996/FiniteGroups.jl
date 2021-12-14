struct Permutation{T <: Integer}
    M::Matrix{T}
    function Permutation(M::AbstractMatrix{<:Integer})
        row = [i for i = 1:size(M, 1) if M[i,1] ≠ M[i,2]]
        dup = [i for i = 2:length(row) if M[row[i-1], 1] == M[row[i], 1]]
        deleteat!(row, dup)
        new{eltype(M)}(Array(M[row, :]))
    end
end

struct PermutationGroup <: AbstractFiniteGroup
    name::String
    multab::Matrix{Int64}
    inv::Vector{Int64}
    cls::Vector{Vector{Int64}}
    clsv::Vector{Int64}
    mult::Vector{Int64}
    operations::Vector{Permutation}
end

function Base.display(g::PointGroup)
    println("Permutation group : $(name(g))")
    println("Group order       : $(order(g))")
    println("Classes           : $(length(class(g)))")
end

operation(g::PermutationGroup) = g.operations
operation(g::PermutationGroup, i) = g.operations[i]

export permutationgroup
function permutationgroup(gens::AbstractVector{<:Permutation}; name::String="Unnamed Group")
    multab, eles = generate_group(gens)
    ginv = group_inverse(multab)
    cls, clsv = conjugate_class(multab, ginv)
    mult = length.(cls)
    PermutationGroup(name, multab, ginv, cls, clsv, mult, eles)
end

permutationgroup(n::Integer) = permutationgroup([cycles(1:2), cycles(1:n)], name="S$n")


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

export permutation
function permutation(M::AbstractMatrix{<:Integer})
    p = sortperm(view(M, :, 1))
    M = M[p, 1:2]
    row = [i for i = 1:size(M, 1) if M[i,1] ≠ M[i,2]]
    dup = [i for i = 2:length(row) if M[row[i-1], 1] == M[row[i], 1]]
    deleteat!(row, dup)
    Permutation(Array(M[row, :]))
end

export cycles
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

export permutation
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

function generate_group(gens::AbstractVector)
    N = 0
    eles = gens 
    while length(eles) > N
        N = length(eles)
        eles = vcat(eles, ([gen * ele for ele in eles] for gen in gens)...) |> sort
        delete_duplicate!(eles)
    end
    tb = Matrix{Int64}(undef, N, N)
    @threads for i = 1:N 
        for j=1:N
            _, k = binary_search(eles, eles[i] * eles[j])
            tb[i, j] = k
        end
    end
    tb, eles
end

#-------------------------------------------------------------------------------
# Helpers
#-------------------------------------------------------------------------------
"""
    binary_search
"""
function binary_search(list::AbstractVector, i)
    l, r = 1, length(list)
    if iszero(r) || i < list[l]
        return false, 1
    elseif i > list[r]
        return false, r+1
    elseif isequal(i, list[l])
        return true, l
    elseif isequal(i, list[r])
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
        if isequal(v[i], v[i-1])
            deleteat!(v, i)
            N -= 1
        else
            i += 1
        end
    end
    v
end
