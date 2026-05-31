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

function Base.show(io::IO, ::MIME"text/plain", g::PermutationGroup)
    println(io, "Permutation group : $(name(g))")
    println(io, "Group order       : $(order(g))")
    print(io,   "Classes           : $(length(class(g)))")
end

operation(g::PermutationGroup) = g.operations
operation(g::PermutationGroup, i) = g.operations[i]

export permutationgroup
"""
    permutationgroup(gens::AbstractVector{<:Permutation}; name="Unnamed Group")
    permutationgroup(n::Integer)

Build a permutation group: either close a set of generating permutations (see
[`cycles`](@ref) and [`permutation`](@ref)) into a full group, or, given an
integer `n`, construct the symmetric group `Sₙ`.

# Examples
```julia
julia> permutationgroup(4)          # S₄
Permutation group : S4
Group order       : 24
Classes           : 5
```
"""
function permutationgroup(gens::AbstractVector{<:Permutation}; name::String="Unnamed Group")
    multab, eles = generate_group(gens)
    ginv = group_inverse(multab)
    cls, clsv = conjugate_class(multab, ginv)
    mult = length.(cls)
    PermutationGroup(name, multab, ginv, cls, clsv, mult, eles)
end

function permutationgroup(n::Integer)
    n ≥ 0 || throw(ArgumentError("permutationgroup(n) needs n ≥ 0; got n = $n."))
    # S₀ and S₁ are trivial: a transposition (1 2) makes sense only for n ≥ 2,
    # otherwise the generating set below would spuriously build an order-2 group.
    n ≥ 2 || return permutationgroup([Permutation(zeros(Int, 0, 2))], name="S$n")
    permutationgroup([cycles(1:2), cycles(1:n)], name="S$n")
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
Base.show(io::IO, p::Permutation) = print(io, string(p))

Base.isequal(p1::Permutation, p2::Permutation) = isequal(p1.M, p2.M)
# Equal permutations share the same canonical `M`, so `==`/`hash` track `isequal`
# and make `Permutation`s behave correctly as `Set`/`Dict` keys.
==(p1::Permutation, p2::Permutation) = p1.M == p2.M
Base.hash(p::Permutation, h::UInt) = hash(p.M, hash(:Permutation, h))
Base.isless(p1::Permutation, p2::Permutation) = isless(reshape(p1.M, :), reshape(p2.M, :))

export permutation
"""
    permutation(M::AbstractMatrix{<:Integer})

Construct a permutation from a two-column matrix `M`: each row `[a b]` means
`a ↦ b`, and fixed points may be omitted. See [`cycles`](@ref) for
disjoint-cycle notation.
"""
function permutation(M::AbstractMatrix{<:Integer})
    size(M, 2) == 2 || throw(ArgumentError("permutation expects a two-column matrix [a b] (a ↦ b); got size $(size(M))."))
    p = sortperm(view(M, :, 1))
    M = M[p, 1:2]
    row = [i for i = 1:size(M, 1) if M[i,1] ≠ M[i,2]]   # ignore fixed points a ↦ a
    src = M[row, 1]                                       # sorted ascending by construction
    dst = M[row, 2]
    # A genuine permutation maps its moved points bijectively onto themselves:
    # sources must be distinct and equal, as a set, to the targets. Without this
    # an invalid map (e.g. `[1 2]`, i.e. 1 ↦ 2 with nothing mapping to 1) is built
    # and later loops forever in `tocycles`/display.
    (allunique(src) && sort(dst) == src) ||
        throw(ArgumentError("not a valid permutation: $(Array(M[row, :])) (rows [a b] mean a ↦ b) is not a bijection on its moved points."))
    Permutation(Array(M[row, :]))
end

export cycles
"""
    cycles(cyc...)

Construct a permutation from disjoint cycles, each given as a vector or tuple of
point labels; e.g. `cycles((1,2,3), (4,5))` is the permutation `(1 2 3)(4 5)`.
With no cycles the identity permutation is returned. See [`permutation`](@ref)
for the two-column mapping form and [`permutationgroup`](@ref) to close a set of
generators into a group.
"""
function cycles(cyc::Union{AbstractVector{<:Integer}, Tuple}...)
    n = sum(length(c) for c in cyc; init=0)   # init=0 so `cycles()` is the identity
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
