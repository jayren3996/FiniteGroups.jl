export FiniteGroup
"""
    FiniteGroup(multab::AbstractMatrix{<:Integer}, name="Unnamed group")

A finite group defined by its multiplication table `multab`, where `multab[i, j]`
is the index of the product of elements `i` and `j` (element `1` is the identity).
Element inverses, conjugacy classes and class sizes are computed on construction.
Validate an untrusted table first with [`check_group`](@ref).

# Examples
```julia
julia> FiniteGroup([1 2 3; 2 3 1; 3 1 2], "C3")   # cyclic group of order 3
Finite group : C3
Group order  : 3
Classes      : 3
```
"""
struct FiniteGroup{T <: Integer} <: AbstractFiniteGroup
    name::String
    multab::Matrix{T}
    inv::Vector{T}
    cls::Vector{Vector{T}}
    clsv::Vector{T}
    mult::Vector{Int64}
end

function FiniteGroup(
    multab::AbstractMatrix{<:Integer}, 
    name::String="Unnamed group"
)
    ginv = group_inverse(multab)
    cls, clsv = conjugate_class(multab, ginv)
    mult = length.(cls)
    FiniteGroup(name, multab, ginv, cls, clsv, mult)
end

"""
    name(g) -> String

Name of group `g`; `name(g, i)` is the label of element/operation `i`.
"""
name(g::AbstractFiniteGroup) = g.name
name(::AbstractFiniteGroup, i::Integer) = string(i)

Base.getindex(g::AbstractFiniteGroup, i, j) = g.multab[i, j]

"""
    inv(g, i)

Index of the inverse of element `i` in group `g`.
"""
Base.inv(g::AbstractFiniteGroup, i) = g.inv[i]

"""
    class(g)    -> Vector{Vector{Int}}
    class(g, i) -> Vector{Int}

Conjugacy classes of `g`: all classes, or the elements of class `i`.
See also [`inclass`](@ref) and [`mult`](@ref).
"""
class(g::AbstractFiniteGroup) = g.cls
class(g::AbstractFiniteGroup, i) = g.cls[i]

"""
    inclass(g, i)

Index of the conjugacy class containing element `i`.
"""
inclass(g::AbstractFiniteGroup, i) = g.clsv[i]

"""
    mult(g)    -> Vector{Int}
    mult(g, i) -> Int

Size of each conjugacy class, or of class `i` (the class multiplicity). This is
*not* element multiplication — use `g[i, j]` for the product of elements `i`, `j`.
"""
mult(g::AbstractFiniteGroup) = g.mult
mult(g::AbstractFiniteGroup, i) = g.mult[i]

"""
    order(g)    -> Int
    order(g, i) -> Int

Order of the group `g` (number of elements), or the order of element `i`
(the smallest `k > 0` with `iᵏ = 1`).
"""
order(g::AbstractFiniteGroup) = length(g.inv)
function order(g::AbstractFiniteGroup, i::Integer)
    ord = 1
    gi = i
    while gi ≠ 1
        ord += 1
        gi = g[i, gi]
    end
    ord
end

"""
    group_inverse(multable::AbstractMatrix{<:Integer})

Compute the inverse of the group elements.

Outputs:
--------
ginv : Vector of indices.
"""
function group_inverse(multable::AbstractMatrix{<:Integer})
    n = size(multable, 1)
    ginv = Vector{eltype(multable)}(undef, n)
    ginv[1] = 1
    for j = 2:n, i = j:n
        if isone(multable[i, j])
            ginv[j] = i
            ginv[i] = j
        end
    end
    ginv
end

"""
    conjugate_class(multable::AbstractMatrix, ginv::Union{AbstractVector,Nothing}=nothing)

Compute the conjugate class.

Outputs:
--------
classes   : Vector of vector which contains the element in that class.
class_vec : Vector indicating the class each group element is in.
"""
function conjugate_class(
    multable::AbstractMatrix{<:Integer}, 
    ginv::Union{AbstractVector{<:Integer},Nothing}=nothing
)
    isnothing(ginv) && (ginv = group_inverse(multable))
    n = size(multable, 1)
    class_vec = zeros(eltype(multable), n)
    class_vec[1] = 1
    NC = 1
    for i = 2:n
        iszero(class_vec[i]) || continue
        NC += 1
        for j = 2:n
            k = multable[j, multable[i, ginv[j]]]
            class_vec[k] = NC
        end
    end
    classes = collect_group(class_vec, NC)
    classes, class_vec
end

"""
Given a vector specifying the classes each element is in, collect the elements in each class.
"""
function collect_group(vec::AbstractVector{<:Integer}, NG::Integer)
    groups = [eltype(vec)[] for i = 1:NG]
    for i = 1:length(vec)
        push!(groups[vec[i]], i)
    end
    groups
end

export check_group
"""
    check_group(multab) -> Bool

Validate that `multab` is a genuine group multiplication table: element `1` is the
identity, the operation is associative, and every row and column is a permutation.
Throws an informative error on the first violation; returns `true` otherwise.
"""
function check_group(multab)
    n = size(multab, 1)
    multab[1, :] == 1:n || error("Element 1 not identity")
    multab[:, 1] == 1:n || error("Element 1 not identity")
    for i=1:n, j=1:n, k=1:n
        multab[multab[i,j],k] == multab[i,multab[j,k]] || error("Elements $i, $j, $k not associative")
    end
    for i = 1:n
        sort(multab[i, :]) == 1:n || error("Not a permutation")
    end
    true
end

# Precompile
precompile(FiniteGroup, (Matrix{Int64}, String))
precompile(group_inverse, (Matrix{Int64},))
precompile(conjugate_class, (Matrix{Int64}, Vector{Int64}))
precompile(collect_group, (Vector{Int64}, Int64))
