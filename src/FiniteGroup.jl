struct FiniteGroup{T <: Integer} <: AbstractFiniteGroup
    name::String
    multab::Matrix{T}
    inv::Vector{T}
    cls::Vector{Vector{T}}
    clsv::Vector{T}
    mult::Vector{Int64}
end

export finitegroup
function finitegroup(
    multab::AbstractMatrix{<:Integer}, 
    name::String="Unnamed group"
)
    ginv = group_inverse(multab)
    cls, clsv = conjugate_class(multab, ginv)
    mult = length.(cls)
    FiniteGroup(name, multab, ginv, cls, clsv, mult)
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