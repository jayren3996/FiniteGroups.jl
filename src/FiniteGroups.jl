# Finite Group Theory Routines
module FiniteGroups
using LinearAlgebra, DelimitedFiles, SparseArrays


export name, order, class, inclass, mult
"""
    AbstractFiniteGroup

Abstract type of finite groups.
Group elements are labeled by integers.
The identity element are 1 by default.
"""
abstract type AbstractFiniteGroup end
name(g::AbstractFiniteGroup) = g.name
Base.getindex(g::AbstractFiniteGroup, i::Integer, j::Integer) = g.multab[i, j]
Base.inv(g::AbstractFiniteGroup, i) = g.inv[i]
class(g::AbstractFiniteGroup) = g.cls
class(g::AbstractFiniteGroup, i::Integer) = g.cls[i]
inclass(g::AbstractFiniteGroup, i::Integer) = g.clsv[i]
mult(g::AbstractFiniteGroup) = g.mult
mult(g::AbstractFiniteGroup, i::Integer) = g.mult[i]
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

function Base.display(g::AbstractFiniteGroup)
    println("Finite group : $(name(g))")
    println("Group order  : $(order(g))")
    println("Classes      : $(length(class(g)))")
end

include("FiniteGroup.jl")
include("Character.jl")
include("Dixon.jl")
include("Representation.jl")
include("PointGroups/PointGroups.jl")


end
