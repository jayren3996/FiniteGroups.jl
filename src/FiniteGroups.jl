# Finite Group Theory Routines
module FiniteGroups

using LinearAlgebra, SparseArrays, Base.Threads
import Base: ==, *
import LinearAlgebra: dot
export name, order, class, inclass, mult

"""
    AbstractFiniteGroup

Abstract type of finite groups.
Group elements are labeled by integers.
The identity element are 1 by default.
"""
abstract type AbstractFiniteGroup end

function Base.display(g::AbstractFiniteGroup)
    println("Finite group : $(name(g))")
    println("Group order  : $(order(g))")
    println("Classes      : $(length(class(g)))")
end

include("FiniteGroup.jl")
include("Character.jl")
include("Representation.jl")
include("ProjReps.jl")
include("Decomposition.jl")
include("PointGroups.jl")
include("Permutation.jl")
include("Precompile.jl")


end
