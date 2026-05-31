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

function Base.show(io::IO, ::MIME"text/plain", g::AbstractFiniteGroup)
    println(io, "Finite group : $(name(g))")
    println(io, "Group order  : $(order(g))")
    print(io,   "Classes      : $(length(class(g)))")
end
Base.show(io::IO, g::AbstractFiniteGroup) = print(io, nameof(typeof(g)), "(\"", name(g), "\", order ", order(g), ")")

include("FiniteGroup.jl")
include("Character.jl")
include("Representation.jl")
include("ProjReps.jl")
include("Decomposition.jl")
include("Dixon.jl")
include("PointGroups.jl")
include("Permutation.jl")
include("Precompile.jl")


end
