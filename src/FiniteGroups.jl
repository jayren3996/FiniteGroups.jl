# Finite Group Theory Routines
module FiniteGroups

using LinearAlgebra, SparseArrays, Base.Threads
import Base.:*
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
include("PointGroups.jl")
include("Character.jl")
include("Representation.jl")
include("ProjReps.jl")
include("Permutation.jl")
include("Decomposition.jl")

# Precompilation
C6h = pointgroup(23)
C6h_ct = charactertable(C6h)
irreps(C6h, C6h_ct)
irreps(C6h, R=true)

Oh = pointgroup("Oh")
irreps(Oh)


end
