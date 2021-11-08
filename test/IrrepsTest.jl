include("$(@__DIR__)/../src/FiniteGroups.jl")
using .FiniteGroups, LinearAlgebra, Test



for i = 1:32
    g = pointgroup(i)
    println("\nGroup : $(name(g))")
    reps = irreps(g)
    for rep in reps
        display(rep)
    end
end