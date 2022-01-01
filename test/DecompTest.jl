include("$(@__DIR__)/../src/FiniteGroups.jl")
using .FiniteGroups, LinearAlgebra, Test
import .FiniteGroups: proj_to_irreps, block_decomposition

g = pointgroup("Oh")
rreps = irreps(g, R=true)

rep = oplus(
    rreps[1], 
    rreps[2], rreps[2], rreps[2], rreps[2], 
    rreps[5], 
    rreps[6], 
    rreps[8], rreps[8], 
    rreps[9], 
    rreps[10]
)
ru = Symmetric(rand(size(rep[1])...)) |> eigvecs
rep = transform_rep(ru', rep, ru)
res = block_decomposition(rep)
