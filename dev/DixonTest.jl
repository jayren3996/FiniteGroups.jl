include("$(@__DIR__)/../src/FiniteGroups.jl")
using .FiniteGroups, LinearAlgebra, Test
import .FiniteGroups: findprime, primitive_root, vec_split, normalize_modp!, maptocomplex, eigvecs_modp
import .FiniteGroups: hermiteform!, column_reduce!, kill_column!, mod!, set_prime_column!
g = pointgroup(32)
h = class_multab(g)
n, p = findprime(g)
Z = primitive_root(n, p)
ξ = exp(1im * 2π/n)
NC = size(h, 1)
v = [I(NC)]

for i = 1:NC
    hi = h[i, :, :]
    v = vcat([vec_split(hi, vi, p) for vi in v]...)
    all(isone(size(vi, 2)) for vi in v) && break
end
"""
mat = vcat(transpose.(v)...)
for i = 1:size(mat, 1)
    normalize_modp!(g, view(mat, i, :), p)
end
maptocomplex.(mat, Z, n, ξ)
"""