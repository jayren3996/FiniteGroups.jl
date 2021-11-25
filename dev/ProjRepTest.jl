include("$(@__DIR__)/../src/FiniteGroups.jl")
using .FiniteGroups, LinearAlgebra, Test

g = pointgroup("Oh")

function get_coeff(g::FiniteGroup, χ)
    n = order(g)
    ct = Matrix{Int64}(undef, 2n, 2n)
    ct[1:2n,1:n] .= 0
    for i = 1:n
        p = χ[inclass(g, i)]
        ct[i, n+1:2n] .= (1-p)÷2
        ct[n+i, n+1:2n] .= (1-p)÷2
    end
    ct
end

function bigger_group(g::FiniteGroup)
    n = order(g)
    mt = Matrix{Int64}(undef, 2n, 2n)
    for i=1:n, j=1:n
        k = g[i,j]
        mt[i, j] = k
        mt[n+i, j] = n+k
        mt[i, n+j] = n+k
        mt[n+i, n+j] = k
    end
    FiniteGroup(mt)
end

function check_group(g::FiniteGroup)
    n = order(g)
    g.multab[1, :] = 1:n
    g.multab[:, 1] = 1:n
    for i=1:n, j=1:n, k=1:n
        g[g[i,j],k] == g[i,g[j,k]] || error("$i, $j, $k not associate")
    end
    for i = 1:n
        sort(g.multab[i, :]) == 1:n || error("Not cycle")
    end
    true
end

ctab = charactertable(g); display(ctab)
χ = ctab[4,:];
coeff = get_coeff(g, χ);
cg = cover_group(bigger_group(g), coeff, 2);
reps = irreps(cg);
check = [check_proj_coeff(bigger_group(g), rep, coeff, 2, tol=1e-7) for rep in reps]

ct = charactertable(cg)
reg = regular_rep(cg)
p = proj_operator(ct[])
