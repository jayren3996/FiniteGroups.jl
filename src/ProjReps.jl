export cover_group
function cover_group(
    g::AbstractFiniteGroup, 
    coeff::AbstractMatrix{<:Integer}, 
    p::Integer
)
    n = order(g)
    mt = Matrix{Int64}(undef, n * p, n * p)
    for i = 1:n, j = 1:n
        k = g[i, j]
        c = coeff[i, j]
        for a = 0:p-1, b = 0:p-1
            row = i + n * a
            col = j + n * b
            val = k + n * mod(a + b + c, p)
            mt[row, col] = val
        end
    end
    FiniteGroup(mt, "Cover Group for $(name(g))")
end

export check_proj_coeff
function check_proj_coeff(
    g::AbstractFiniteGroup, 
    r::AbstractVector{<:AbstractMatrix},
    coeff::AbstractMatrix{<:Integer}, 
    p::Integer;
    tol::Real=1e-7
)
    n = order(g)
    for i = 1:n, j = 1:n
        k = g[i,j]
        err = r[i] * r[j] - exp(1im * 2π/p * coeff[i, j]) * r[k]
        norm(err) > tol && (return false)
    end
    true
end

export proj_reps
function proj_reps(
    g::AbstractFiniteGroup, 
    coeff::AbstractMatrix{<:Integer}, 
    p::Integer;
    tol::Real=1e-7
)
    cg = cover_group(g, coeff, p)
    reps = irreps(cg)
    check = [check_proj_coeff(g, rep, coeff, p, tol=tol) for rep in reps]
    [reps[i][1:order(g)] for i=1:length(reps) if check[i]]
end