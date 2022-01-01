export proj_reps, cover_group, check_proj_coeff
"""
    proj_reps(g, coeff, p; R, tol)

Calculate the projective representation of `g` with coefficients `coeff`.

Inputs:
-------
g     : Finite group object.
coeff : Coefficients represented by an integer matrix.
p     : Modulus of the the integer.
R     : Real or not, default = `false`.
tol   : Tolerance, default = `1e-7`.
"""
function proj_reps(
    g::AbstractFiniteGroup, 
    coeff::AbstractMatrix{<:Integer}, 
    p::Integer;
    R::Bool=false,
    tol::Real=1e-7
)
    cg = cover_group(g, coeff, p)
    reps = irreps(cg, R=R)
    check = [check_proj_coeff(g, rep, coeff, p, tol=tol) for rep in reps]
    [reps[i][1:order(g)] for i=1:length(reps) if check[i]]
end
#-------------------------------------------------------------------------------
"""
    cover_group(g::AbstractFiniteGroup, coeff::AbstractMatrix{<:Integer}, p::Integer)

Construct the cover group of `g`, given the projective coefficient `coeff`. 
The projective coefficient is unitary and the phase is module p.

Inputs:
-------
g     : Finite group.
coeff : D(gᵢ) * D(gⱼ) = exp(i2π/p * coeff[i,j]) * D(gᵢ⋅gⱼ).
p     : Module of coefficient.
"""
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
#-------------------------------------------------------------------------------
"""
Check whether a linear representation of the cover group of `g` is indeed the 
projective representation obeying the projective coefficients.
"""
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

#-------------------------------------------------------------------------------
# Chiral projective representation.
#-------------------------------------------------------------------------------
export chiral_proj_reps
"""

"""
function chiral_proj_reps(
    g::AbstractFiniteGroup, 
    χ::AbstractVector{<:Integer};
    R::Bool=true
)
    n = order(g)
    bg, coeff = double_group(g, χ)
    reps = proj_reps(bg, coeff, 2, R=R)
    out = []
    for i = 1:length(reps)
        rep = reps[i]
        S = rep[n+1]
        e, v = eigen(Symmetric(S))
        dim = length(e)÷2
        @assert norm(e[1:dim] .+ 1) < 1e-7 "Chiral operator not unitary, eigvals = $e."
        @assert norm(e[dim+1:2dim] .- 1) < 1e-7 "Chiral operator not unitary, eigvals = $e."
        U = hcat(svd(v[:,dim+1:2dim]).U, svd(v[:, 1:dim]).U)
        @assert size(U) == (2dim, 2dim) "Dimension net correct, got $(size(U)), expect $((2dim, 2dim))"
        Ui = U'
        push!(out, [Ui * rep[i] * U for i = 1:n])
    end
    out
end
#-------------------------------------------------------------------------------
"""
    double_group(g::AbstractFiniteGroup, χ::AbstractVector{<:Integer})

Construct a double group with element g⊕S⋅g, satisfying:
    gᵢS⋅gⱼ = S⋅gᵢgⱼ,
The projective coefficients is:
    D(gᵢ)D(S⋅gⱼ) = χᵢD(S⋅gᵢgⱼ)

The function takes Finitegroup `g` and coefficent `χ`, and return the doubled 
group together with the projective coefficients. 
"""
function double_group(g::AbstractFiniteGroup, χ::AbstractVector{<:Integer})
    @assert length(χ) == length(g.cls) "Length of χ should be equal to classes of g, while length(χ) = $(length(χ))."
    n = order(g)
    coeff = begin
        ct = zeros(Int64, 2n, 2n)
        for i = 1:n
            p = χ[inclass(g, i)]
            if p == -1
                ct[i, n+1:2n] .= 1
                ct[n+i, n+1:2n] .= 1
            elseif p != 1
                error("p should be ±1, get $p.")
            end
        end
        ct
    end
    mt = Matrix{Int64}(undef, 2n, 2n)
    for i=1:n, j=1:n
        k = g[i,j]
        mt[i, j] = k
        mt[n+i, j] = n+k
        mt[i, n+j] = n+k
        mt[n+i, n+j] = k
    end
    FiniteGroup(mt), coeff
end