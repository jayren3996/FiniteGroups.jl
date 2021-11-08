export charactertable
"""
    charactertable(g::AbstractFiniteGroup; method)

Return the character table of group `g`.
Use Burnside's method by defualt.
"""
function charactertable(g::AbstractFiniteGroup; method::String="Burnside")
    h = class_multab(g)
    if method == "Burnside"
        burnside(g, h)
    end
end

#-------------------------------------------------------------------------------
# Burnside's Mothod
#-------------------------------------------------------------------------------
"""
    burnside(g::FiniteGroup, h::AbstractArray{<:Integer, 3}; tol::Real=1e-7)

Get character table using Burnside's method.
"""
function burnside(
    g::AbstractFiniteGroup, 
    h::AbstractArray{<:Integer, 3}; 
    tol::Real=1e-7
)
    NC = size(h, 1)
    v = [I(NC)]
    while true
        rh = sum(rand() * h[i,:,:] for i = 1:NC)
        v = vcat([vec_split(rh, vi, tol=tol) for vi in v]...)
        all(isone(size(vi, 2)) for vi in v) && break
    end
    for i = 1:NC
        normalize_chi!(g, v[i])
    end
    χ = sort_chi(vcat(transpose.(v)...))
    beautify.(χ, tol=tol)
end

function normalize_chi!(
    g::AbstractFiniteGroup, 
    χ::AbstractArray
)
    NC = length(χ)
    χ ./= χ[1]
    d2 = order(g) / sum( abs2(χ[j])/mult(g,j) for j=1:NC )
    for i = 1:NC
        χ[i] *= sqrt(d2)/mult(g,i)
    end
    χ
end

function sort_chi(table)
    inds = [i for i = 1:size(table,1)]
    for i = 1:length(inds), j=i+1:length(inds)
        if !order_chi(table[inds[i], :], table[inds[j], :])
            temp = inds[i]
            inds[i] = inds[j]
            inds[j] = temp
        end
    end
    table[inds, :]
end

function order_chi(v1::AbstractVector, v2::AbstractVector; tol::Real=1e-7)
    for i = 1:length(v1)
        a, b = log(Complex(v1[i])), log(Complex(v2[i]))
        dr = real(a) - real(b)
        if dr > tol
            return false
        elseif dr < -tol
            return true
        else
            ai = imag(a) .< -tol ? 2π + imag(a) : imag(a)
            bi = imag(b) .< -tol ? 2π + imag(b) : imag(b)
            di = ai - bi
            if di > tol
                return false
            elseif di < -tol
                return true
            end
        end
    end
    println("Warning: two vector the same.")
    false
end

function vec_split(
    h::AbstractMatrix, 
    vs::AbstractVecOrMat; 
    tol::Real=1e-7
)
    isone(size(vs, 2)) && return [vs]
    php = vs' * h * vs 
    e, v = eigen(php)
    e_split = spectrum_split(e, tol=tol)
    # zero-energy degeneracy is numerically unstable!!!
    isone(length(e_split)) ? [vs] : [vs * v[:, ei] for ei in e_split]
end

"""
Group the spectrum into degenerate sets.
"""
function spectrum_split(vals::AbstractVector{<:Number}; tol::Real=1e-7)
    n = length(vals)
    NG = 0
    gvec = zeros(Int, n)
    for i = 1:n
        if iszero(gvec[i])
            NG += 1
            gvec[i] = NG
        else
            continue
        end
        target = vals[i]
        # Pay attention to the sorting of complex number:
        # a - b im, a + b im, ã - b̃ im, ã + b̃ im
        tally = 0
        for j = i+1:n
            if abs(vals[j] - target) < tol
                gvec[j] = NG
                tally = 0
            else
                tally += 1
            end
            (tally > 1) && break
        end
    end
    collect_group(gvec, NG)
end

#-------------------------------------------------------------------------------
# Helper Functions
#-------------------------------------------------------------------------------
export class_multab
"""
    class_multab(multable, classes, class_vec)

Compute class multiplicity:

Cᵢ × Cⱼ = hᵢⱼₖ Cₖ

Outputs:
--------
h::Array{Int64, 3} : Class multiplicity.
"""
function class_multab(g::AbstractFiniteGroup)
    NC = class(g) |> length
    h = zeros(Int64, NC, NC, NC)
    for i = 1:NC
        Ci = class(g, i)
        for a in Ci, b in Ci
            k = inclass(g, g[a, b])
            h[i, i, k] += 1
        end
        for j = i+1:NC
            Cj = class(g, j)
            for a in Ci, b in Cj 
                k = inclass(g, g[a, b])
                h[i, j, k] += 1
                h[j, i, k] += 1
            end
        end
    end
    for k=1:NC
        h[:, :, k] .÷= mult(g, k)
    end
    h
end

beautify(a::Integer; tol) = a
function beautify(a::Real; tol=1e-7)
    a_int = round(Int64, a)
    abs(a-a_int) < tol ? a_int : a
end
function beautify(a::Number; tol=1e-7)
    if abs(imag(a)) < tol
        beautify(real(a))
    elseif abs(real(a)) < tol
        1im * beautify(imag(a))
    else
        a
    end
end