#-------------------------------------------------------------------------------
# Type: Characters
#-------------------------------------------------------------------------------
struct Characters{T <: Number, Tg <: AbstractFiniteGroup}
    χ::AbstractVector{T}
    R::Int
    g::Tg
end

function Characters(g::AbstractFiniteGroup, χ::AbstractVector{<:Number})
    R = check_real_rep(g, χ)
    Characters(isone(R) ? real.(χ) : χ, R, g)
end
#-------------------------------------------------------------------------------
function Base.display(c::Characters)
    reality = if c.R == 1
        "Real"
    elseif c.R == 0 Complex
        "Complex"
    else
        "Pseudo-real"
    end
    println("Characters of $reality representation of $(name(c.g)):")
    println(c.χ)
end
#-------------------------------------------------------------------------------
Base.getindex(c::Characters, i) = c.χ[i]
Base.length(c::Characters) = length(c.χ)
Base.eltype(c::Characters{T,Tg}) where {T,Tg} = T
Base.conj(c::Characters) = Characters(conj(c.χ), c.R, c.g)
function Base.isless(c1::Characters, c2::Characters)
    @assert isequal(c1.g, c2.g) "Compared characters should belong to tha same group."
    order_chi(c1.χ, c2.χ)
end
function ==(c1::Characters, c2::Characters)
    isequal(c1.g, c2.g) && norm(c1.χ - c2.χ) < 1e-10
end
function dot(c1::Characters, c2::Characters)
    @assert isequal(c1.g, c2.g) "Dotted characters should belong to tha same group."
    g = c1.g
    p = sum(conj(c1.χ[k]) * c2.χ[k] * mult(g, k) for k = 1:length(class(g))) / order(g)
    pint = round(Int64, p)
    @assert abs(p-pint) < 1e-10 "Inner product = $p, not integer."
    pint
end
#-------------------------------------------------------------------------------
# Type: CharacterTable
#-------------------------------------------------------------------------------
struct CharacterTable 
    tab::Vector{Characters}
    names::Vector{String}
end

function CharacterTable(tab::Vector{Characters})
    names = Vector{String}(undef, length(tab))
    n = 0
    for i = 1:length(tab)
        isassigned(names, i) && continue
        n += 1
        if iszero(tab[i].R)
            χc = conj(tab[i])
            for j = i+1:length(tab)
                if tab[j] == χc
                    names[i] = "$(i)+"
                    names[j] = "$(j)-"
                    break
                end
            end
        else
            names[i] = "$(i)"
        end
    end
    CharacterTable(tab, names)
end
#-------------------------------------------------------------------------------
function Base.display(ct::CharacterTable)
    g = ct.tab[1].g
    mats = Matrix(ct)
    n = size(mats, 1)+1
    out = Matrix{Any}(undef, n, n)
    out[1, 1] = ""
    out[1, 2:end] .= [name(g, cls[1]) for cls in class(g)]
    out[2:end, 1] .= ct.names
    out[2:end, 2:end] .= beautify.(mats)
    display(out)
end
Base.Matrix(ct::CharacterTable) = vcat((transpose(c.χ) for c in ct.tab)...)
Base.length(ct::CharacterTable) = length(ct.tab)
Base.size(ct::CharacterTable) = (length(ct.tab), length(ct.tab[1]))
Base.size(ct::CharacterTable, i::Integer) = isone(i) ? length(ct.tab) : (i == 2) ? length(ct.tab[1]) : 1
Base.getindex(ct::CharacterTable, i::Integer, j) = ct.tab[i][j]
Base.getindex(ct::CharacterTable, inds::AbstractVector{<:Integer}, j) = vcat((transpose(ct.tab[i][j]) for i in inds)...)
Base.getindex(ct::CharacterTable, c::Colon, j) = vcat((transpose(ct.tab[i][j]) for i in 1:length(ct))...)
Base.getindex(ct::CharacterTable, i) = ct.tab[i]
#-------------------------------------------------------------------------------
# Main output function
#-------------------------------------------------------------------------------
export charactertable
"""
    charactertable(g::AbstractFiniteGroup; method)

Return the character table of group `g`.
Use Burnside's method.
"""
function charactertable(g::AbstractFiniteGroup; tol::Real=1e-7)
    h = class_multab(g)
    burnside(g, h, tol=tol)
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
    characters = Vector{Characters}(undef, NC)
    @threads for i = 1:NC
        χ = beautify.(reshape(normalize_chi!(g, v[i]), :), tol=tol)
        characters[i] = Characters(g, χ)
    end
    tab = sort(characters)
    CharacterTable(tab)
end
#-------------------------------------------------------------------------------
function normalize_chi!(g::AbstractFiniteGroup, χ::AbstractArray)
    NC = length(χ)
    χ ./= χ[1]
    d2 = order(g) / sum( abs2(χ[j])/mult(g,j) for j=1:NC )
    for i = 1:NC
        χ[i] *= sqrt(d2)/mult(g,i)
    end
    χ
end
#-------------------------------------------------------------------------------
function order_chi(v1::AbstractVector, v2::AbstractVector; tol::Real=1e-7)
    for i = 1:length(v1)
        a, b = log(Complex(v1[i])), log(Complex(v2[i]))
        (real(a) < -10) && (real(b) < -10) && continue
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
#-------------------------------------------------------------------------------
"""
Split the vector space according to the split of the eigen values.
"""
function vec_split(h::AbstractMatrix, vs::AbstractVecOrMat; tol::Real=1e-7)
    isone(size(vs, 2)) && return [vs]
    php = vs' * h * vs 
    e, v = eigen(php)
    e_split = spectrum_split(e, tol=tol)
    # zero-energy degeneracy is numerically unstable!!!
    isone(length(e_split)) ? [vs] : [vs * v[:, ei] for ei in e_split]
end
#-------------------------------------------------------------------------------
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
        #tally = 0
        for j = i+1:n
            if abs(vals[j] - target) < tol
                gvec[j] = NG
            #    tally = 0
            #else
            #    tally += 1
            end
            #(tally > 2) && break
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
#-------------------------------------------------------------------------------

function beautify(a::Real; tol::Real=1e-7)
    a_int = round(a)
    abs(a-a_int) < tol ? a_int : a
end
function beautify(a::Number; tol::Real=1e-7)
    if abs(imag(a)) < tol
        convert(typeof(a), beautify(real(a), tol=tol))
    elseif abs(real(a)) < tol
        convert(typeof(a), 1im * beautify(imag(a), tol=tol))
    else
        a
    end
end

