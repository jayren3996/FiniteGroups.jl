struct Representation{G<:AbstractFiniteGroup, T<:Number, M<:AbstractVector{<:AbstractMatrix}}
    g::G
    m::M
    χ::Vector{T}
    reduced::Bool
end

function Representation(
    g::AbstractFiniteGroup,
    m::AbstractVector{<:AbstractMatrix};
    χ::Union{Nothing, AbstractVector{<:Number}}=nothing,
    reduced::Union{Nothing, Bool}=nothing
)
    vχ = if isnothing(χ)
        [tr(m[i]) for i = 1:length(m)]
    else
        Array(χ)
    end
    isnothing(reduced) && ( reduced = abs(sum(χ)-order(g)) < 1e-7 )
    Representation(g, m, vχ, reduced)
end

function Base.display(r::Representation)
    println("$(size(r.m[1], 1))-d $(r.reduced ? "irreducible" : "reducible") representation of $(name(r.g)): ")
    for i = 1:length(r.m)
        println("\nElement $i:")
        display(r.m[i])
    end
end

Base.getindex(r::Representation, i) = r.m[i]
Base.length(r::Representation) = length(r.m)
Base.iterate(r::Representation) = (r[1], 1)
Base.iterate(r::Representation, i::Integer) = i == length(r.m) ? nothing : (r[i+1], i+1)


@inline character(r::Representation) = r.χ

export regular_rep
"""
    regular_rep(multab::AbstractMatrix{<:Integer})

Return regular representation:
"""
function regular_rep(g::AbstractFiniteGroup)
    n = order(g)
    m = [sparse([g[i,j] for j=1:n], 1:n, fill(1, n), n, n) for i=1:n]
    reduced = order(g) == 1
    Representation(g, m, reduced=reduced)
end

export irreps
function irreps(
    g::AbstractFiniteGroup;
    χ::Union{Nothing, AbstractMatrix}=nothing,
    method::String="Burnside"
)
    isnothing(χ) && ( χ = charactertable(g, method=method) )
    reg = regular_rep(g)
    get_irreps(g, χ, reg)
end

"""
    get_irreps(g::FiniteGroup, χ::AbstractMatrix, reg::AbstractArray{<:Integer, 3})

Get the irreducible representations.
"""
function get_irreps(
    g::AbstractFiniteGroup, 
    χ::AbstractMatrix, 
    reg::Representation
)
    irreps = [prep(g, χ[i, :], reg) for i = 1:size(χ, 2)]
    irreps
end


function prep(
    g::AbstractFiniteGroup, 
    χ::AbstractVector{<:Number}, 
    reg::Representation;
    tol::Real=1e-7
)
    D = round(Int64, real(χ[1]))
    preg = begin
        vs = if order(g) > 3000
            eigs(proj_operator(χ, reg), nev=D^2)[2]
        else
            e, v = eigen(proj_operator(χ, reg))
            v[:, e .> 1e-2]
        end
        @assert size(vs, 2) == D^2 "Dimension error: exprect $(D^2), get $(size(vs, 2))"
        [vs' * regmat * vs for regmat in reg]
    end
    isone(D) && return Representation(g, [beautify.(mi) for mi in preg], χ=χ, reduced=true)
    v0 = nothing
    for i = 1:length(class(g))
        rep = preg[class(g,i)[1]]
        e, v = eigen(rep)
        e_split = spectrum_split(e, tol=tol)
        for ei in e_split
            if length(ei) == D
                v0 = v[:,ei[1]]
                #println(i)
                break
            end
        end
        isnothing(v0) || break
    end
    m = if isnothing(v0) 
        # No nondegenerate eigenvector
        min_degen = order(g)
        element = nothing
        for i = 1:length(class(g))
            rep = preg[class(g,i)[1]]
            e, v = eigen(rep)
            e_split = spectrum_split(e, tol=tol)
            degen = minimum(length.(e_split))
            if degen < min_degen
                min_degen = degen
                element = class(g,i)[1]
            end
        end
        error("No nondegenerate eigenvector. Minimum degeneracy: $min_degen for $element'th element.")
    else
        vs = krylov_space(preg[2:order(g)], v0, D)
        [vs' * pregmat * vs for pregmat in preg]
    end
    Representation(g, m, χ=χ, reduced=true)
end

export proj_operator
function proj_operator(r::Representation)
    D = size(r[1],1)
    m = zeros(eltype(r.χ) <: Real ? Float64 : ComplexF64, D, D)
    for i = 1:order(r.g)
        m .+= conj(r.χ[inclass(r.g, i)]) * r[i]
    end
    Hermitian(m)
end
function proj_operator(χ::AbstractVector{<:Number}, r::Representation)
    m = sum(conj(χ[inclass(r.g, i)]) * r[i] for i=1:order(r.g))
    Hermitian(Array(m))
end

function krylov_space(
    mats::AbstractVector{<:AbstractMatrix}, 
    v0::AbstractVector, 
    n::Integer; 
    tol::Real=1e-3
)
    r, i, nmat = 1, 0, length(mats)
    vs = v0
    while r < n
        i = mod(i, nmat) + 1
        vs = hcat(vs, mats[i] * vs)
        U, S = svd(vs)
        pos = S .> tol
        vs = U[:, pos]
        r = size(vs, 2)
    end
    @assert r == n "Dimension of Krylov space not match: expect d = $n, get d = $r"
    vs
end
