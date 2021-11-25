"""
32 point groips
"""
const PointGroups = [
    "C1" , "Ci" , "C2" , "Cs" , "C2h", "D2" , "C2v", "D2h", 
    "C4" , "S4" , "C4h", "D4" , "C4v", "D2d", "D4h", "C3" , 
    "C3i", "D3" , "C3v", "D3d", "C6" , "C3h", "C6h", "D6" , 
    "C6v", "D3h", "D6h", "T"  , "Th" , "O"  , "Td" , "Oh" 
]

struct PointGroup <: AbstractFiniteGroup
    name::String
    multab::Matrix{Int64}
    inv::Vector{Int64}
    cls::Vector{Vector{Int64}}
    clsv::Vector{Int64}
    mult::Vector{Int64}
    operations::Vector{SymOperation{3}}
end

function Base.display(g::PointGroup)
    println("Point group : $(name(g))")
    println("Group order : $(order(g))")
    println("Classes     : $(length(class(g)))")
end

export pointgroup
function pointgroup(i::Integer)
    g = Crystalline.pointgroup(i)
    multab = Crystalline.MultTable(g).table
    ginv = group_inverse(multab)
    cls, clsv = conjugate_class(multab, ginv)
    mult = length.(cls)
    PointGroup(
        PointGroups[i], 
        multab,
        ginv,
        cls,
        clsv,
        mult,
        g.operations
    )
end

pointgroup(s::String) = findfirst(x->x==s, PointGroups) |> pointgroup

export operation
operation(g::PointGroup) = g.operations
operation(g::PointGroup, i) = g.operations[i]
