"""
32 point groips
"""
const PointGroups = [
    "C1" , "Ci" , "C2" , "Cs" , "C2h", "D2" , "C2v", "D2h", 
    "C4" , "S4" , "C4h", "D4" , "C4v", "D2d", "D4h", "C3" , 
    "C3i", "D3" , "C3v", "D3d", "C6" , "C3h", "C6h", "D6" , 
    "C6v", "D3h", "D6h", "T"  , "Th" , "O"  , "Td" , "Oh" 
]

export pointgroup 
function pointgroup(i::Integer)
    data = readdlm("$(@__DIR__)/Tables/$(PointGroups[i]).dat", '\t', Int)
    multab = reshape(data, size(data, 1), size(data, 2))
    FiniteGroup(multab, PointGroups[i])
end

function pointgroup(g::String)
    data = readdlm("$(@__DIR__)/Tables/$(g).dat", '\t', Int)
    multab = reshape(data, size(data, 1), size(data, 2))
    FiniteGroup(multab, g)
end
