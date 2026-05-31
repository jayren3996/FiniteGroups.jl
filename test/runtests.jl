using FiniteGroups, LinearAlgebra, Test

@testset "FiniteGroups.jl" begin
    include("BurnsideTest.jl")
    include("DecompTest.jl")
    include("ProjRepTest.jl")
    include("DixonTest.jl")
    include("ReviewFixesTest.jl")
end
