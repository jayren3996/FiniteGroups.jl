precompile(FiniteGroup, (Matrix{Int64}, String))
precompile(group_inverse, (Matrix{Int64},))
precompile(conjugate_class, (Matrix{Int64}, Vector{Int64}))
precompile(collect_group, (Vector{Int64}, Int64))
precompile(beautify, (Float64,))
precompile(beautify, (ComplexF64,))
precompile(class_multab, (FiniteGroup,))
precompile(class_multab, (PermutationGroup,))
precompile(spectrum_split, (Vector{Float64},))
precompile(spectrum_split, (Vector{ComplexF64,}))
precompile(vec_split, (Matrix{Float64}, Matrix{Float64}))
precompile(vec_split, (Matrix{Float64}, Matrix{ComplexF64}))
precompile(vec_split, (Matrix{ComplexF64}, Matrix{Float64}))
precompile(vec_split, (Matrix{ComplexF64}, Matrix{ComplexF64}))
precompile(normalize_chi!, (FiniteGroup, Vector{Float64}))
precompile(normalize_chi!, (FiniteGroup, Vector{ComplexF64}))
precompile(normalize_chi!, (PermutationGroup, Vector{Float64}))
precompile(normalize_chi!, (PermutationGroup, Vector{ComplexF64}))
precompile(burnside, (FiniteGroup, Array{Int64, 3}))
precompile(charactertable, (FiniteGroup,))
precompile(charactertable, (PointGroup,))
precompile(charactertable, (PermutationGroup,))
precompile(regular_rep, (FiniteGroup,))
precompile(regular_rep, (PointGroup,))
precompile(regular_rep, (PermutationGroup,))
precompile(unitaryeigen, (Matrix{Float64},))
precompile(unitaryeigen, (Matrix{ComplexF64},))
precompile(unitarysqrt, (Matrix{Float64},))
precompile(unitarysqrt, (Matrix{ComplexF64},))
precompile(transform_rep, (Matrix{Float64}, Vector{Matrix{Float64}}, Matrix{Float64}))
precompile(transform_rep, (Matrix{Float64}, Vector{Matrix{ComplexF64}}, Matrix{Float64}))
precompile(transform_rep, (Matrix{ComplexF64}, Vector{Matrix{ComplexF64}}, Matrix{ComplexF64}))
precompile(transform_rep, (Matrix{ComplexF64}, Vector{Matrix{Float64}}, Matrix{ComplexF64}))
precompile(krylov_space, (Vector{Matrix{Float64}}, Vector{Float64}, Int64))
precompile(krylov_space, (Vector{Matrix{ComplexF64}}, Vector{ComplexF64}, Int64))
precompile(find_least_degen, (Vector{Matrix{Float64}}, Int64))
precompile(find_least_degen, (Vector{Matrix{ComplexF64}}, Int64))