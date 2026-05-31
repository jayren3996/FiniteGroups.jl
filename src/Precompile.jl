# Precompilation workload.
#
# Run a small but representative set of calls during precompilation only
# (`jl_generating_output() == 1`) so the actual hot signatures — including the
# sparse-integer regular-representation path used by `irreps` — are baked into
# the precompile cache. This replaces a hand-maintained `precompile(...)` list
# whose concrete signatures had drifted from the real call sites (e.g. it
# precompiled dense `proj_operator` while the hot call is sparse-integer).
if ccall(:jl_generating_output, Cint, ()) == 1
    let
        # Point-group path: bundled data, Burnside table, projected irreps.
        gp = pointgroup("D3")
        ctp = charactertable(gp)
        irreps(gp, ctp)
        irreps(gp; R=true)

        # Generic FiniteGroup path: Burnside from a multiplication table.
        gf = FiniteGroup([1 2 3; 2 3 1; 3 1 2], "C3")
        irreps(gf)

        # Permutation-group path.
        gs = permutationgroup(3)
        irreps(gs)
    end
end
