# Representations

## Irreducible representations

[`irreps`](@ref) returns the irreducible matrix representations of a group: a
vector with one entry per irrep, each entry itself a vector of matrices (one per
group element). By default these are complex unitary representations.

```@example reps
using FiniteGroups

g = pointgroup("D3")
rs = irreps(g)
[size(r[1], 1) for r in rs]       # irrep dimensions — D₃ has irreps of dim 1, 1, 2
```

The matrices really are representations: `rs[k][i] * rs[k][j]` equals
`rs[k][m]` whenever ``g_i g_j = g_m`` in the group.

## Real and complex realizations

Whether an irrep can be written with real matrices is governed by its
Frobenius–Schur indicator: `+1` for a real representation, `-1` for a
pseudo-real one, and `0` for a genuinely complex one. [`check_real_rep`](@ref)
computes the indicator for a representation, [`real_rep`](@ref) realifies a
single representation, and [`real_irreps`](@ref) returns *all* the irreps in
real form:

```@example reps
ri = real_irreps(g)
[eltype(r[1]) for r in ri]        # real-valued matrices
```

## The regular representation and projectors

[`regular_rep`](@ref) builds the (sparse, integer) regular representation;
[`proj_operator`](@ref) forms the group-algebra projector ``\sum_i c_i\, r_i``
from a character; and [`block_decomposition`](@ref) / [`proj_to_irrep`](@ref)
decompose a given representation into its irreducible blocks. This is the
machinery `irreps` itself uses to split the regular representation:

```julia
reg    = regular_rep(g)                 # |G| sparse permutation matrices
blocks = block_decomposition(reg, g)    # split into irreducible blocks
```

## Projective and double-group representations

For spin-1/2 and double-group problems the package builds projective
representations, where group multiplication is realized only up to a phase
(factor system). [`cover_group`](@ref) forms the relevant central extension,
[`proj_reps`](@ref) computes the projective irreps for a given factor system,
[`chiral_proj_reps`](@ref) handles the chiral case, and
[`check_proj_coeff`](@ref) verifies a factor system. The chiral routine takes a
one-dimensional real irrep (a character with `χ[1] == 1`) as its phase factor:

```julia
ct = charactertable(g)
# use the one-dimensional A₂ irrep (row 2) as the chiral factor
rreps = chiral_proj_reps(g, Int.(ct[2, :]); R=true)
```

## Manipulating representations

A few utilities round out the toolkit: [`transform_rep`](@ref) conjugates a
representation by a change of basis, [`oplus`](@ref) forms the direct sum of
representations, [`unitary_rep`](@ref) returns an equivalent unitary
representation, and [`equivalent_transform`](@ref) finds the similarity
transformation relating two equivalent representations.
