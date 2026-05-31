# Constructing Groups

Three constructors cover the common cases. Whatever you build is an
`AbstractFiniteGroup` and supports the same accessors — [`name`](@ref),
[`order`](@ref), [`class`](@ref), [`inclass`](@ref), [`mult`](@ref) — and the
same character-table and representation machinery.

## Crystallographic point groups

[`pointgroup`](@ref) returns one of the 32 crystallographic point groups, by
Schoenflies symbol or by index `1:32`:

```@example groups
using FiniteGroups

g = pointgroup("Oh")        # full octahedral group
(; group_name=name(g), group_order=order(g), num_classes=length(class(g)))
```

The valid symbols, in index order, are

```
C1, Ci, C2, Cs, C2h, D2, C2v, D2h, C4, S4, C4h, D4, C4v, D2d, D4h,
C3, C3i, D3, C3v, D3d, C6, C3h, C6h, D6, C6v, D3h, D6h, T, Th, O, Td, Oh
```

Point groups carry bundled reference data: the standard character table (used
by [`charactertable`](@ref) unless you request another method), conventional
irrep labels via [`repname`](@ref), and the symmetry operations' rotation
matrices via [`matrix`](@ref) (crystal-axis basis) and [`rotation`](@ref)
(Cartesian). [`groupindex`](@ref) and [`operation`](@ref) recover a group's
index and the list of operation symbols.

## Symmetric and permutation groups

[`permutationgroup`](@ref) builds the symmetric group ``S_n`` from an integer,
or closes a set of generating permutations into a group. Permutations are
written with [`cycles`](@ref) (disjoint-cycle notation) or [`permutation`](@ref)
(an explicit ``a \mapsto b`` table):

```@example groups
s4 = permutationgroup(4)        # the symmetric group S₄
```

```@example groups
# the same group, generated from a transposition and a 4-cycle
g = permutationgroup([cycles((1, 2)), cycles((1, 2, 3, 4))], name="S4")
order(g)
```

## From a multiplication table

For any finite group, hand [`FiniteGroup`](@ref) its multiplication table — a
matrix whose entry `[i, j]` is the index of ``g_i g_j`` — together with a name.
Here is the cyclic group ``C_4``, whose characters are genuinely complex:

```@example groups
multab = [1 2 3 4; 2 3 4 1; 3 4 1 2; 4 1 2 3]
c4 = FiniteGroup(multab, "C4")
charactertable(c4)
```

[`check_group`](@ref) validates that such a multiplication table really defines
a group — a two-sided identity, associativity, and a permutation in every row
and column. It takes the table itself:

```@example groups
check_group(multab)
```
