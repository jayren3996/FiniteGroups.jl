# Getting Started

## Installation

`FiniteGroups` is registered in the General registry, so it installs the usual
way:

```julia
using Pkg
Pkg.add("FiniteGroups")
```

or, from the Pkg REPL (press `]`):

```
pkg> add FiniteGroups
```

## A first character table

```@example gs
using FiniteGroups

g = pointgroup("C3v")
charactertable(g)
```

The rows are the irreducible representations, labelled with their conventional
Mulliken names, and the columns are the conjugacy classes. To get the plain
numeric matrix, wrap the table in `Matrix`:

```@example gs
Matrix(charactertable(g))
```

## Group properties

However a group was built, it answers the same small set of queries:

```@example gs
(; group_name=name(g), group_order=order(g),
   num_classes=length(class(g)), class_sizes=mult(g))
```

- [`name`](@ref), [`order`](@ref) — the group's name and order ``|G|``.
- [`class`](@ref) — the conjugacy classes, each a vector of element indices.
- [`inclass`](@ref) — the class index of a given element.
- [`mult`](@ref) — the size of each conjugacy class.

## Irreducible representations

`charactertable` gives the characters; [`irreps`](@ref) gives the
representation matrices themselves. Each entry of the returned vector is one
irrep — itself a vector of matrices, one matrix per group element.

```@example gs
rs = irreps(g)
[size(r[1], 1) for r in rs]      # the dimension of each irrep
```

So `C₃ᵥ` has two one-dimensional irreps and one two-dimensional irrep, and
`rs[3][k]` is the ``2 \times 2`` matrix representing the ``k``-th group element
in the `E` representation.

## Where to next

- [Constructing Groups](groups.md) — point groups, symmetric groups, and raw
  multiplication tables.
- [Character Tables](character-tables.md) — the `:table`, `:burnside`, and
  `:dixon` methods.
- [Representations](representations.md) — real/complex/projective
  representations and decomposition.
