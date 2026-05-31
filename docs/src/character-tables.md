# Character Tables

[`charactertable`](@ref) is the main entry point. It returns a `CharacterTable`
whose rows are the irreducible characters and whose columns are the conjugacy
classes. Indexing the table gives entries, rows, or columns, and `Matrix(ct)`
extracts the plain numeric array.

```@example ct
using FiniteGroups

g = pointgroup("C4v")
ct = charactertable(g)
```

```@example ct
ct[2, :]            # the second irrep's character, across the classes
```

## Choosing a method

Three algorithms are available through the `method` keyword:

- **`:table`** — the standard table bundled with the package, with conventional
  irrep labels. This is the default for point groups, and is only available for
  them.
- **`:dixon`** — Dixon's method: the table is computed *exactly* by working in a
  finite field ``\mathbb{F}_p`` and lifting the result back to characteristic
  zero. This is the default for every group other than the point groups.
- **`:burnside`** — a floating-point method that simultaneously diagonalizes the
  class algebra. Entries are `Float64` / `ComplexF64` and carry small numerical
  error; pass `method=:burnside` to use it instead of the exact default.

The bundled and computed tables agree (up to the ordering of the irreps). Here
is the same group via Dixon's exact method:

```@example ct
charactertable(g; method=:dixon)
```

For a group given only by its multiplication table, the bundled `:table` does
not apply, but `:dixon` (default) and `:burnside` both do:

```@example ct
s5 = permutationgroup(5)            # 7 classes ⇒ 7 irreducible characters
charactertable(s5; method=:dixon)
```

The standalone [`dixon`](@ref) function returns the raw exact character matrix
(as a `Matrix`) without wrapping it in a `CharacterTable`, which is convenient
when you only want the numbers.

## The class multiplication table

Both solvers are driven by the class-algebra structure constants. These are
exposed through [`class_multab`](@ref), which returns the structure constants of
the conjugacy-class algebra; you can use them to build custom solvers or to
inspect the algebra directly.
