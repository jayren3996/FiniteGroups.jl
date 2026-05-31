<p align="center">
  <img src="docs/src/assets/logo.svg" alt="FiniteGroups.jl" width="180"/>
</p>

<h1 align="center">FiniteGroups.jl</h1>

<p align="center">
  <em>Character tables and representations of finite groups, in Julia.</em>
</p>

<p align="center">
  <a href="https://jayren3996.github.io/FiniteGroups.jl/dev/"><img alt="Documentation (dev)" src="https://img.shields.io/badge/docs-dev-blue.svg"/></a>
  <a href="https://github.com/jayren3996/FiniteGroups.jl/actions/workflows/CI.yml"><img alt="CI status" src="https://github.com/jayren3996/FiniteGroups.jl/actions/workflows/CI.yml/badge.svg"/></a>
  <a href="https://github.com/jayren3996/FiniteGroups.jl/actions/workflows/documentation.yml"><img alt="Documentation build status" src="https://github.com/jayren3996/FiniteGroups.jl/actions/workflows/documentation.yml/badge.svg"/></a>
  <a href="https://julialang.org/"><img alt="Julia" src="https://img.shields.io/badge/made%20with-Julia-9558B2.svg?logo=julia"/></a>
  <a href="LICENSE"><img alt="License: MIT" src="https://img.shields.io/badge/license-MIT-yellow.svg"/></a>
</p>

---

`FiniteGroups.jl` computes the representation theory of finite groups: conjugacy
classes, character tables, and explicit matrix representations — including the
real and projective (double-group) representations that appear in
condensed-matter and molecular-symmetry problems. Hand it a crystallographic
point group, a symmetric group, or any group given as a multiplication table.

The 32 crystallographic point groups ship with reference data — standard
character tables, conventional Mulliken labels, and rotation matrices — so for
those you get publication-style tables immediately. For an arbitrary group
everything is computed from the multiplication table: a fast floating-point
method (Burnside) by default, or an exact finite-field method (Dixon) when you
want exact character values rather than floating point.

## Installation

`FiniteGroups` is registered in Julia's General registry:

```julia
pkg> add FiniteGroups
```

Or directly from GitHub:

```julia
pkg> add https://github.com/jayren3996/FiniteGroups.jl
```

Then load it with `using FiniteGroups`.

## Quick start

```julia
julia> using FiniteGroups

julia> g = pointgroup("D3")          # a crystallographic point group, by Schoenflies symbol
Point group : D3
Group order : 6
Classes     : 3

julia> charactertable(g)             # standard table, with conventional Mulliken labels
4×4 Matrix{Any}:
 ""     "1"    "3₀₀₁⁺"    "2₋₁₁₀"
 "A1"  1      1          1
 "A2"  1      1         -1
 "E"   2     -1          0

julia> reps = irreps(g);             # explicit representation matrices, one vector per irrep

julia> length(reps), [size(r[1], 1) for r in reps]   # 3 irreps, of dimension 1, 1, 2
(3, [1, 1, 2])
```

## Features

- **Built-in groups** — all 32 crystallographic point groups (`pointgroup`) with
  standard character tables, Mulliken labels (`repname`), and rotation matrices
  (`matrix`, `rotation`); symmetric groups `Sₙ` and generated permutation groups
  (`permutationgroup`, `cycles`, `permutation`).
- **Any finite group** — build one from a multiplication table with
  `FiniteGroup`, and validate an untrusted table with `check_group`.
- **Character tables, two ways** — `charactertable` via a fast floating-point
  class-algebra method (Burnside), or `method=:dixon` for an exact result over a
  finite field; `dixon` returns the raw exact character matrix.
- **Representations** — `irreps` for complex unitary irreps, `real_irreps` /
  `real_rep` for real forms (with the Frobenius–Schur indicator from
  `check_real_rep`), and `regular_rep` / `proj_operator` / `block_decomposition`
  for building and decomposing representations.
- **Projective & double-group reps** — `proj_reps`, `chiral_proj_reps`,
  `cover_group`, and `check_proj_coeff` for spin-½ / double-group problems.

## Two character-table methods

`charactertable(g)` uses **Burnside's method** by default: it simultaneously
diagonalizes the conjugacy-class algebra and reads off the irreducible
characters. It is fast, but the entries are floating point and carry a small
numerical tolerance.

Passing `method=:dixon` switches to **Dixon's method**, which runs the same
computation in a finite field 𝔽ₚ — chosen so the representations are defined
there — and lifts the result back to characteristic zero. The characters are
then obtained exactly, with no floating-point eigenvalue tolerance:

```julia
julia> charactertable(pointgroup("D3"); method=:dixon)   # same table, computed exactly
4×4 Matrix{Any}:
 ""    "1"    "3₀₀₁⁺"    "2₋₁₁₀"
 "1"  1.0    1.0        1.0
 "2"  1.0    1.0       -1.0
 "3"  2.0   -1.0        0.0
```

For the 32 point groups the default is the bundled standard table
(`method=:table`); `:burnside` and `:dixon` recompute it from scratch.

## Documentation

Full manual at <https://jayren3996.github.io/FiniteGroups.jl/dev/>.

- [Getting Started](https://jayren3996.github.io/FiniteGroups.jl/dev/getting-started/) — installation, the first character table, accessors, irreps.
- [Constructing Groups](https://jayren3996.github.io/FiniteGroups.jl/dev/groups/) — point groups, symmetric groups, and raw multiplication tables.
- [Character Tables](https://jayren3996.github.io/FiniteGroups.jl/dev/character-tables/) — the `:table`, `:burnside`, and `:dixon` methods.
- [Representations](https://jayren3996.github.io/FiniteGroups.jl/dev/representations/) — real / complex / projective representations and decomposition.
- [API Reference](https://jayren3996.github.io/FiniteGroups.jl/dev/api/) — generated from the docstrings.

## License

MIT — see [LICENSE](LICENSE).
