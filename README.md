<div align="center">

<img src="docs/src/assets/logo.svg" alt="FiniteGroups.jl logo" width="170"/>

# FiniteGroups.jl

**Character tables and representations of finite groups, in Julia.**

From the 32 crystallographic point groups to any group given as a multiplication table.

[![Docs](https://img.shields.io/badge/docs-dev-9558B2.svg)](https://jayren3996.github.io/FiniteGroups.jl/dev/) [![CI](https://github.com/jayren3996/FiniteGroups.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/jayren3996/FiniteGroups.jl/actions/workflows/CI.yml) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE) [![Julia](https://img.shields.io/badge/Julia-1.10%2B-389826.svg)](https://julialang.org)

</div>

---

FiniteGroups.jl computes the representation theory of finite groups: conjugacy classes,
character tables, and explicit matrix representations — including the real and projective
(double-group) representations that appear in condensed-matter and molecular-symmetry problems.
It is written for researchers who want a precise Julia interface with exact, documented
conventions. Hand it a crystallographic point group, a symmetric group, or any group given as a
multiplication table.

## ✨ Features

|  |  |
| --- | --- |
| 🧱 **Built-in groups** | All 32 crystallographic point groups (`pointgroup`) with standard tables, Mulliken labels (`repname`), and rotation matrices (`matrix`, `rotation`); symmetric groups `Sₙ` and generated permutation groups (`permutationgroup`, `cycles`, `permutation`). |
| 🔢 **Any finite group** | Build one from a multiplication table with `FiniteGroup`, and validate an untrusted table with `check_group`. |
| 📊 **Character tables, two ways** | `charactertable` via an exact finite-field method (Dixon) by default, or `method=:burnside` for a floating-point alternative; `dixon` returns the raw exact character matrix. |
| 🎛️ **Representations** | `irreps` for complex unitary irreps, `real_irreps` / `real_rep` for real forms (with the Frobenius–Schur indicator from `check_real_rep`), and `regular_rep` / `proj_operator` / `block_decomposition` for building and decomposing them. |
| 🌀 **Projective & double-group reps** | `proj_reps`, `chiral_proj_reps`, `cover_group`, and `check_proj_coeff` for spin-½ / double-group problems. |
| 🎯 **Exact by default** | Dixon's finite-field method gives character values with no eigenvalue tolerance; the 32 point groups ship with publication-style reference data. |

## 📦 Installation

FiniteGroups is registered in Julia's General registry:

```julia
pkg> add FiniteGroups
```

Or directly from GitHub:

```julia
pkg> add https://github.com/jayren3996/FiniteGroups.jl
```

## 🚀 Quick Start

Build a crystallographic point group, print its character table, and pull out the irrep matrices:

```julia
using FiniteGroups

g = pointgroup("D3")          # a point group, by Schoenflies symbol
charactertable(g)             # standard table, with conventional Mulliken labels

reps = irreps(g)              # explicit representation matrices, one vector per irrep
length(reps), [size(r[1], 1) for r in reps]   # 3 irreps, of dimension 1, 1, 2
```

```julia
julia> charactertable(pointgroup("D3"))
4×4 Matrix{Any}:
 ""     "1"    "3₀₀₁⁺"    "2₋₁₁₀"
 "A1"  1      1          1
 "A2"  1      1         -1
 "E"   2     -1          0
```

## 🧭 Choosing What to Call

| Modeling need | Use |
| --- | --- |
| A crystallographic point group | `pointgroup("D3")` |
| A symmetric or permutation group | `permutationgroup`, `cycles`, `permutation` |
| Any group from a multiplication table | `FiniteGroup` (validate with `check_group`) |
| A character table (exact) | `charactertable(g)` — Dixon by default |
| A character table (fast, floating-point) | `charactertable(g; method=:burnside)` |
| The raw exact character matrix | `dixon(g)` |
| Complex unitary irreps | `irreps(g)` |
| Real irreps / Frobenius–Schur indicator | `real_irreps`, `real_rep`, `check_real_rep` |
| Projective / double-group reps | `proj_reps`, `chiral_proj_reps`, `cover_group` |

For the 32 point groups, `charactertable` returns the bundled standard table (`method=:table`)
by default; `:burnside` and `:dixon` recompute it from scratch.

## 🛠 Usage

<details>
<summary><b>Built-in point groups</b></summary>

<br>

The 32 crystallographic point groups ship with reference data — standard character tables,
conventional Mulliken labels, and rotation matrices — so you get publication-style output
immediately:

```julia
using FiniteGroups

g = pointgroup("D3")
g                       # summary: name, order, number of classes
charactertable(g)       # standard table with Mulliken labels (repname)
matrix(g)               # rotation matrices of the group elements
```

</details>

<details>
<summary><b>Any finite group from a multiplication table</b></summary>

<br>

For an arbitrary group, everything is computed from the multiplication table. The character
table comes from an exact finite-field class-algebra method (Dixon) by default:

```julia
using FiniteGroups

g = FiniteGroup(mult_table)    # mult_table[i, j] = index of gᵢ * gⱼ
check_group(g)                 # validate an untrusted table

charactertable(g)              # exact, via Dixon's method
charactertable(g; method=:burnside)   # floating-point alternative, faster on large groups
dixon(g)                       # raw exact character matrix
```

Dixon's method diagonalizes the conjugacy-class algebra in a finite field 𝔽ₚ — chosen so the
representations are defined there — and recovers each character as an exact integer combination
of roots of unity, with no floating-point eigenvalue tolerance. Burnside's method performs the
same diagonalization with floating-point eigenvalues: faster on large groups, at the cost of a
small numerical tolerance.

</details>

<details>
<summary><b>Representations: complex, real, and decomposition</b></summary>

<br>

`irreps` returns the complex unitary irreps; real forms and the Frobenius–Schur indicator are
one call away, and the projection / block-decomposition tools build and reduce representations:

```julia
using FiniteGroups

g = pointgroup("D3")

reps = irreps(g)               # complex unitary irreps, one matrix per element
real_irreps(g)                 # real forms of the irreps
check_real_rep(g)              # Frobenius–Schur indicators

regular_rep(g)                 # the regular representation
block_decomposition(g, ρ)      # reduce a representation ρ into irreducible blocks
```

</details>

<details>
<summary><b>Projective & double-group representations</b></summary>

<br>

Spin-½ and double-group problems need projective representations. FiniteGroups builds them from
the covering group and lets you check the cocycle:

```julia
using FiniteGroups

g = pointgroup("D3")

proj_reps(g)                   # projective (double-group) representations
chiral_proj_reps(g)            # chiral projective representations
cover_group(g)                 # the covering group
check_proj_coeff(g, ω)         # verify a projective phase factor (cocycle)
```

</details>

## 📚 Documentation

Full documentation lives at **[jayren3996.github.io/FiniteGroups.jl](https://jayren3996.github.io/FiniteGroups.jl/dev/)**.

- [Getting Started](https://jayren3996.github.io/FiniteGroups.jl/dev/getting-started/) — installation, the first character table, accessors, irreps.
- [Constructing Groups](https://jayren3996.github.io/FiniteGroups.jl/dev/groups/) — point groups, symmetric groups, and raw multiplication tables.
- [Character Tables](https://jayren3996.github.io/FiniteGroups.jl/dev/character-tables/) — the `:table`, `:burnside`, and `:dixon` methods.
- [Representations](https://jayren3996.github.io/FiniteGroups.jl/dev/representations/) — real / complex / projective representations and decomposition.
- [API Reference](https://jayren3996.github.io/FiniteGroups.jl/dev/api/) — generated from the docstrings.

## 📄 License

[MIT](LICENSE) © 2021 JieRen and contributors.
