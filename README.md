# FiniteGroups.jl

A Julia package for finite group theory: conjugacy classes, character tables,
(numerical) irreducible representations, realification, and projective /
double-group representations, with the 32 crystallographic point groups and the
symmetric groups built in.

## Installation

In the Julia `REPL`, run:

```julia
using Pkg
Pkg.add("FiniteGroups")
```

Or install directly from the GitHub URL:

```julia
using Pkg
Pkg.add(url="https://github.com/jayren3996/FiniteGroups.jl")
```

## Examples

All examples assume the package is loaded:

```julia
using FiniteGroups
```

### Create a finite group

Create a point group by its number or its Schoenflies name. For example:

```julia
julia> g = pointgroup(32)
Point group : Oh
Group order : 48
Classes     : 10

julia> g = pointgroup("Th")
Point group : Th
Group order : 24
Classes     : 8
```

The symmetric group `Sₙ` is available via `permutationgroup`:

```julia
julia> permutationgroup(4)
Permutation group : S4
Group order       : 24
Classes           : 5
```

In general, given a multiplication table `multab`, create the group with the
`FiniteGroup` constructor (call `check_group` first to validate an untrusted
table):

```julia
# Multiplication table of point group D3:
multab = [
    1  2  3  4  5  6
    2  3  1  6  4  5
    3  1  2  5  6  4
    4  5  6  1  2  3
    5  6  4  3  1  2
    6  4  5  2  3  1
]
g = FiniteGroup(multab, "D3")
```

### Character table

Compute the character table of a group with `charactertable`:

```julia
ctable = charactertable(g)
```

For the point group `Oh` the result is:

```julia
julia> charactertable(pointgroup("Oh"))
11×11 Matrix{Any}:
 ""      "1"    "2₀₀₁"    "3₁₁₁⁺"    "2₁₁₀"    "4₀₀₁⁻"    "-1"    "m₀₀₁"    "-3₁₁₁⁺"    "m₁₁₀"    "-4₀₀₁⁻"
 "A1g"  1      1         1          1         1          1       1         1           1         1
 "A1u"  1      1         1          1         1         -1      -1        -1          -1        -1
 "A2g"  1      1         1         -1        -1          1       1         1          -1        -1
 "A2u"  1      1         1         -1        -1         -1      -1        -1           1         1
 "Eg"   2      2        -1          0         0          2       2        -1           0         0
 "Eu"   2      2        -1          0         0         -2      -2         1           0         0
 "T2g"  3     -1         0          1        -1          3      -1         0           1        -1
 "T2u"  3     -1         0          1        -1         -3       1         0          -1         1
 "T1g"  3     -1         0         -1         1          3      -1         0          -1         1
 "T1u"  3     -1         0         -1         1         -3       1         0           1        -1
```

The returned `ctable` is a `CharacterTable`. Indexing with a single integer
returns the `Characters` of one irrep (which also records whether it is real,
complex, or pseudo-real):

```julia
julia> ctable = charactertable(pointgroup("Oh"));

julia> ctable[10]
Characters of Real representation of Oh:
[3, -1, 0, -1, 1, -3, 1, 0, 1, -1]
```

It can also be sliced like a matrix; rows and submatrices come back as plain
arrays:

```julia
julia> ctable[3, :]
10-element Vector{Int64}:
  1
  1
  1
 -1
 -1
  1
  1
  1
 -1
 -1

julia> ctable[1:3, :]
3×10 Matrix{Int64}:
 1  1  1   1   1   1   1   1   1   1
 1  1  1   1   1  -1  -1  -1  -1  -1
 1  1  1  -1  -1   1   1   1  -1  -1
```

Groups with complex irreducible representations return complex characters, e.g.
the point group `T`:

```julia
julia> charactertable(pointgroup("T"))
5×5 Matrix{Any}:
 ""        "1"         "2₀₀₁"       "3₁₁₁⁺"          "3₁₁₁⁻"
 "A"   1.0+0.0im   1.0+0.0im    1.0+0.0im        1.0+0.0im
 "1E"  1.0+0.0im   1.0+0.0im   -0.5-0.866025im  -0.5+0.866025im
 "2E"  1.0+0.0im   1.0+0.0im   -0.5+0.866025im  -0.5-0.866025im
 "T"   3.0+0.0im  -1.0+0.0im    0.0+0.0im        0.0+0.0im
```

### Irreducible representations

Compute all irreducible representations of a group with `irreps`. It returns one
entry per irrep; each entry is a vector of representation matrices indexed by
group element:

```julia
julia> g = pointgroup("T");

julia> reps = irreps(g);

julia> length(reps)          # T has 4 irreducible representations
4

julia> T = reps[end];        # the 3-dimensional irrep

julia> length(T), size(T[1]) # 12 group elements, each a 3×3 matrix
(12, (3, 3))

julia> T[1]                  # the identity element
3×3 Matrix{Float64}:
 1.0  0.0  0.0
 0.0  1.0  0.0
 0.0  0.0  1.0

julia> T[5]
3×3 Matrix{Float64}:
 0.0  0.0  1.0
 1.0  0.0  0.0
 0.0  1.0  0.0
```

Pass `R=true` to obtain real (orthogonal) irreps, projecting complex ones into
their real form:

```julia
reps = irreps(g; R=true)
```

A reducible representation can be split back into irreducible blocks with
`block_decomposition`.
