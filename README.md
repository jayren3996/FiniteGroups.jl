# FiniteGroups.jl
 Julia package for finite group theory calculation.

## Installation
In julia `REPL`, run the following script:

```julia
using Pkg
Pkg.add("FiniteGroups")
```

Or, install the package directly from the GitHub URL:

```julia
using Pkg
Pkg.add(url="https://github.com/jayren3996/FiniteGroups.jl")
```

## Examples

### Create a Finite Group

We can creat a point group using the group number or group name. For example, the following command:

```julia
julia> g = pointgroup(32)
Point group : Oh
Group order : 48
Classes     : 10
```

or use the point group name (for example Th group): 

```julia
julia> g = pointgroup("Th")
Point group : Th
Group order : 24
Classes     : 8
```

In general, given the multiplication table `multab` of the group, we can create the group object using the command:

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
g = FiniteGroup(multab)
```

### Character Table

We can calculate the character table of a finite group using the command

```julia
tab = character(g)
```

If group `g` is chosen to be the pointgroup Oh, the displayed result is:

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

The `chartable` is of type `CharacterTable`, from which we can extract a specific set of characters:

```julia
julia> ctable[10]
Characters of Real representation of Oh:
[3.0, -1.0, 0.0, -1.0, 1.0, -3.0, 1.0, 0.0, 1.0, -1.0]
```

The `CharacterTable` can be sliced as a matrix:

```julia
julia> ctable[3,:]
10-element Vector{Float64}:
  1.0
  1.0
  1.0
 -1.0
 -1.0
  1.0
  1.0
  1.0
 -1.0
 -1.0
julia> ctable[1:3,:]
3×10 Matrix{Float64}:
 1.0  1.0  1.0   1.0   1.0   1.0   1.0   1.0   1.0   1.0
 1.0  1.0  1.0   1.0   1.0  -1.0  -1.0  -1.0  -1.0  -1.0
 1.0  1.0  1.0  -1.0  -1.0   1.0   1.0   1.0  -1.0  -1.0
```

### Irreducible Representation

We can also compute all irreducible representations of a finite group `g`, using the command `irreps`. For example, the character table for point group T is:

```julia
julia> g = pointgroup("T"); charactertable(g)
5×5 Matrix{Any}:
 ""        "1"         "2₀₀₁"       "3₁₁₁⁺"          "3₁₁₁⁻"
 "A"   1.0+0.0im   1.0+0.0im    1.0+0.0im        1.0+0.0im
 "1E"  1.0+0.0im   1.0+0.0im   -0.5-0.866025im  -0.5+0.866025im
 "2E"  1.0+0.0im   1.0+0.0im   -0.5+0.866025im  -0.5-0.866025im
 "T"   3.0+0.0im  -1.0+0.0im    0.0+0.0im        0.0+0.0im
```

We see there is a three-dimensional T representation. To obtain the representation matrices, simply using the following command:

```julia
rep = irreps(g)[end]
```

The output is a list of matrices:

```julia
julia> display.(rep)
3×3 Matrix{Float64}:
 1.0  0.0  0.0
 0.0  1.0  0.0
 0.0  0.0  1.0
3×3 Matrix{Float64}:
 1.0   0.0   0.0
 0.0  -1.0   0.0
 0.0   0.0  -1.0
3×3 Matrix{Float64}:
 -1.0   0.0  0.0
  0.0  -1.0  0.0
  0.0   0.0  1.0
3×3 Matrix{Float64}:
 -1.0  0.0   0.0
  0.0  1.0   0.0
  0.0  0.0  -1.0
3×3 Matrix{Float64}:
 0.0  0.0  1.0
 1.0  0.0  0.0
 0.0  1.0  0.0
3×3 Matrix{Float64}:
 0.0   0.0  -1.0
 1.0   0.0   0.0
 0.0  -1.0   0.0
3×3 Matrix{Float64}:
  0.0   0.0  1.0
 -1.0   0.0  0.0
  0.0  -1.0  0.0
3×3 Matrix{Float64}:
  0.0  0.0  -1.0
 -1.0  0.0   0.0
  0.0  1.0   0.0
3×3 Matrix{Float64}:
 0.0  1.0  0.0
 0.0  0.0  1.0
 1.0  0.0  0.0
3×3 Matrix{Float64}:
 0.0  -1.0   0.0
 0.0   0.0  -1.0
 1.0   0.0   0.0
3×3 Matrix{Float64}:
  0.0  -1.0  0.0
  0.0   0.0  1.0
 -1.0   0.0  0.0
3×3 Matrix{Float64}:
  0.0  1.0   0.0
  0.0  0.0  -1.0
 -1.0  0.0   0.0
```

