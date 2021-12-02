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

If group `g` is chosen to be the pointgroup Oh, the result is:

```julia
julia> charactertable(g)
10×10 Matrix{Int64}:
 1   1   1   1   1   1   1   1   1   1
 1   1   1  -1  -1  -1  -1   1   1  -1
 1  -1   1   1  -1  -1   1   1  -1   1
 1  -1   1  -1   1   1  -1   1  -1  -1
 2   0  -1   2   0   0  -1   2   0   2
 2   0  -1  -2   0   0   1   2   0  -2
 3   1   0   1  -1   1   0  -1  -1  -3
 3   1   0  -1   1  -1   0  -1  -1   3
 3  -1   0  -1  -1   1   0  -1   1   3
 3  -1   0   1   1  -1   0  -1   1  -3
```

### Irreducible Representation

We can also compute all irreducible representations of a finite group `g`, using the command:

```julia
representations = irreps(g)
```

