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
Finite group : Oh
Group order  : 48
Classes      : 10
```

or use the point group name (for example Th group):

```julia
julia> g = pointgroup("Th")
Finite group : Th
Group order  : 24
Classes      : 8
```

In general, given the multiplication table `multab` of the group, we can create the group object using the command:

```julia
g = finitegroup(multab)
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

For example, consider the group D4:

```julia
g = pointgroup("D4")
representations = irreps(g)
```

The result is a 5-element vector, of type `Representation`. The group D4 has 4 1-d irreducible representations and one 2-d irreducible representation. We can display the 2-d one by:

```julia
julia> representations[5]
2-d irreducible representation of D4: 

Element 1:
2×2 Matrix{Int64}:
 1  0
 0  1

Element 2:
2×2 Matrix{Number}:
 0-1im   0
  0     0+1im

Element 3:
2×2 Matrix{Int64}:
 -1   0
  0  -1

Element 4:
2×2 Matrix{Number}:
 0+1im   0
  0     0-1im

Element 5:
2×2 Matrix{Int64}:
 0  1
 1  0

Element 6:
2×2 Matrix{Number}:
  0     0-1im
 0+1im   0

Element 7:
2×2 Matrix{Int64}:
  0  -1
 -1   0

Element 8:
2×2 Matrix{Number}:
  0     0+1im
 0-1im   0
```

