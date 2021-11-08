# FiniteGroups.jl
 Julia package for finite group theory calculation.

## Installation
In julia `REPL`, run the following script:

```julia
using Pkg
Pkg.add("FiniteGroups")
```

## Examples

### Create a Finite Group

We can creat a point group using the group number or group name. For example, the following command:

```julia
g = pointgroup(32)
```

or

```julia
g = pointgroup("Oh")
```

will all return an Oh group object.

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
10×10 Matrix{Float64}:
 1.0   1.0   1.0   1.0   1.0   1.0   1.0   1.0   1.0   1.0
 1.0   1.0   1.0  -1.0  -1.0  -1.0  -1.0   1.0   1.0  -1.0
 1.0  -1.0   1.0   1.0  -1.0  -1.0   1.0   1.0  -1.0   1.0
 1.0  -1.0   1.0  -1.0   1.0   1.0  -1.0   1.0  -1.0  -1.0
 2.0   0.0  -1.0   2.0   0.0   0.0  -1.0   2.0   0.0   2.0
 2.0  -0.0  -1.0  -2.0  -0.0  -0.0   1.0   2.0   0.0  -2.0
 3.0   1.0   0.0   1.0  -1.0   1.0  -0.0  -1.0  -1.0  -3.0
 3.0   1.0   0.0  -1.0   1.0  -1.0  -0.0  -1.0  -1.0   3.0
 3.0  -1.0   0.0  -1.0  -1.0   1.0   0.0  -1.0   1.0   3.0
 3.0  -1.0  -0.0   1.0   1.0  -1.0  -0.0  -1.0   1.0  -3.0
```

### Irreducible Representation

We can also compute all irreducible representations of a finite group `g`, using the command:

```julia
representations = irreps(g)
```

For example, consider the group D4:

```
g = pointgroup("D4")
representations = irreps(g)
```

The result is a 5-element vector, of type `Representation`. The group D4 has 4 1-d irreducible representations and one 2-d irreducible representation. We can display the 2-d one by:

```julia
julia> representations[5]
2-d irreducible representation of D4: 

Element 1:
2×2 Matrix{ComplexF64}:
         1.0+0.0im          1.57009e-16-1.12686e-32im
 1.57009e-16+1.12686e-32im          1.0+0.0im

Element 2:
2×2 Matrix{ComplexF64}:
  6.84228e-49-1.0im          1.12686e-32+1.57009e-16im
 -1.12686e-32+1.57009e-16im          0.0+1.0im

Element 3:
2×2 Matrix{ComplexF64}:
         -1.0+0.0im          -1.57009e-16+1.12686e-32im
 -1.57009e-16-1.12686e-32im          -1.0+3.08149e-33im

Element 4:
2×2 Matrix{ComplexF64}:
         0.0+1.0im          -1.12686e-32-1.57009e-16im
 1.12686e-32-1.57009e-16im  -3.08149e-33-1.0im

Element 5:
2×2 Matrix{ComplexF64}:
 -3.14018e-16+0.0im          -1.0+1.2326e-32im
         -1.0-9.24446e-33im   0.0+0.0im

Element 6:
2×2 Matrix{ComplexF64}:
 -2.25371e-32-1.71057e-49im  1.2326e-32+1.0im
   1.2326e-32-1.0im                 0.0+0.0im

Element 7:
2×2 Matrix{ComplexF64}:
 3.14018e-16+0.0im         1.0-1.2326e-32im
         1.0+1.2326e-32im  0.0+0.0im

Element 8:
2×2 Matrix{ComplexF64}:
  2.25371e-32+1.71057e-49im  -1.2326e-32-1.0im
 -9.24446e-33+1.0im                  0.0+0.0im
```

