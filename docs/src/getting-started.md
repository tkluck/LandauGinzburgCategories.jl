# Getting Started With `LandauGinzburgCategories.jl`

## Installation

Refer to [the Julia website](https://julialang.org/downloads/) for details on
installing Julia. As soon as you have, start it and type `]` to get in
package mode. Then, run
```julia-repl
(v1.1) pkg> add https://github.com/tkluck/LandauGinzburgCategories.jl#release
```
to install `LandauGinzburgCategories` and its dependencies. To test whether it worked,
type

```@repl getting-started
using LandauGinzburgCategories.Library
orbifold_equivalence(TwoVariables.A{9}, TwoVariables.D{6})
```
If you see a matrix factorization, the installation worked!

## Overview

### Composing two matrix factorizations

Let us compose two simple matrix factorizations. For the first, let's take
a unit matrix factorization:

```@repl getting-started
using LandauGinzburgCategories, LandauGinzburgCategories.Library;
@ring! Int[x,y,s,t,u,v];
A = unit_matrix_factorization(x^6 - y^2, x => s, y => t)
```

As you see, we first load the package by writing `using ....`. Then, we declare
the variables `x,y,s,t,u,v` that we will use. The last step calls the function
[`unit_matrix_factorization`](@ref). The result is a matrix factorization:

```@repl getting-started
A^2
```

You can compute the left- and right quantum dimensions using [`quantum_dimensions`](@ref):
```@repl getting-started
quantum_dimensions(A, (x,y), (s,t))
```
As you see, you need to pass the sets of left (`x,y`) and right (`s,t`) variables.

The polynomial ``x^6 - y^2`` is a representative of the ``A_5`` equivalence class.
We now load an orbifold equivalence from the library that we can compose with it.

```@repl getting-started
B = orbifold_equivalence(TwoVariables.A{5}, TwoVariables.D{4}, (s, t), (u, v))
```

The function [`fuse`](@ref) can compute the composition for us. The following
command fuses along the `s` and `t` variables.

```@repl getting-started
fuse(A, B, s, t)
```

The result is a sparse matrix; for easier viewing we can call `collect` on it.

```@repl getting-started
collect(ans)
```

### Computing the central charge of a potential

The function [`centralcharge`](@ref) computes the central charge. You can
use it as follows:

```@repl getting-started
centralcharge(x^6 - x*y^2, x, y)
```

As you can see, it is necessary to also pass the variables with respect to
which you want to compute the central charge -- in this case `x,y`. For
example, this is useful when you want to compute the central charge of a
parametrized family:

```@repl getting-started
centralcharge(x^6 + t*x*y^2, x, y)
```

## Frequently Asked Questions

Be the first to ask a question! Feel free to [open an issue for it](https://github.com/tkluck/LandauGinzburgCategories.jl/issues/new).
