# Getting Started With `LandauGinzburgCategories.jl`

## Installation

Refer to [the Julia website](https://julialang.org/downloads/) for details on
installing Julia. As soon as you have, start it and type `[` to get in
package mode. Then, run
```julia-repl
(v1.1) pkg> add https://github.com/tkluck/LandauGinzburgCategories.jl
```
to install `LandauGinzburgCategories` and its dependencies. To test whether it worked,
type

```@repl getting-started
using LandauGinzburgCategories.Library
orbifold_equivalence(TwoVariables.A{9}, TwoVariables.D{6})
```
If you see a matrix factorization, the installation worked!

## Overview



## Frequently Asked Questions

Be the first to ask a question! Feel free to [open an issue for it](https://github.com/tkluck/LandauGinzburgCategories.jl/issues/new).
