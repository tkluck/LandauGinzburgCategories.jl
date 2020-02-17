# `LandauGinzburgCategories.jl` Documentation

## Features

This package aims to help you experiment with matrix factorizations
in the bicategory of Landau-Ginzburg models. It contains implementations
of operations such as _unit_, _composition_, _fusion of variables_,
_duality_, and _left/right quantum dimension_. Future additions may
include _(co)evaluation maps_.

In addition, a library of well-known named potentials, as well as
known orbifold equivalences, is provided.

For getting started, have a look
at [Getting Started With `LandauGinzburgCategories.jl`](@ref).

## Limitations

**Expensive computations.** Many known orbifold equivalences are presented
by way of an _Ansatz_ together with the verifiable statement that the
resulting equations have a solution. We represent this by placing these
coefficients in a quotient ring. Every operation in this quotient ring
necessarily involves a reduction step, and this can be noticably slow.


## Table of contents
```@contents
Pages = [
    # keep in sync with make.jl
    "index.md",
    "getting-started.md",
    "functions.md",
    "reference.md",
]
Depth = 3
```
