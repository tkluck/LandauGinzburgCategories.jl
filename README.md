# Landau-Ginzburg categories for Julia

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://tkluck.github.io/LandauGinzburgCategories.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://tkluck.github.io/LandauGinzburgCategories.jl/latest)


| **Build Status**                                                | **Test coverage**                                       |
|:---------------------------------------------------------------:|:-------------------------------------------------------:|
| [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] | [![Coverage Status][codecov-img]][codecov-url]      |

A library for algebra inside the bicategory of Landau-Ginzburg models.

## Usage

Composition and fusion:

```julia
julia> using LandauGinzburgCategories, PolynomialRings

julia> @ring! ℚ[x,y,z]

julia> A = unit_matrix_factorization(x^2, x = y)
2×2 Arrax{ℚ[x,y,z],2}:
 0       x^2 + x*y + y^2
 -x + y  0

julia> B = unit_matrix_factorization(y^2, y = z)
2×2 Array{ℚ[x,y,z],2}:
 0       y^2 + y*z + z^2
 -y + z  0

julia> A ⨶ B
4×4 Array{ℚ[x,y,z],2}:
 0       0                   y^2 + y*z + z^2  x^2 + x*y + y^2
 0       0                   -x + y           y + -z         
 -y + z  x^2 + x*y + y^2     0                0              
 -x + y  -y^2 + -y*z + -z^2  0                0              


julia> fuse(A ⨶ B, :y)
2×2 Arrax{ℚ[x,y,z],2}:
 0       x^2 + x*z + z^2
 -x + z  0
```

Library of named potentials and of known orbifold equivalences between them:

```julia
julia> using LandauGinzburgCategories; LGLib = LandauGinzburgCategories.Library;
```

The `Aₙ`-series of potentials:

```
julia> LGLib.A₅()
x^6 + y^2

julia> LGLib.A₅(x, y)
x^6 + y^2

julia> LGLib.A(5, x, y)
x^6 + y^2
```
Exceptional unimodular singularities:
```
julia> LGLib.E₆(x, y)
x^3 + y^4
```

Et cetera.

Known orbifold equivalences:

```julia
julia> LGLib.orbifold_equivalence(LGLib.A5, LGLib.A2A2)
.....

```


## Status

This library has not been released yet and should therefore be considered alpha-quality software.

## Citation

If this library has been useful for your work, please cite it as https://arxiv.org/abs/1901.09019.

[travis-img]: https://travis-ci.org/tkluck/LandauGinzburgCategories.jl.svg?branch=master
[travis-url]: https://travis-ci.org/tkluck/LandauGinzburgCategories.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/4g6ax1ni7ijx3rn4?svg=true
[appveyor-url]: https://ci.appveyor.com/project/tkluck/landauginzburgcategories-jl

[codecov-img]: https://codecov.io/gh/tkluck/LandauGinzburgCategories.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/tkluck/LandauGinzburgCategories.jl
