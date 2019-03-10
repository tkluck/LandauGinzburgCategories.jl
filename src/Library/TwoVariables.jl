module TwoVariables

import ...Library: Potential, subscript, leftvars, rightvars
import ...Library: orbifold_equivalence

struct A{N} <: Potential{2}
    A{n}()     where n = A{n}(leftvars(2)...)
    A{n}(x, y) where n = x^(n+1) - y^2
end
(::Type{A})(n, vars...) = A{n}(vars...)

for n in 2:30
    @eval const $( subscript(:A, n) ) = A{$n}
end

struct D{N} <: Potential{2}
    D{n}()     where n = D{n}(leftvars(2)...)
    D{n}(x, y) where n = x^(n-1) - x*y^2
end
(::Type{D})(n, vars...) = D{n}(vars...)

for n in 2:30
    @eval const $( subscript(:D, n) ) = D{$n}
end

struct E₆ <: Potential{2}
    E₆() = E₆(leftvars(2)...)
    E₆(x, y) = x^3 + y^4
end

struct E₇ <: Potential{2}
    E₇() = E₇(leftvars(2)...)
    E₇(x, y) = x^3 + x*y^3
end

struct E₈ <: Potential{2}
    E₈() = E₈(leftvars(2)...)
    E₈(x, y) = x^3 + y^5
end

struct A₂A₂ <: Potential{2}
    A₂A₂() = A₂A₂(leftvars(2)...)
    A₂A₂(x, y) = x^3 + y^3
end

end
