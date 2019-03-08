module Library

import PolynomialRings: @ring!

import ..Operations: unit_matrix_factorization

function leftvars(n)
    if n == 1
        @ring! ℤ[x]
        return (x,)
    elseif n == 2
        @ring! ℤ[x,y]
        return x, y
    elseif n == 3
        @ring! ℤ[x,y,z]
        return x, y, z
    else
        @ring! ℤ[x[1:n]]
        return x[1:n]
    end
end
function rightvars(n)
    if n == 1
        @ring! ℤ[u]
        return (u,)
    elseif n == 2
        @ring! ℤ[u, v]
        return u, v
    elseif n == 3
        @ring! ℤ[u, v, w]
        return u, v, w
    else
        @ring! ℤ[u[1:n]]
        return u[1:n]
    end
end

abstract type Potential{NumVars} end
numvars(::Type{Potential{NumVars}}) where NumVars = NumVars
numvars(f::Potential) = numvars(typeof(f))


module TwoVariables
    import ...Library: Potential

    struct A{N} <: Potential{2} end
    (::Type{A})(n, vars...) = A{n}()(vars...)
    (::A{n})(x, y) where n = x^(n+1) + y^2
    (Aₙ::A)() = Aₙ(leftvars(2)...)

    A₂ = A{2}(); A₃ = A{3}(); A₄ = A{4}(); A₅ = A{5}(); A₆ = A{6}()
    A₇ = A{7}(); A₈ = A{8}(); A₉ = A{9}()

    struct D{N} <: Potential{2} end
    (::Type{D})(n, vars...) = D{n}()(vars...)
    (::D{n})(x, y) where n = x^(n-1) + x*y^2
    (Dₙ::D)() = Dₙ(leftvars(2)...)

    D₂ = D{2}(); D₃ = D{3}(); D₄ = D{4}(); D₅ = D{5}(); D₆ = D{6}()
    D₇ = D{7}(); D₈ = D{8}(); D₉ = D{9}()

    struct _E6 <: Potential{2} end
    (::_E6)() = _E6(leftvars(2)...)
    (::_E6)(x, y) = x^3 + y^4
    E6 = E₆ = _E6() # yuck

    struct _E7 <: Potential{2} end
    (::_E7)() = _E7(leftvars(2)...)
    (::_E7)(x, y) = x^3 + x*y^3
    E7 = E₇ = _E7() # yuck

    struct _E8 <: Potential{2} end
    (::_E8)() = _E7(leftvars(2)...)
    (::_E8)(x, y) = x^3 + y^5
    E8 = E₈ = _E8() # yuck
end


function orbifold_equivalence(f::Potential, g::Potential)
    return orbifold_equivalence(f, g, leftvars(numvars(f)), rightvars(numvars(g)))
end

function orbifold_equivalence(f::Potential, g::Potential, left_vars, right_vars)
    error("No orbifold equivalence between $f and $g available in this library")
end

function orbifold_equivalence(::P, ::P, left_vars, right_vars) where P <: Potential
    # self-equivalence
    return unit_matrix_factorization(P(left_vars...), left_vars, right_vars)
end

end
