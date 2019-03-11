module Library

import LinearAlgebra: I

import PolynomialRings: @ring!, @ring, polynomial_ring

import ..QuasiHomogeneous: centralcharge
import ..Operations: unit_matrix_factorization, dual

abstract type Potential{NumVars} end
numvars(::Type{<:Potential{NumVars}}) where NumVars = NumVars
numvars(f::Potential) = numvars(typeof(f))

const _subscripts = Dict(0:9 .=> collect("₀₁₂₃₄₅₆₇₈₉"))
subscript(a) = Symbol(a)
subscript(a...) = Symbol(join(map(subscript, a)))
subscript(n::Integer) = Symbol(reverse(join(_subscripts[i] for i in digits(n))))

leftsyms(f::Type{<:Potential})  = leftsyms(numvars(f))
rightsyms(f::Type{<:Potential}) = rightsyms(numvars(f))

function leftsyms(n::Integer)
    n == 1 && return (:x,)
    n == 2 && return (:x, :y)
    n == 3 && return (:x, :y, :z)
    error("Not implemented")
end

function rightsyms(n::Integer)
    n == 1 && return (:u,)
    n == 2 && return (:u, :v)
    n == 3 && return (:u, :v, :w)
    error("Not implemented")
end

function leftvars(spec)
    _, vars = polynomial_ring(leftsyms(spec)...)
    return vars
end

function rightvars(spec)
    _, vars = polynomial_ring(rightsyms(spec)...)
    return vars
end



"""
    orbifold_equivalence(f, g)

Return a matrix representing an orbifold equivalence between `f` and `g`,
if one is available in the library. Return `false` if it is known that
f and g are not orbifold equivalent. Return `missing` if this is not known
by this library.
"""
function orbifold_equivalence end
function _orbifold_equivalence_def end

_orbifold_equivalence_def(f, g, left_vars, right_vars) = missing

function orbifold_equivalence(f::Type{<:Potential}, g::Type{<:Potential})
    return orbifold_equivalence(f, g, leftvars(numvars(f)), rightvars(numvars(g)))
end

function orbifold_equivalence(f::Type{<:Potential}, g::Type{<:Potential}, left_vars, right_vars)
    res = _orbifold_equivalence_def(f, g, left_vars, right_vars)
    if ismissing(res)
        res = _orbifold_equivalence_def(g, f, right_vars, left_vars)
        if !ismissing(res) && res !== false
            res = dual(res)
        end
    end
    if ismissing(res)
        if centralcharge(f(), leftsyms(f)...) != centralcharge(g(), leftsyms(g)...)
            res = false
        end
    end
    @assert ismissing(res) || res === false ||
            res^2 == (g(right_vars...) - f(left_vars...))*I "Library contains invalid orbifold equivalence for $f and $g. Please open an issue at github.com/tkluck/LandauGinzburgCategory.jl"

    return res
end

function _orbifold_equivalence_def(::Type{P}, ::Type{P}, left_vars, right_vars) where P <: Potential
    # self-equivalence
    n = numvars(P)
    xsyms = [Symbol("x", i) for i in 1:n]
    _, xvals = polynomial_ring(xsyms...)
    res = unit_matrix_factorization(P(xvals...); (xsyms .=> right_vars)...)
    res = res(; (xsyms .=> left_vars)...)
end

include("Library/TwoVariables.jl")
include("Library/CarquevilleRunkel2012.jl")
include("Library/RecknagelWeinreb2017.jl")

export orbifold_equivalence
export TwoVariables

end
