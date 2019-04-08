module Library

import LinearAlgebra: I

import PolynomialRings: @ring!, @ring, polynomial_ring

import ..QuasiHomogeneous: centralcharge
import ..Operations: unit_matrix_factorization, dual

"""
    Potential{NumVars}

A type representing a potential through a classification scheme.
For example, the potential ``x^4 - y^2`` is called ``A_3`` through
the ADE classification of unimodular potentials, and it is
represented in Julia by the type `TwoVariables.A₃`. This is a subtype
of `Potential{2}`.
"""
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
        if centralcharge(f) != centralcharge(g)
            res = false
        end
    end
    if !ismissing(res) && res !== false
        res_sq = res^2
        target = eltype(res_sq)(g(right_vars...) - f(left_vars...))
        if res_sq != target * I
            @error "Library contains invalid orbifold equivalence for $f and $g. Please open an issue at https://github.com/tkluck/LandauGinzburgCategories.jl/issues/new" res_sq[1,1] target
            error()
        end
    end

    return res
end

centralcharge(::Type{P}) where P <: Potential = centralcharge(P(leftvars(P)...), leftsyms(P)...)

function _orbifold_equivalence_def(::Type{P}, ::Type{P}, left_vars, right_vars) where P <: Potential
    # self-equivalence
    n = numvars(P)
    xsyms = [Symbol("x", i) for i in 1:n] # TODO: needs gensym() or some such
    _, xvals = polynomial_ring(xsyms...)
    res = unit_matrix_factorization(P(xvals...); (xsyms .=> right_vars)...)
    res = res(; (xsyms .=> left_vars)...)
end

include("Library/TwoVariables.jl")
include("Library/ThreeVariables.jl")
include("Library/CarquevilleRunkel2012.jl")
include("Library/RecknagelWeinreb2017.jl")
include("Library/NewtonRosCamacho2015.jl")
include("Library/Interactive.jl")

export orbifold_equivalence
export TwoVariables, ThreeVariables
export choose_equivalence

end
