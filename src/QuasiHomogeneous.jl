module QuasiHomogeneous

import PolynomialRings: Polynomial, AbstractMonomial, expand, monomialtype, polynomial_ring
import PolynomialRings.Expansions: expansionorder
import PolynomialRings.NamingSchemes: NamingScheme, Named, namingscheme, num_variables, variablesymbols

struct Gradings{Scheme <: NamingScheme, N, I<:Integer}
    scheme :: Scheme
    grades :: NTuple{N, I}
end

namingscheme(g::Gradings) = g.scheme
Base.values(g::Gradings) = g.grades
function Base.show(io::IO, g::Gradings)
    print(io, "Gradings(", g.scheme, ", ")
    join(io, g.grades, ", ")
    print(io, ")")
end
function Base.show(io::IO, g::Gradings{<:Named})
    ntup = NamedTuple{variablesymbols(namingscheme(g))}(values(g))
    print(io, "Gradings", ntup)
end

function Gradings(scheme::NamingScheme, grades::Integer...)
    num_variables(scheme) == length(grades) || error("Gradings: need as many gradings as variables")
    return Gradings(scheme, promote(grades...))
end

function Gradings(; kwds...)
    scheme = namingscheme(keys(kwds.data)...)
    return Gradings(scheme, values(kwds.data)...)
end

function forgradedmonomials(f, total_grading, g::Gradings)
    M = monomialtype(namingscheme(g))
    f′(e...) = f(exp(M, e))
    _forgradedmonomials(f′, total_grading, values(g))
end

function _forgradedmonomials(f, total_grading, g::NTuple{N, <:Integer}) where N
    total_grading >= 0 || return
    if length(g) < 1
        return
    elseif length(g) == 1
        if total_grading % g[1] == 0
            f(total_grading ÷ g[1])
        end
    else
        for exp in 0:div(total_grading, g[1])
            remaining = total_grading - g[1] * exp
            _forgradedmonomials(remaining, g[2:end]) do (e...)
                f(exp, e...)
            end
        end
    end
end

"""
    gradings = find_quasihomogeneous_degrees(f, vars...)

Find the assignment of gradings to the variables `vars` such that `f` is a
quasi-homogeneous polynomial with respect to these gradings, if such an
assignment exists, and return it. Otherwise, raise an exception.

# Example
```jldoctest
julia> find_quasihomogeneous_degrees(x^4*y + x*y^9, :x, :y)
Gradings(x=8, y=3)
```
"""
function find_quasihomogeneous_degrees(f::Polynomial, vars...)
    exps = [e_i for (e,c) in expand(f, vars...) for e_i in e]
    exps = reshape(exps, (length(vars), div(length(exps), length(vars))))'
    exps = exps // 1

    gradings = exps \ [1 for _=1:size(exps,1)]
    k = lcm(map(denominator, gradings)...)
    gradings = numerator.(k .* gradings)

    Gradings(namingscheme(expansionorder(vars...)), gradings...)
end

"""
    d = quasidegree(f::Polynomial, g::Gradings)

The quasidegree of `f`, with variable gradings specified by `g`.

# Example
```jldoctest
julia> using PolynomialRings, LandauGinzburgCategories

julia> @ring! Int[x, y];

julia> quasidegree(x^2 * y^3, (x=2, y=1))
7
```
"""
function quasidegree(f::Polynomial, g::Gradings)
    iszero(f) && return -1
    maximum( sum(prod, zip(w, values(g))) for (w, p) in expand(f, namingscheme(g)) )
end

quasidegree(f::AbstractMonomial, g::Gradings) = sum(exponents(f, namingscheme(g)) .* values(g))

"""
    c = centralcharge(f, vars...)

Return the central charge of `f` with respect to the variables `vars`.

This is the value of the expression

``\\sum_{i=1}^n 1 - q_i``

where ``q_i`` is the grading of the ``i``th variable under
a (ℚ-valued) grading for which `f` is homogeneous of degree 2.

# Example
```jldoctest
julia> using PolynomialRings, LandauGinzburgCategories

julia> @ring! Int[x, y];

julia> centralcharge(x^5 + y^2, x, y)
7
```
"""
function centralcharge(f, vars...)
    degs = find_quasihomogeneous_degrees(f, vars...)
    scale = 2//quasidegree(f, degs)
    return sum(1 - scale*d for d in values(degs))
end

function complete_quasihomogeneous_polynomial(f, gr::Gradings, next_coeff, d = quasidegree(f, gr))
    P = promote_type(typeof(f), eltype(next_coeff))
    result = P(f)
    forgradedmonomials(d, gr) do m
        if iszero(f[m])
            result += next_coeff() * m
        end
    end
    return result
end

function generic_quasihomogeneous_polynomial(grade::Integer, gr::Gradings, next_coeff)
    R = polynomial_ring(namingscheme(gr), basering=Int)
    return complete_quasihomogeneous_polynomial(zero(R), gr, next_coeff, grade)
end

function generic_quasihomogeneous_array(gradings::Array{<:Integer}, gr::Gradings, next_coeff)
    return [ generic_quasihomogeneous_polynomial(grade, gr, next_coeff) for grade in gradings ]
end

end
