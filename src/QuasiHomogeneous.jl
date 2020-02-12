module QuasiHomogeneous

import PolynomialRings: Polynomial, expand, monomialtype, polynomial_ring

Gradings{I<:Integer} = NamedTuple{Names, NTuple{N, I}} where {Names, N}

function forgradedmonomials(f, total_grading, g::Gradings)
    M = monomialtype(keys(g)...)
    f′(e...) = f(M(e))
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
```
julia> find_quasihomogeneous_degrees(x^4*y + x*y^9, :x, :y)
(x = 8, y = 3)
```
"""
function find_quasihomogeneous_degrees(f::Polynomial, vars::Symbol...)
    exps = [e_i for (e,c) in expand(f, vars...) for e_i in e]
    exps = reshape(exps, (length(vars), div(length(exps), length(vars))))'
    exps = exps // 1

    gradings = exps \ [1 for _=1:size(exps,1)]
    k = lcm(map(denominator, gradings)...)
    gradings = numerator.(k .* gradings)

    NamedTuple{vars}(tuple(gradings...))
end

function quasidegree(f::Polynomial, g::Gradings)
    iszero(f) && return -1
    maximum( sum(prod, zip(w, values(g))) for (w, p) in expand(f, keys(g)...) )
end

"""
    c = centralcharge(f, vars...)

Return the central charge of `f` with respect to the variables `vars`.

This is the value of the expression

``\\sum_{i=1}^n 1 - q_i``

where ``q_i`` is the grading of the ``i``th variable under
a (ℚ-valued) grading for which `f` is homogeneous of degree 2.
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
    R,_ = polynomial_ring(keys(gr)..., basering=Int)
    return complete_quasihomogeneous_polynomial(zero(R), gr, next_coeff, grade)
end

function generic_quasihomogeneous_array(gradings::Array{<:Integer}, gr::Gradings, next_coeff)
    return [ generic_quasihomogeneous_polynomial(grade, gr, next_coeff) for grade in gradings ]
end

end
