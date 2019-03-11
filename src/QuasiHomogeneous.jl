module QuasiHomogeneous

import PolynomialRings: Polynomial, expansion

Gradings{I<:Integer} = NamedTuple{Names, NTuple{N, I}} where {Names, N}


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
    exps = [e_i for (e,c) in expansion(f, vars...) for e_i in e]
    exps = reshape(exps, (length(vars), div(length(exps), length(vars))))'
    exps = exps // 1

    gradings = exps \ [1 for _=1:size(exps,1)]
    k = lcm(map(denominator, gradings)...)
    gradings = numerator.(k .* gradings)

    NamedTuple{vars}(tuple(gradings...))
end

function quasidegree(f::Polynomial, g::Gradings)
    iszero(f) && return -1
    maximum( sum(prod, zip(w, values(g))) for (w, p) in expansion(f, keys(g)...) )
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

end
