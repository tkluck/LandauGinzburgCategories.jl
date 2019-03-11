module QuasiHomogeneous

import PolynomialRings: Polynomial, expansion

Gradings{I<:Integer} = NamedTuple{Names, NTuple{N, I}} where {Names, N}

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

function centralcharge(f, vars...)
    degs = find_quasihomogeneous_degrees(f, vars...)
    scale = 2//quasidegree(f, degs)
    return sum(1 - scale*d for d in values(degs))
end

end
