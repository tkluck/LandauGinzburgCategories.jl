module OrbifoldEquivalence

import SparseArrays: spzeros
import LinearAlgebra: det, diagind, I

import PolynomialRings: base_extend, coefficient, gröbner_transformation
import PolynomialRings: constant_coefficient, flat_coefficients, Ideal

import ..Operations: supertrace, getpotential

function multivariate_residue(g, f, vars...)
    R = base_extend(eltype(f))
    G, tr = gröbner_transformation(f)

    # TODO: compute that R/G is finite dimensional; otherwise, this computation
    # does not terminate
    M = spzeros(R, length(vars), length(f))
    for (row, v) in enumerate(vars)
        x = convert(R,v)
        factors, x_red = divrem(x, G)
        while !iszero(x_red)
            x = x^2
            factors, x_red = divrem(x, G)
        end
        M[row,:] = factors*tr
    end

    f_transformed = M * f
    g_transformed = g * det(M)

    term, r = divrem(prod(f_transformed), prod(convert(R, v) for v in vars))
    @assert(iszero(r))

    return coefficient(g_transformed, term, vars...)
end

"""
    lqdim, rqdim = quantum_dimensions(Q, left_vars, right_vars)

The left and right quantum dimensions of Q.
"""
function quantum_dimensions(Q::AbstractMatrix, left_vars, right_vars)
    W = getpotential(Q)
    return quantum_dimensions(Q, W, left_vars, right_vars)
end

"""
### A note on signs
The right quantum dimension of Q is, up to sign, equal to

    quantum_dimension(Q, right_vars, left_vars)

(note the swapped `left_vars` and `right_vars`). The difference in the explicit
signs at the front of the formula is

    (-1)^(binomial(n + 1, 2) - binomial(m + 1, 2)
    =
    (-1)^(n * (n + 1)/2 - m * (m + 1)/2)
    =
    (-1)^((n^2 + n - m^2 - m)/2)

As for the sign difference from the cyclic action inside the supertrace,
it is

    (-1)^(m * (n + m - 1))
    =
    (-1)^(m*n + m^2 - m)

Do they cancel? Not likely, as the first one depends on parity mod 4
whereas the second one depends on parity mod 2.
"""
function quantum_dimensions(Q::AbstractMatrix, W, left_vars, right_vars)
    lqdim = quantum_dimension(Q, W, left_vars, right_vars)
    rqdim = quantum_dimension(Q, W, right_vars, left_vars)

    m = length(left_vars)
    n = length(right_vars)
    ϵ₁ = (-1)^(m * (n + m - 1))
    ϵ₂ = (-1)^div(n^2 + n - m^2 - m, 2)

    return lqdim, ϵ₁*ϵ₂*rqdim
end

function quantum_dimension(Q::AbstractMatrix, W, left_vars, right_vars)
    g = supertrace(prod( diff(Q, v) for v in left_vars) * prod( diff(Q, v) for v in right_vars) )
    f = [diff(W,v) for v in left_vars]

    ϵ = iseven(binomial(length(left_vars)+1, 2)) ? 1 : -1
    return ϵ * constant_coefficient(multivariate_residue(g, f, left_vars...), right_vars...)
end

function materialize_ansatz(Q, f, vars...)
    shouldvanish = Q^2 - f*I

    iszero(shouldvanish) && return Q

    eqns = flat_coefficients(shouldvanish, vars...)
    P = eltype(eqns)
    J = Ideal(eqns)
    R = P/J

    if iszero(one(R))
        error("Cannot materialize ansatz: it yields no solutions")
    end

    return base_extend(Q, R)
end

export supertrace, multivariate_residue, quantum_dimension, quantum_dimensions
export materialize_ansatz


end
