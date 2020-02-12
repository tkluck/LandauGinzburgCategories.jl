
"""
Orbifold equivalences from

> Timo Kluck, Ana Ros Camacho
> "Computational aspects of orbifold equivalence"
> https://arxiv.org/pdf/1901.09019.pdf
"""

import LinearAlgebra: det, I
import SparseArrays: spzeros

import PolynomialRings: @flat_coefficients, deg, max_variable_index, linear_coefficients
import PolynomialRings: @polyvar
import PolynomialRings.Util.LinAlgUtil: echelon
import PolynomialRings.NamingSchemes: @namingscheme

import ..QuasiHomogeneous: find_quasihomogeneous_degrees, quasidegree, forgradedmonomials, Gradings, complete_quasihomogeneous_polynomial
import ..OrbifoldEquivalence: quantum_dimensions

complete_with_generic_monomials(Q, g::Gradings, next_coeff) = complete_quasihomogeneous_polynomial.(Q, Ref(g), Ref(next_coeff))

function adjugate(M)
    cofactor(M, i, j) = det(@view M[[1:i-1; i+1:end], [1:j-1;j+1:end]])
    return [
        (-1)^(i+j) * cofactor(M, j, i)
        for i in axes(M, 2), j in axes(M, 1)
    ]
end

function _nullspace(M)
    M_aug = echelon(M)
    zero_cols = findall(j -> iszero(@view M_aug[1:size(M, 1), j]), axes(M_aug, 2))
    return M_aug[(size(M, 1) + 1):end, zero_cols]
end

function reduce_linear_equations(eqns, namingscheme)
    eqns1 = filter(f -> 1 == deg(f, namingscheme), unique(eqns))

    n = max_variable_index(eqns1)
    m = length(eqns1)
    M = spzeros(Int, m, n)

    for (i, eqn) in enumerate(eqns1)
        v = linear_coefficients(eqn, namingscheme)
        M[i, 1:length(v)] .= v
    end

    K = _nullspace(M)
    @show size(K)
    u, v = size(K)

    @polyvar d[]
    substitutions = K * [d[i] for i in axes(K, 2)]
    return filter(!iszero, eqns(c = i -> i <= u ? substitutions[i] : d[i - u + v + 1]))
end

function _orbifold_equivalence_def(f::Type{ThreeVariables.E₁₈}, g::Type{ThreeVariables.Q₁₂{:v2}}) #, left_vars, right_vars)
    @ring! Int[c[]]
    @ring! ℚ[x,y,z,u,v,w]
    gr = (x=3, y=10, z=15, u=6, v=10, w=12)

    dQ♯ = [
                      z     0  y - v              u^2
                      0     z   -u^3 -v^2 - y*v - y^2
        v^2 + y*v + y^2   u^2     -z                0
                   -u^3 v - y      0               -z
    ]

    function safediv(f, g)
        d, r = divrem(f, g)
        @assert iszero(r)
        return d
    end
    dQ♭ = safediv.(adjugate(dQ♯), u^5 + v^3 - y^3 - z^2)
    dQ = [
        zero(dQ♯)        dQ♯
              dQ♭  zero(dQ♭)
    ]
    dQ = complete_with_generic_monomials(dQ, gr, c)
    c_l, c_r = c(), c()

    W = g(u,v,w) - f(x,y,z)
    eqns = @flat_coefficients dQ^2 - W*I x y z u v w
    q_l, q_r = quantum_dimensions(dQ, W, [:x,:y,:z], [:u,:v,:w])
    push!(eqns, 1 - c_l * q_l)
    push!(eqns, 1 - c_r * q_r)

    return reduce_linear_equations(eqns, @namingscheme c[])
end

function _orbifold_equivalence_def(f::Type{ThreeVariables.Q₁₈}, g::Type{ThreeVariables.E₃₀}) #, left_vars, right_vars)
    @ring! Int[c[]]
    @ring! ℚ[x,y,z,u,v,w]
    gr = (x=6, y=16, z=21, u=3, v=16, w=24)

    dQ♯ = [
              z v^2 + y*v + y^2 x^4 + w                 0
          y - v            -x*z       0           x^4 + w
        x^4 - w               0    -x*z -(v^2 + y*v + y^2)
              0         x^4 - w   v - y                 z
    ]

    function safediv(f, g)
        d, r = divrem(f, g)
        @assert iszero(r)
        return d
    end
    dQ♭ = safediv.(adjugate(dQ♯), v^3 + w^2 - x^8 - y^3 - x*z^2)
    dQ = [
        zero(dQ♯)        dQ♯
              dQ♭  zero(dQ♭)
    ]
    dQ = complete_with_generic_monomials(dQ, gr, c)
    c_l, c_r = c(), c()

    W = g(u,v,w) - f(x,y,z)
    eqns = @flat_coefficients dQ^2 - W*I x y z u v w
    q_l, q_r = quantum_dimensions(dQ, W, [:x,:y,:z], [:u,:v,:w])
    push!(eqns, 1 - c_l * q_l)
    push!(eqns, 1 - c_r * q_r)

    return reduce_linear_equations(eqns, @namingscheme c[])
end
