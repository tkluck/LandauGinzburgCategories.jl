module OrbifoldEquivalence

import SparseArrays: spzeros
import LinearAlgebra: det, diagind, I

import Combinatorics: with_replacement_combinations

import PolynomialRings: base_extend, coefficient, gröbner_transformation
import PolynomialRings: constant_coefficient, flat_coefficients, Ideal
import PolynomialRings: expansion, formal_coefficients, xdivrem, ofminring, minring
import PolynomialRings: to_dense_monomials
import PolynomialRings.AbstractMonomials: any_divisor
import PolynomialRings.Expansions: expansionorder
import PolynomialRings.NamingSchemes: namingscheme, @namingscheme

import ..Operations: supertrace, getpotential
import ..QuasiHomogeneous: Gradings, find_quasihomogeneous_degrees, quasidegree, generic_quasihomogeneous_array, centralcharge, forgradedmonomials

"""
    r = multivariate_residue(g, f, vars...)

Compute the multivariate residue

``r = \\mathrm{res} \\left ( \\frac{g \\mathrm{d}x_1 \\wedge \\cdots \\wedge \\mathrm{d}x_n }{f_1,\\cdots,f_n}\\right )``

where `n = length(f)` and `vars = [x1, ...., xn]`. See

> Lipman, Joseph, Residues and traces of differential forms via Hochschild
homology vol. 61, (American Mathematical Soc., 1987).

# Example
```jldoctest
julia> using PolynomialRings, LandauGinzburgCategories;

julia> @ring! Int[x,y];

julia> multivariate_residue(3x^2*y^2, [x^3, y^3], x, y)
3

"""
function multivariate_residue(g, f, vars...)
    R = mapreduce(minring, promote_type, f)
    R = base_extend(R)
    f = R.(f)
    vars = R.(vars)
    G, tr = gröbner_transformation(f)

    # TODO: compute that R/G is finite dimensional; otherwise, this computation
    # does not terminate
    M = spzeros(R, length(vars), length(f))
    for (row, x) in enumerate(vars)
        factors, x_red = divrem(x, G)
        while !iszero(x_red)
            x = x^2
            factors, x_red = divrem(x, G)
        end
        M[row, :] = factors * tr
    end

    f_transformed = M * f
    g_transformed = g * det(M)

    term, r = divrem(prod(f_transformed), prod(vars))
    @assert(iszero(r))

    return coefficient(g_transformed, term, vars...)
end

"""
    lqdim, rqdim = quantum_dimensions(Q, left_vars, right_vars, W=getpotential(Q))

The left and right quantum dimensions of Q.

# Examples
```jldoctest
julia> using PolynomialRings, LandauGinzburgCategories;

julia> @ring! Int[x,y,u,v];

julia> Q = unit_matrix_factorization(x^5 - y^2, x => u, y => v);

julia> quantum_dimensions(Q, (x,y), (u,v))
(1, -1)

"""
function quantum_dimensions(Q::AbstractMatrix, left_vars, right_vars, W = getpotential(Q))
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

function materialize_ansatz(Q, W, left_vars, right_vars, extracoeff)
    shouldvanish = Q^2 - W*I

    qdim1, qdim2 = quantum_dimensions(Q, left_vars, right_vars, W)

    if iszero(shouldvanish)
        if !iszero(qdim1) && !iszero(qdim2)
            return Q
        else
            return nothing
        end
    end

    qdim_eqn = ofminring(qdim1 * qdim2 * extracoeff - 1)
    other_eqns = flat_coefficients(shouldvanish, left_vars..., right_vars...)

    eqns = [qdim_eqn; other_eqns]
    P = eltype(eqns)
    J = Ideal(eqns)
    R = P/J

    iszero(one(R)) && return nothing
    return base_extend(Q, R)
end

"""
    for gradings in weight_split_criterion_gradings(W, N, vars...)

Iterate over all grading matrices satisfying the weight split criterion
from Recknagel-Weinreb (2017).
"""
function weight_split_criterion_gradings(W, N, vars...)
    gr = find_quasihomogeneous_degrees(W, vars...)
    Wgrading = quasidegree(W, gr)

    if isodd(Wgrading)
        # ensure the quasidegree of a matrix factorization of
        # W is an integer as well.
        gr = Gradings(namingscheme(gr), (2 .* values(gr))...)
        Wgrading = quasidegree(W, gr)
    end

    monomials_in_W = [m for (m, c) in expansion(W, vars...)]

    gradings_of_divisors = map(monomials_in_W) do m
        gradings = Int[]
        # abuse any_divisor for just looping over the divisors
        any_divisor(m, namingscheme(expansionorder(vars...))) do divisor
            if !isone(divisor) && m != divisor
                d = quasidegree(divisor, gr)
                push!(gradings, d)
            end
            false
        end
        gradings
    end

    max_grading_diff_subsequent_rows = Wgrading ÷ 2
    nontrivial_grades = filter(1 : N * max_grading_diff_subsequent_rows) do i
        any_graded_monomial = false
        forgradedmonomials(i, gr) do _
            any_graded_monomial = true
        end
        any_graded_monomial
    end

    admissible_rows = Vector{Vector{Int}}()
    for c in with_replacement_combinations(nontrivial_grades, N)
        if all(any(x in c for x in g) for g in gradings_of_divisors)
            push!(admissible_rows, c)
        end
    end

    Iterators.Filter(!isnothing, Iterators.flatten(Base.Generator(admissible_rows) do row
        Base.Generator(admissible_rows) do col
            col = reverse(col)
            if all( (row .+ (col[k] - col[1])) in admissible_rows for k = 2:N)
                m = [ row[j] + (col[k] - col[1]) for j=1:N, k=1:N ]
                n = Wgrading .- transpose(m)
                z = fill(-1, size(m))
                ([ z m; n z], gr)
            else
                nothing
            end
        end
    end))
end

"""
    search_orbifold_equivalence(f, g, left_vars, right_vars; max_rank=10)

Perform a brute-force search for an orbifold equivalence between f and g.

This uses the algorithm made popular by

> Recknagel, Andreas and Weinreb, Paul, "Orbifold equivalence: structure and
> new examples", arXiv preprint arXiv:1708.08359 (2017).

# Example
```
julia> @ring! Int[x,y];

julia> search_orbifold_equivalence(x^3, 2y^3, (x,), (y,))

"""
function search_orbifold_equivalence(f, g, left_vars, right_vars; max_rank=10)
    if centralcharge(f, left_vars...) != centralcharge(g, right_vars...)
        @info "Central charges disagree; no orbifold equivalence exists"
        return nothing
    end
    W = g - f
    for N in 1 : max_rank
        for (grading_matrix, gr) in weight_split_criterion_gradings(W, N, left_vars..., right_vars...)
            next_coeff = formal_coefficients(typeof(W), :c)
            c1 = next_coeff()

            Q = generic_quasihomogeneous_array(grading_matrix, gr, next_coeff)

            qdim1, qdim2 = quantum_dimensions(Q, left_vars, right_vars, W)
            if iszero(qdim1) || iszero(qdim2)
                @info("Found an admissible grading distribution, but its quantum dimension vanishes identically")
                continue
            end

            @info "Found an admissible grading distribution. See if we can materialize it..." Q
            if (Q = materialize_ansatz(Q, W, left_vars, right_vars, c1)) |> !isnothing
                @info "Done"
                return Q
            end
            @info "We can't."
        end
    end
    return nothing
end

export supertrace, multivariate_residue, quantum_dimension, quantum_dimensions
export materialize_ansatz


end
