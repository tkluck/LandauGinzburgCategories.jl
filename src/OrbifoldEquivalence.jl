module OrbifoldEquivalence

import SparseArrays: spzeros
import LinearAlgebra: det, diagind, I

import PolynomialRings: base_extend, coefficient, gröbner_transformation
import PolynomialRings: constant_coefficient, flat_coefficients, Ideal
import PolynomialRings: expansion, formal_coefficients
import PolynomialRings.AbstractMonomials: any_divisor
import PolynomialRings.NamingSchemes: namingscheme

import ..Operations: supertrace, getpotential
import ..QuasiHomogeneous: find_quasihomogeneous_degrees, quasidegree, generic_quasihomogeneous_array, centralcharge

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

function materialize_ansatz(Q, W, vars...)
    shouldvanish = Q^2 - W*I

    iszero(shouldvanish) && return Q

    eqns = flat_coefficients(shouldvanish, vars...)
    P = eltype(eqns)
    J = Ideal(eqns)
    R = P/J

    iszero(one(R)) && return nothing
    return base_extend(Q, R)
end

"""
    for gradings in weight_split_criterion_gradings(W, N, vars...)
    end

Iterate over all grading matrices satisfying the weight split criterion
from Recknagel-Weinreb (2017).
"""
function weight_split_criterion_gradings(W, N, vars...)
    gr = find_quasihomogeneous_degrees(W, vars...)
    Wgrading = quasidegree(W, gr)

    monomials_in_W = [m for (m, c) in expansion(W, vars...)]

    gradings_of_divisors = map(monomials_in_W) do m
        gradings = Int[]
        # abuse any_divisor for just looping over the divisors
        any_divisor(m, namingscheme(vars...)) do divisor
            if !isone(divisor) && m != divisor
                d = sum(values(gr) .* divisor.e)
                push!(gradings, d)
            end
            false
        end
        gradings
    end

    admissible_rows = Set{Vector{Int}}()

    for divisor_grading_for_each_monomial in Iterators.product(gradings_of_divisors...)
        unique_divisor_gradings = Set(divisor_grading_for_each_monomial)
        if N == length(unique_divisor_gradings)
            divisor_gradings = sort(collect(unique_divisor_gradings))
            push!(admissible_rows, divisor_gradings)
        elseif N > length(unique_divisor_gradings)
            # FIXME: skipping this case for now
        else
            # not possitble
        end
    end

    Iterators.Filter(!isnothing, Base.Generator(admissible_rows) do row
        col = ntuple(i -> Wgrading - row[i], N)
        if all( (row .+ (col[k] - col[1])) in admissible_rows for k = 2:N)
            m = [ row[j] + (col[k] - col[1]) for j=1:N, k=1:N ]
            n = [ col[k] + (row[j] - row[1]) for j=1:N, k=1:N ]
            z = fill(-1, size(m))
            [ z m; n z]
        else
            nothing
        end
    end)
end

function search_orbifold_equivalence(f, g, left_vars, right_vars; max_rank=10)
    centralcharge(f, left_vars...) == centralcharge(g, right_vars...) || error("Cannot search orbifold equivalence if central charges disagree")
    W = g - f
    vgr = find_quasihomogeneous_degrees(W, left_vars..., right_vars...)
    for N in 1 : max_rank
        for gr in weight_split_criterion_gradings(W, N, left_vars..., right_vars...)
            next_coeff = formal_coefficients(typeof(W), :c)

            Q = generic_quasihomogeneous_array(gr, next_coeff)

            qdim1, qdim2 = quantum_dimensions(Q, W, left_vars, right_vars)

            if iszero(qdim1) || iszero(qdim2)
                @info("Found an admissible grading distribution, but its quantum dimension vanishes identically")
                continue
            end


            @info "Found an admissible grading distribution. Computing its Gröbner basis..."
            if (Q = materialize_ansatz(Q, W, left_vars..., right_vars...)) |> !isnothing
                return Q
            end
        end
    end
    return nothing
end

export supertrace, multivariate_residue, quantum_dimension, quantum_dimensions
export materialize_ansatz


end
