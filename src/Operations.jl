module Operations

import Combinatorics: permutations, parity
import LinearAlgebra: I, diagind, UpperTriangular
import SparseArrays: sparse, spzeros, issparse

import PolynomialRings: Polynomial, polynomial_ring
import PolynomialRings: constant_coefficient, gröbner_transformation
import PolynomialRings: map_coefficients, base_extend
import PolynomialRings: ofminring, minring
import PolynomialRings: gröbner_basis, lift, syzygies
import PolynomialRings.Expansions: expansiontypes, expand
import PolynomialRings.QuotientRings: QuotientRing
import PolynomialRings.MonomialOrderings: MonomialOrder
import PolynomialRings.NamingSchemes: Named
import PolynomialRings.Reductions: mingenerators
import PolynomialRings.Solve: matrix_solve_affine
import PolynomialRings.Util: nzpairs, @showprogress

import ..MatrixUtil: sweepscalars!, triangularperm, exactinv, exactsqrt, columns

# for clarity when calling collect() for its "side-effect"
# of returning a dense vector/matrix.
const todense = collect

"""
    M = flatten_blocks(X)

Construct a matrix from a matrix of matrices by concatenating
them horizontally and vertically.
"""
flatten_blocks(X) = vcat([hcat(X[i,:]...) for i=1:size(X,1)]...)

function from_alternating_grades(M::AbstractMatrix)
    rows, cols = size(M)
    rowperm = [1:2:rows; 2:2:rows]
    colperm = [1:2:cols; 2:2:cols]
    M[rowperm, colperm]
end

function to_alternating_grades(M::AbstractMatrix)
    rows, cols = size(M)
    rowperm = [1:2:rows; 2:2:rows]
    colperm = [1:2:cols; 2:2:cols]
    M[invperm(rowperm), invperm(colperm)]
end

"""
    A⨷B

Graded tensor product of ℤ/2 graded block matrices. The result consists
of even/odd blocks just like the input and the operation is associative.

_Graded_ means that we introduce Koszul signs. Writing ``A = a^i_j e_i e^j`` and
``B = b^k_ℓ f_k f^ℓ`` for some (co)basis vectors, we represent the tensor
product as

``A⨷B = (-1)^{|j|(|k|+|ℓ|)} p^i_j q^k_ℓ e_i f_k f^ℓ e^j``

where the sign comes from commuting ``e^j`` with ``f_k f^ℓ``. This ensures that,
upon matrix multiplication between two such tensor products, the contractions
between ``e^j`` and ``e_j``, and between ``f^ℓ`` and ``f_ℓ``, are already
adjacent and do not introduce more signs.

Mapping the multi-indices ``(i,k)`` to rows and ``(ℓ,j)`` to columns requires
some care to ensure that we end up with even/odd blocks. A comment in the code
explains how this is done. We first interlace the rows/columns to get alternating
grades in rows/columns (`to_alternating_grades` and `from_alternating_grades`).
Then we need only to reverse some row/column orders. Finally, we de-interlace
to get even/odd blocks. (I wonder if there is a way to avoid making the detour
through the alternating signs representation, but it seems hard to maintain
associativity.)
"""
function ⨷(A,B)
    n,m = size(A)
    k,l = size(B)
    all(iseven, (n,m,k,l)) || throw(ValueError("Need ℤ/2 graded matrices, even rank and corank"))

    A,B = to_alternating_grades.((A,B))

    res = [a*B for a in A]

    for row in axes(res,1), col in axes(res,2)
        inner = res[row,col]
        for inner_row in axes(inner,1), inner_col in axes(inner,2)
            # Koszul signs; see formula in the documentation string
            if iseven(col) && isodd(inner_row + inner_col)
                inner[inner_row, inner_col] *= -1
            end
        end
        # make sure gradings alternate in the output matrix: they alternate
        # in both inner and outer, and the result grading is their sum. At a
        # boundary between two 'inner' matrices, both signs change, so the
        # adjacent rows/columns have the same sign. We can fix this by reversing
        # the rows/columns within an inner matrix, depending on the sign of the
        # outer row/column.
        if iseven(row)
            inner[:,:] = inner[end:-1:1,:]
        end
        if iseven(col)
            inner[:,:] = inner[:,end:-1:1]
        end
    end
    return from_alternating_grades(flatten_blocks(res))
end

"""
    t = supertrace(Q::AbstractMatrix)

Trace of a of ℤ/2-graded block matrix with the Koszul sign convention.
"""
function supertrace(Q::AbstractMatrix)
    n,m = size(Q)
    (n == m && iseven(n)) || throw(ArgumentError("Cannot compute supertrace of $n x $m matrix"))

    firsthalf = diagind(Q)[1:end÷2]
    secondhalf = diagind(Q)[end÷2+1:end]
    return sum(i->Q[i], firsthalf) - sum(i->Q[i], secondhalf)
end

"""
    A⨶B

Tensor product of matrix factorizations.

This is the composition operation in the bicategory or
Landau Ginzburg models.
"""
⨶(A,B) = A⨷one(B) + one(A)⨷B

"""
    X⊗quotient_ring

Matrix-representation of the operator obtained from ``X``, a matrix acting on a
module over a ring ``R``, by tensoring (over ``R``) with a quotient of ``R``.

Note that `eltype(X)` is not necessarily equal to ``R``; it may also be
an ``R``-algebra. For example, ``R=k[y]`` and `eltype(X) == @ring(k[x,y])` works.
"""
function ⊗(X::AbstractMatrix{<:Polynomial}, quotient_ring::Type{<:QuotientRing})
    X_inflated = representation_matrix.(quotient_ring, X)
    flatten_blocks(X_inflated)
end


"""
    Q,ϵ = ⨶(A,B,W,vars...)

Finite-rank homotopy representation of ``A⨶B``, where we remove the variables
`vars` from the result.

We return a matrix factorization ``Q`` together with an idempotent ϵ representing
the direct summand of ``Q`` that represents ``A⨶B``.

See the pushforward paper by Dykerhoff&Murfet.
"""
function ⨶(A,B,W,vars...)
    R,_ = polynomial_ring(vars...;basering=basering(W))
    ∇W = diff.(W, vars)
    Jacobian = R/Ideal(∇W...)

    Q = (A⨶B) ⊗ Jacobian
end

"""
    unit_matrix_factorization(f; source_to_target...)

A ℤ/2-graded matrix that squares to `f(;source_to_target...) - f` times
the identity matrix.

The source for this formulation is

> Adjunctions and defects in Landau-Ginzburg models, Nils Carqueville and Daniel Murfet
"""
function unit_matrix_factorization(f; source_to_target...)
    source_to_target = collect(source_to_target)
    function ∂(f, n)
        f₊ = f(; source_to_target[1:n - 1]...)
        f₋ = f(; source_to_target[1:n    ]...)
        x, y = source_to_target[n]
        num = f₊ - f₋
        return div(num, typeof(num)(x) - typeof(num)(y))
    end

    # x represents, through its bit-representation, a basis element of the exterior
    # algebra. To be precise, x represents the element theta_i_1 \wedge ... \wedge theta_i_n
    # where i_1 ... i_n are the bits set in x.
    #
    # The use of 'gray code' (see wikipedia) ensures that subsequent elements differ by
    # exactly one bit. This way, rows/columns of our result matrix have _alternating_ signs.
    N = length(source_to_target)
    R = promote_type(typeof(f), typeof(f(;source_to_target...)))
    gray_code(x) = xor(x, x>>1)
    permutation = map(n->gray_code(n)+1, 0:2^N-1)
    inv_perm = invperm(permutation)
    to_index(x) = inv_perm[x+1]

    function wedge_product_matrix(T, i)
        result = zeros(T, 2^N,2^N)
        for j in 0:2^N-1
            j&(1<<(i-1)) != 0 && continue

            k = j | (1 << (i-1))
            sign = (-1)^count_ones(j & (1<<i - 1))
            result[to_index(j), to_index(k)] = T(sign)
        end
        return result
    end

    function lift_matrix(T, i)
        result = zeros(T, 2^N, 2^N)
        for j in 0:2^N-1
            j&(1<<(i-1)) == 0 && continue

            k = j & ~(1<<(i-1))
            sign = (-1)^count_ones(j & (1<<i - 1))
            result[to_index(j), to_index(k)] = T(sign)
        end
        return result
    end

    delta_plus = sum(∂(f, i) * wedge_product_matrix(R, i)
                     for i=1:N)
    delta_minus = sum((R(x) - R(y)) * lift_matrix(R, i)
                      for (i,(x,y)) in enumerate(source_to_target))
    return from_alternating_grades(delta_plus + delta_minus)
end

"""
    D, A = block_diagonalization(X)

Decompose the matrix factorization X into a direct sum of irreducible
matrix factorizations, and represent this direct sum as a block-diagonal
matrix factorization `D` such that ``A^{-1} D A = X``.

NOTE: for now, this function only returns D, and not yet A!
"""
function block_diagonalization(X)
    D = copy(X)
    A = one(X)
    top_right = @view D[1:end÷2, end÷2+1:end]

    # the following functions take parameters indexing `top_right`, but they
    # operate simultaneously on the bottom left of `D` as well. This happens
    # in such a way that X^2 remains the same.
    function rowop(i, factor, j)
        D[i,:] += factor * D[j,:]
        D[:,j] -= factor * D[:,i]
    end
    function colop(i, factor, j)
        D[:,end÷2 + i] += factor * D[:,end÷2 + j]
        D[end÷2 + j,:] -= factor * D[end÷2 + i,:]
    end

    for _=1:100
        for row in axes(top_right,1), col in axes(top_right,2)
            for row2 in axes(top_right,1)
                row2 == row && continue
                iszero(top_right[row2,col]) && continue
                (d,),r = divrem(top_right[row2,col], [top_right[row,col]])
                if iszero(r)
                    rowop(row2, -d, row)
                    @goto more_loops
                end
            end
            for col2 in axes(top_right,2)
                col2 == col && continue
                iszero(top_right[row,col2]) && continue
                (d,),r = divrem(top_right[row,col2], [top_right[row,col]])
                if iszero(r)
                    colop(col2, -d, col)
                    @goto more_loops
                end
            end
        end
        break
        @label more_loops
    end

    # most inefficient algorithm I can think of
    blocks = []
    for row in axes(top_right,1), col in axes(top_right,2)
        if !iszero(top_right[row, col])
            push!(blocks, (Set([row]), Set([col])))
        end
    end
    for i in 1:length(blocks)
        j = i+1
        while j <= length(blocks)
            if !isempty(blocks[i][1] ∩ blocks[j][1]) || !isempty(blocks[i][2] ∩ blocks[j][2])
                blocks[i] = (blocks[i][1] ∪ blocks[j][1], blocks[i][2] ∪ blocks[j][2])
                deleteat!(blocks, j)
                j = i+1
            else
                j += 1
            end
        end
    end

    rowperm = vcat(sort.(collect.(getindex.(blocks, 1)))...)
    colperm = vcat(sort.(collect.(getindex.(blocks, 2)))...)
    D[1:end÷2,:] = D[rowperm,:]
    D[:,end÷2+1:end] = D[:,colperm .+ end÷2]

    D[end÷2+1:end,:] = D[colperm .+ end÷2,:]
    D[:,1:end÷2] = D[:,rowperm]

    #D,A
    D
end

struct RowSwap
    row1::Int
    row2::Int
end

function (op::RowSwap)(M::AbstractMatrix)
    res = copy(M)
    res[[op.row1,op.row2],:] = res[[op.row2,op.row1],:]
    res[:,[op.row1,op.row2]] = res[:,[op.row2,op.row1]]
    res
end

struct ColSwap
    col1::Int
    col2::Int
end

function (op::ColSwap)(M::AbstractMatrix)
    res = copy(M)
    res[:,[op.col1,op.col2]] = res[:,[op.col2,op.col1]]
    res[[op.col1,op.col2],:] = res[[op.col2,op.col1],:]
    res
end

"""
    Q′ = dual(Q)

Return the dual (i.e. the adjoint) matrix factorization of `Q`.

If `Q` factors the potential `f`, then this is a matrix factorization that
factors `-f`.
"""
function dual(M::AbstractMatrix)
    n, m = size(M)
    (n == m && iseven(n)) || throw(ArgumentError("Need square, even-rank matrix for applying MatrixFactorizations.dual"))

    [
       M[1:end÷2, 1:end÷2]     M[end÷2+1:end, 1:end÷2]    ;
      -M[1:end÷2, end÷2+1:end] M[end÷2+1:end, end÷2+1:end];
    ]
end

function inflatepowers(M::AbstractMatrix, var, exp)
    varval = eltype(M)(var)
    res = spzeros(eltype(M), (exp .* size(M))...)
    for row = axes(M, 1), col = axes(M, 2)
        curblock = @view res[(row-1)*exp+1:row*exp, (col-1)*exp+1:col*exp]
        for ((n,),c) in expand(M[row, col], var)
            for i = 0:exp-1
                j = mod(i+n, exp)
                curblock[j+1, i+1] += c*varval^((n + i - j)÷exp)
            end
        end
    end
    return res
end

function getpotential(A::AbstractMatrix)
    A_sq = A^2
    f = A_sq[1,1]
    A_sq == f*I || throw(ArgumentError("getpotential() needs a matrix factorization"))
    return f
end

function iscomposable(A::AbstractMatrix, B::AbstractMatrix, var_to_fuse, vars_to_fuse...)
    vars_to_fuse = [var_to_fuse; vars_to_fuse...]

    W = getpotential(A)
    V = getpotential(B)

    return all(iszero(diff(W + V, v)) for v in vars_to_fuse)
end

signedpermutations(A) = (((-1)^parity(σ), σ) for σ in permutations(eachindex(A)))


"""
    docstring goes here
"""
function fuse_abstract(A::AbstractMatrix, B::AbstractMatrix, var_to_fuse, vars_to_fuse...)
    vars_to_fuse = [var_to_fuse; vars_to_fuse...]
    Q = A⨶B

    iscomposable(A, B, vars_to_fuse...) || error("These matrix factorizations are not composable along $(join(vars_to_fuse, ","))")

    W = getpotential(A)
    V = getpotential(B)
    f = ofminring(V - constant_coefficient(V, vars_to_fuse...))

    R = base_extend(typeof(f))
    ∇f = [diff(f, v) for v in vars_to_fuse]
    ∇B = [diff(B, v) for v in vars_to_fuse]
    gr, tr = gröbner_transformation(∇f)

    var_data = map(vars_to_fuse) do v
        pow = 1
        while !iszero(rem(R(v)^pow, gr))
            pow += 1
            pow > 30 && error("Power for computing Jacobian of $R/$f is too high; exiting")
        end
        lift = div(R(v)^pow, gr)*tr
        λ = one(A)⨷sum(prod, zip(lift, ∇B))
        t = gensym()
        (v, pow, t, λ)
    end

    function inflate(M)
        for (v, pow, t, λ) in var_data
            M = inflatepowers(M, v, pow)
        end
        M
    end
    p(x) = constant_coefficient(x, vars_to_fuse...)

    QQ = inflate(Q)
    At = -sum(ε * prod(diff(QQ, v) for (v, pow, t, λ) in var_data[σ])
              for (ε, σ) in signedpermutations(var_data))//factorial(length(var_data))
    λλ = inflate(prod(λ for (v, pow, t, λ) in var_data))
    e = (-1)^length(var_data) * p(λλ * At)
    QQQ = p(QQ)

    @assert QQQ*e == e*QQQ "Internal error in LandauGinzburgCategories: computed non-morphism for e"
    #h = matrix_solve_affine(h->QQQ*h + h*QQQ, e^2 - e, size(QQQ))
    #@assert h != nothing "Internal error in LandauGinzburgCategories: computed a non-idempotent morphism for e"
    return QQQ, e
end

"""
    A_B = fuse(A, B, vars_to_fuse...)

Return a finite-rank representative of A⨶B by fusing the variables vars_to_fuse.
"""
function fuse(A::AbstractMatrix, B::AbstractMatrix, var_to_fuse, vars_to_fuse...)
    Q, e = fuse_abstract(A, B, var_to_fuse, vars_to_fuse...)

    sweepscalars!(Q, e, :u, :v, :x, :y)
    relevantrows = filter(1 : size(Q, 1)) do i
        count(!iszero, Q[i, :]) > 1 || count(!iszero, Q[:, i]) > 1
    end
    Q = Q[relevantrows, relevantrows]
    e = e[relevantrows, relevantrows]

    N = e^2 - e
    J = triangularperm(N)
    J′ = invperm(J)

    f = (I - exactinv(exactsqrt(UpperTriangular(I + 4N[J,J]))).data)[J′, J′]/2;
    E = e + f * (I - 2e)
    @assert iszero(E^2 - E)

    # J is reverse-sorted by degree. Experimentally, it's been useful
    # to do the Gröbner basis computation with the lowest degree monomials
    # as 'leading coefficients' for the module, so we use that order from
    # now on.
    E = E[reverse(J), reverse(J)]
    Q = Q[reverse(J), reverse(J)]

    im = mingenerators(columns(E))
    G = hcat(im...)
    d = size(G, 2)
    AB = matrix_solve_affine(AB -> G * AB, Q * G, (d, d))
    return AB
end

function ⊕(A::AbstractMatrix, B::AbstractMatrix)
    m1,n1 = size(A)
    m2,n2 = size(B)
    T = promote_type(eltype(A), eltype(B))
    res = [
        A                zeros(T, m1, n2);
        zeros(T, m2, n1) B;
    ]
end

function ⊞(A::AbstractMatrix, B::AbstractMatrix)
    C = to_alternating_grades(A)
    D = to_alternating_grades(B)
    return from_alternating_grades(C ⊕ D)
end

export ⨷, ⨶, ⊞, ⊕
export unit_matrix_factorization
export block_diagonalization
export rows, columns, topleft, topright, bottomleft, bottomright
export RowOp, ColOp, RowSwap, ColSwap
export dual


end
