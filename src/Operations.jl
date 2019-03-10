module Operations

import Combinatorics: permutations, parity
import LinearAlgebra: I
import SparseArrays: sparse

import PolynomialRings: Polynomial, polynomial_ring
import PolynomialRings: constant_coefficient, gröbner_transformation
import PolynomialRings: map_coefficients
import PolynomialRings.Expansions: _expansion_types, expansion
import PolynomialRings.QuotientRings: QuotientRing
import PolynomialRings.MonomialOrderings: MonomialOrder
import PolynomialRings.NamingSchemes: Named
import PolynomialRings.Solve: matrix_solve_affine

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
    n%2 == m%2 == k%2 == l%2 == 0 || throw(ValueError("Need ℤ/2 graded matrices, even rank and corank"))

    A,B = to_alternating_grades.((A,B))

    res = [a*B for a in A]

    for row in axes(res,1), col in axes(res,2)
        inner = res[row,col]
        for inner_row in axes(inner,1), inner_col in axes(inner,2)
            # Koszul signs; see formula in the documentation string
            if col%2 == 0 && (inner_row+inner_col)%2 == 1
                inner[inner_row, inner_col] *= -1
            end
        end
        # make sure gradings alternate in the output matrix: they alternate
        # in both inner and outer, and the result grading is their sum. At a
        # boundary between two 'inner' matrices, both signs change, so the
        # adjacent rows/columns have the same sign. We can fix this by reversing
        # the rows/columns within an inner matrix, depending on the sign of the
        # outer row/column.
        if row%2 == 0
            inner[:,:] = inner[end:-1:1,:]
        end
        if col%2 == 0
            inner[:,:] = inner[:,end:-1:1]
        end
    end
    return from_alternating_grades(flatten_blocks(res))
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

struct RowOp
    target_row::Int
    target_factor
    source_row::Int
    source_factor
end


function (op::RowOp)(M::AbstractMatrix)
    n, m = size(M)
    (n == m && n%2 == 0) || throw(ArgumentError("Need square, even-rank matrix for applying MatrixFactorizations.RowOp"))
    res = copy(M)
    if op.source_row == 0
        res[op.target_row,:] *= op.target_factor
        res[:,op.target_row] //= op.target_factor
    elseif (op.target_row <= n÷2 && op.source_row <= n÷2) || (op.target_row > n÷2 && op.source_row > n÷2)
        # row operation
        target_topright = @view res[[op.target_row,op.source_row],:]
        t = [op.target_factor op.source_factor; 0 1]
        target_topright .= t * target_topright
        # corresonding column operation in the other block
        target_bottomleft = @view res[:,[op.target_row,op.source_row]]
        t_inv = [1 -op.source_factor; 0 op.target_factor]//op.target_factor
        target_bottomleft .= target_bottomleft * t_inv
    else
        throw(ArgumentError("MatrixFactorizations.RowOp needs to have both rows in the same graded block"))
    end
    return res
end

struct ColOp
    target_col::Int
    target_factor
    source_col::Int
    source_factor
end

function (op::ColOp)(M::AbstractMatrix)
    n, m = size(M)
    (n == m && n%2 == 0) || throw(ArgumentError("Need square, even-rank matrix for applying MatrixFactorizations.ColOp"))
    res = copy(M)
    if op.source_col == 0
        res[:,op.target_col] *= op.target_factor
        res[op.target_col,:] //= op.target_factor
    elseif (op.target_col <= n÷2 && op.source_col <= n÷2) || (op.target_col > n÷2 && op.source_col > n÷2)
        # column operation
        target_topright = @view res[:,[op.target_col,op.source_col]]
        t = [op.target_factor 0; op.source_factor 1]
        target_topright .= target_topright * t
        # corresonding column operation in the other block
        target_bottomleft = @view res[[op.target_col,op.source_col],:]
        t_inv = [1 0; -op.source_factor op.target_factor]//op.target_factor
        target_bottomleft .= t_inv * target_bottomleft
    else
        throw(ArgumentError("MatrixFactorizations.ColOp needs to have both rows in the same graded block"))
    end
    return res
end

RowOp(target_row, factor) = RowOp(target_row, factor, 0, 0)
ColOp(target_col, factor) = ColOp(target_col, factor, 0, 0)
RowOp(target_row, source_factor, source_row) = RowOp(target_row, 1, source_row, source_factor)
ColOp(target_col, source_factor, source_col) = ColOp(target_col, 1, source_col, source_factor)

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

function dual(M::AbstractMatrix)
    n, m = size(M)
    (n == m && n%2 == 0) || throw(ArgumentError("Need square, even-rank matrix for applying MatrixFactorizations.dual"))

    [
       M[1:end÷2, 1:end÷2]     M[end÷2+1:end, 1:end÷2]    ;
      -M[1:end÷2, end÷2+1:end] M[end÷2+1:end, end÷2+1:end];
    ]
end

function matrix_over_subring(M::AbstractMatrix, var, exp, substitution_var)
    _, P = _expansion_types(eltype(M), MonomialOrder{:degrevlex, Named{(var,)}}())
    R, (substitution_var_val,) = polynomial_ring(substitution_var, basering=Int)
    S = promote_type(P, R)
    blocks = map(M) do f
        res = zeros(S, exp, exp)
        for ((n,),c) in expansion(f, var)
            for i = 0:exp-1
                j = mod(i+n, exp)
                res[j+1,i+1] += c*substitution_var_val^((n + i - j)÷exp)
            end
        end
        res
    end
    return flatten_blocks(blocks)
end

function getpotential(A::AbstractMatrix)
    A_sq = A^2
    f = A_sq[1,1]
    A_sq == f*I || throw(ArgumentError("getpotential() needs a matrix factorization"))
    return f
end

small_fractions = [0//1;[a//b for a=-2000:2000 for b = 1:2000 if gcd(a,b) == 1]]
sort!(small_fractions)
function round_to_small_fraction(a::Rational)
    ix = searchsortedfirst(small_fractions, a)
    (ix <= 1 || ix >= length(small_fractions)) && throw(ErrorException())
    below, above = small_fractions[ix-1:ix]
    δ₂ = above - a
    δ₁ = a - below
    return δ₂ < δ₁ ? above : below
end
function find_zero_in_small_fractions(f, df, x0)
    max_loops = 7
    x = x0
    for _ in 1:max_loops
        f_x = f(x)
        iszero(f_x) && return x
        # Newton's method
        δx = matrix_solve_affine(δx->df(x,δx), f_x, size(x))
        x = x - δx

        # after a sufficient number of iterations, we should be close enough
        # that rounding to a small fraction gives a zero.
        x_rounded = map_coefficients.(round_to_small_fraction, x)
        iszero(f(x_rounded)) && return x_rounded
    end
    return nothing
end

function find_analytic_equality(f, y, df, x0)
    max_loops = 7
    x = x0
    for _ in 1:max_loops
        f_x = f(x)
        f_x ≈ y && return x
        # Newton's method
        δx = matrix_solve_affine(δx->df(x,δx), f_x - y, size(x))
        x = x - δx
    end
    return nothing
end

signedpermutations(A) = (((-1)^parity(σ), σ) for σ in permutations(eachindex(A)))

function fuse_abstract(A::AbstractMatrix, B::AbstractMatrix, var_to_fuse, vars_to_fuse...)
    vars_to_fuse = [var_to_fuse; vars_to_fuse...]
    R = promote_type(eltype(A), eltype(B))
    Q = A⨶B

    W, V = getpotential.((A,B))
    f = V - constant_coefficient(V, vars_to_fuse...)
    ∇f = [diff(f, v) for v in vars_to_fuse]
    ∇B = [diff(B, v) for v in vars_to_fuse]
    gr, tr = gröbner_transformation(∇f)

    var_data = map(vars_to_fuse) do v
        pow = 1
        while !iszero(rem(R(v)^pow, gr))
            pow += 1
            pow > 30 && error("Power for computing Jacobian is too high; exiting")
        end
        lift = div(R(v)^pow, gr)*tr
        λ = one(A)⨷sum(prod, zip(lift, ∇B))
        t = gensym()
        (v, pow, t, λ)
    end

    function inflate(M)
        for (v, pow, t, λ) in var_data
            M = matrix_over_subring(M, v, pow, t)
        end
        M
    end
    function p(x)
        for (v, pow, t, λ) in var_data
            x = x(; t => 0)
        end
        x
    end

    QQ = inflate(Q)
    At = -sum(ε * prod(diff(QQ, t) for (v, pow, t, λ) in var_data[σ])
              for (ε, σ) in signedpermutations(var_data))//factorial(length(var_data))
    λλ = inflate(prod(λ for (v, pow, t, λ) in var_data))
    e = (-1)^length(var_data) * p(λλ * At)
    QQQ = p(QQ)

    return QQQ, e
end

function fuse(A::AbstractMatrix, B::AbstractMatrix, var_to_fuse, vars_to_fuse...)
    R = promote_type(eltype(A), eltype(B))
    QQQ, e = fuse_abstract(A, B, var_to_fuse, vars_to_fuse)
    QQQ = sparse(QQQ)
    e = sparse(e)

    h = matrix_solve_affine(h->QQQ*h + h*QQQ, e^2 - e, size(QQQ))
    @assert h != nothing "e should be idempotent up to homotopy"

    f(b) = -b + h + e*b + b*e + (b^2*QQQ + QQQ*b^2)//2 + b*QQQ*b
    df(b, xi) = -xi + e*xi + xi*e + ((b*xi + xi*b)*QQQ + QQQ*(b*xi + xi*b))//2 + xi*QQQ*b + b*QQQ*xi
    b = find_analytic_zero(f, df, zero(h))
    b === nothing && error("Failed to strictify e")
    e_strict = e + QQQ*b + b*QQQ
    @assert e_strict^2 == e_strict

    # naive splitting of the idempotent
    image = hcat(gröbner_basis(columns(e_strict))...)
    d = size(image, 2)
    return todense( matrix_solve_affine(Q -> image*Q, QQQ*image, (d,d)) )
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

rows(M::AbstractMatrix) =    [transpose(M[i,:]) for i=axes(M,1)]
columns(M::AbstractMatrix) = [M[:,i]  for i=axes(M,2)]
topleft(M::AbstractMatrix)     = M[1:end÷2,     1:end÷2]
topright(M::AbstractMatrix)    = M[1:end÷2,     end÷2+1:end]
bottomleft(M::AbstractMatrix)  = M[end÷2+1:end, 1:end÷2]
bottomright(M::AbstractMatrix) = M[end÷2+1:end, end÷2+1:end]

export ⨷, ⨶, ⊞, ⊕
export unit_matrix_factorization
export block_diagonalization
export rows, columns, topleft, topright, bottomleft, bottomright
export RowOp, ColOp, RowSwap, ColSwap
export dual
export matrix_over_subring



end
