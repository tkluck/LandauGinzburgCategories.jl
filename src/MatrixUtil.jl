module MatrixUtil

import Base: inv, transpose, adjoint

import LinearAlgebra: checksquare, UpperTriangular
import SparseArrays: dropzeros!, SparseMatrixCSC

import DataStructures: DefaultDict
import PolynomialRings: Polynomial, basering, deg, checkconstant
import ProgressMeter: Progress, update!, finish!

struct ColOp{T}
    target    :: Int
    source    :: Int
    coeffs    :: Matrix{T}
end

struct RowOp{T}
    target    :: Int
    source    :: Int
    coeffs    :: Matrix{T}
end

apply!(M::AbstractMatrix, op::RowOp) = (tgt = view(M, [op.target, op.source], :); tgt .= op.coeffs * tgt; M)
apply!(M::AbstractMatrix, op::ColOp) = (tgt = view(M, :, [op.target, op.source]); tgt .= tgt * transpose(op.coeffs); M)

apply(M::AbstractMatrix, op::Union{RowOp, ColOp}) = apply!(deepcopy(M), op)

ColOp(col, factor) = ColOp{typeof(factor)}(
    col, 1,
    [
        factor       zero(factor)
        zero(factor) one(factor)
    ],
)
ColOp(target, source_factor, source, target_factor=one(source_factor)) = ColOp{typeof(source_factor)}(
    target, source,
    [
        target_factor        source_factor
        zero(source_factor)  one(source_factor)
    ],
)

RowOp(row, factor) = RowOp{typeof(factor)}(
    row, 1,
    [
        factor       zero(factor)
        zero(factor) one(factor)
    ],
)
RowOp(target, source_factor, source, target_factor=one(source_factor)) = RowOp{typeof(source_factor)}(
    target, source,
    [
        target_factor        source_factor
        zero(source_factor)  one(source_factor)
    ],
)

function inv(op::Union{RowOp, ColOp})
    a, c, b, d = op.coeffs
    coeffs = [d -b; -c a] * inv(a*d - b*c)
    typeof(op)(op.target, op.source, coeffs)
end
transpose(op::RowOp) = ColOp{eltype(op.coeffs)}(op.target, op.source, transpose(op.coeffs))
transpose(op::ColOp) = RowOp{eltype(op.coeffs)}(op.target, op.source, transpose(op.coeffs))
adjoint(op::Union{RowOp, ColOp}) = transpose(inv(op))

function conjugate!(M::AbstractMatrix, op::Union{RowOp, ColOp})
    apply!(M, op)
    apply!(M, adjoint(op))
end

conjugate(M::AbstractMatrix, op::Union{RowOp, ColOp}) = conjugate!(deepcopy(M), op)

gettarget(M, op::ColOp, ix) = M[ix, op.target]
gettarget(M, op::RowOp, ix) = M[op.target, ix]
getsource(M, op::ColOp, ix) = M[ix, op.source]
getsource(M, op::RowOp, ix) = M[op.source, ix]

settarget!(M, op::ColOp, val, ix) = M[ix, op.target] = val
settarget!(M, op::RowOp, val, ix) = M[op.target, ix] = val
setsource!(M, op::ColOp, val, ix) = M[ix, op.source] = val
setsource!(M, op::RowOp, val, ix) = M[op.source, ix] = val

# avoid method ambiguity by separating RowOp/ColOp
apply!(M::SparseMatrixCSC, op::RowOp) = _sparseapply!(M, op)
apply!(M::SparseMatrixCSC, op::ColOp) = _sparseapply!(M, op)
function _sparseapply!(M, op)
    if iszero(op.coeffs[2, 1])
        for i in gettarget(M, op, :).nzind
            gettarget(M, op, i) .*= op.coeffs[1, 1]
        end
        for i in getsource(M, op, :).nzind
            val = gettarget(M, op, i)
            val .+= op.coeffs[1, 2] .* getsource(M, op, i)
            # in case val is zero, the sparse getindex returns
            # a newly constructed zero. We have to assign it back if we
            # operate in-place.
            settarget!(M, op, val, i)
        end
        for i in getsource(M, op, :).nzind
            getsource(M, op, i) .*= op.coeffs[2, 2]
        end
    elseif iszero(op.coeffs[1, 2])
        for i in getsource(M, op, :).nzind
            getsource(M, op, i) .*= op.coeffs[2, 2]
        end
        for i in gettarget(M, op, :).nzind
            val = getsource(M, op, i)
            val .+= op.coeffs[2, 1] .* gettarget(M, op, i)
            # See above
            setsource!(M, op, val, i)
        end
        for i in gettarget(M, op, :).nzind
            gettarget(M, op, i) .*= op.coeffs[1, 1]
        end
    else
        invoke(apply!, Tuple{AbstractMatrix, ColOp}, M, op)
    end
    M
end

function triangularperm(N::AbstractMatrix{<:Polynomial}, vars...)
    n = checksquare(N)

    maxdepth = DefaultDict(() -> 0)

    while true
        any_update = false
        ix = findfirst(!iszero, N)
        while ix != nothing
            i, j = ix[1], ix[2]
            if i != j
                degdiff = deg(N[ix], vars...)
                if (d = maxdepth[j] + degdiff) > maxdepth[i]
                    d > n && error("triangularperm: argument is not triangular")
                    maxdepth[i] = d
                    any_update = true
                end
            end

            ix = findnext(!iszero, N, nextind(N, ix))
        end
        !any_update && break
    end
    I = sort!([1:n;], by = i -> -maxdepth[i])
    return I
end

function sweepscalars!(M::AbstractMatrix{<:Polynomial}, e, vars...)
    p = Progress(size(M, 2), 1, "Sweeping rows/columns containing a scalar")
    ix = findfirst(!iszero, M)
    #Msq = M^2
    while ix != nothing
        i, j = ix[1], ix[2]
        update!(p, j)
        #if deg(M[ix], vars...) == 0
        if length(M[ix].monomials) == 1 && isone(M[ix].monomials[1])
            op = RowOp(i, inv(M[ix]))
            conjugate!(M, op); conjugate!(e, op)
            for k in M[:, j].nzind
                k == i && continue
                op = RowOp(k, -M[k, j], i)
                conjugate!(M, op); conjugate!(e, op)
            end
            for k in M[i, :].nzind
                k == j && continue
                op = ColOp(k, -M[i, k], j)
                conjugate!(M, op); conjugate!(e, op)
            end
            #@assert sum(M[:, j]) == M[ix] == sum(M[i, :])
            #@assert iszero(Msq - M^2)
            #@assert M*e == e*M
            dropzeros!(M)
        end
        ix = findnext(!iszero, M, nextind(M, ix))
    end
    finish!(p)
    return M
end

function exactsqrt(a::Integer)
    m = round(typeof(a), sqrt(a))
    m^2 == a && return m
    (m - 1)^2 == a && return m - 1
    error("Cannot compute square root of $a")
end

exactsqrt(a::Number) = isone(a) ? a : error("Exact sqrt for $(typeof(a)) not implemented")
exactsqrt(a::Rational) = exactsqrt(numerator(a)) // exactsqrt(denominator(a))
exactsqrt(f::Polynomial) = exactsqrt(checkconstant(f))

function exactsqrt(A::UpperTriangular)
    B = A.data
    n = checksquare(B)
    R = zero(B)
    for j in 1:n
        nz_rows = Int[]
        R[j, j] = exactsqrt(B[j, j])
        for i in reverse(1 : j - 1)
            r = B[i, j]
            for k in nz_rows
                r -= R[i, k] * R[k, j]
            end
            if !iszero(r)
                R[i, j] = r / (checkconstant(R[i, i]) + checkconstant(R[j, j]))
                push!(nz_rows, i)
            end
        end
    end
    return UpperTriangular(R)
end

function exactinv(B::UpperTriangular)
    A = one(B.data)
    n = checksquare(B)
    for i = 1:n
        nz_cols = Int[]
        for j = 1:n
            Aij = A[i, j]
            for k in nz_cols
                Aij -= A[i, k] * B.data[k, j]
            end
            if !iszero(Aij)
                A[i, j] = Aij / checkconstant(B.data[j, j])
                push!(nz_cols, j)
            end
        end
    end
    return UpperTriangular(A)
end

rows(M::AbstractMatrix) =    [transpose(M[i,:]) for i=axes(M,1)]
columns(M::AbstractMatrix) = [M[:,i]  for i=axes(M,2)]
topleft(M::AbstractMatrix)     = M[1:end÷2,     1:end÷2]
topright(M::AbstractMatrix)    = M[1:end÷2,     end÷2+1:end]
bottomleft(M::AbstractMatrix)  = M[end÷2+1:end, 1:end÷2]
bottomright(M::AbstractMatrix) = M[end÷2+1:end, end÷2+1:end]



end
