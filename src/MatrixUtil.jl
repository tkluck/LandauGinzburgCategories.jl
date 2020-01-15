module MatrixUtil

import LinearAlgebra: checksquare, UpperTriangular

import DataStructures: DefaultDict
import PolynomialRings: Polynomial, basering

function triangularperm(N::AbstractMatrix)
    n = checksquare(N)

    maxdepth = DefaultDict(() -> 0)

    while true
        any_update = false
        ix = findfirst(!iszero, N)
        while ix != nothing
            i, j = ix[1], ix[2]
            if i != j
                if (d = maxdepth[j] + 1) > maxdepth[i]
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

checkconstant(a::Number) = a
checkconstant(f::Polynomial) = convert(basering(f), f)

function exactsqrt(a::Integer)
    m = round(typeof(a), sqrt(a))
    m^2 == a && return m
    (m - 1)^2 == a && return m - 1
    error("Cannot compute square root of $a")
end

exactsqrt(a::Rational) = exactsqrt(numerator(a)) // exactsqrt(denominator(a))
exactsqrt(f::Polynomial) = exactsqrt(checkconstant(f))

function exactsqrt(A::UpperTriangular)
    B = A.data
    n = checksquare(B)
    R = zero(B)
    for j in 1:n
        R[j, j] = exactsqrt(B[j, j])
        for i in reverse(1 : j - 1)
            r = B[i, j]
            for k in i + 1 : j - 1
                r -= R[i, k] * R[k, j]
            end
            iszero(r) || (R[i, j] = r / (checkconstant(R[i, i]) + checkconstant(R[j, j])))
        end
    end
    return UpperTriangular(R)
end

function exactinv(B::UpperTriangular)
    A = one(B.data)
    n = checksquare(B)
    for i = 1:n
        for j = 1:n
            Aij = A[i, j]
            for k = 1 : j - 1
                Aij -= A[i, k] * B.data[k, j]
            end
            A[i, j] = Aij / checkconstant(B[j, j])
        end
    end
    return UpperTriangular(A)
end


end
