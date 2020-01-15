using Test
using PolynomialRings
using LandauGinzburgCategories
using LandauGinzburgCategories.Library

import LinearAlgebra: I

@testset "Landau-Ginzburg Categories" begin
    @testset "Tensor product fusion" begin
        R = @ring! ℚ[x,y,z]

        A = unit_matrix_factorization(x^3, x = y)
        B = unit_matrix_factorization(y^3, y = z)

        @test A^2 == (y^3 - x^3) * I
        @test B^2 == (z^3 - y^3) * I
        @test quantum_dimensions(A, [:x], [:y]) == (1, -1)

        @test (A ⨶ B)^2 == (z^3 - x^3) * I

        Q, e = fuse_abstract(A, B, :y)
        @test Q^2 == (z^3 - x^3) * I
        @test Q*e == e*Q
        @test fuse(A, B, :y)^2 == (z^3 - x^3) * I
    end

    @testset "Library" begin
        R = @ring! ℚ[x,y,u,v]

        Q = orbifold_equivalence(TwoVariables.A{9}, TwoVariables.D{6}, [x, y], [u, v])
        @test Q^2 == (TwoVariables.D{6}(u, v) - TwoVariables.A{9}(x, y))*I

        @test false == orbifold_equivalence(TwoVariables.A{6}, TwoVariables.D{6})

        for (name, closure) in LandauGinzburgCategories.Library.WellKnownEquivalences
            @testset "$name" begin
                method = methods(closure).ms[1]
                if method.nargs == 1
                    @test orbifold_equivalence(closure()...) isa AbstractMatrix
                else
                    for n = 3:10
                        @test orbifold_equivalence(closure(n)...) isa AbstractMatrix
                    end
                end
            end
        end

    end
end
