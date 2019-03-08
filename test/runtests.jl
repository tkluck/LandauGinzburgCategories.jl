using Test
using PolynomialRings
using LandauGinzburgCategories

import LinearAlgebra: I

@testset "Tensor product fusion" begin
    R = @ring! ℚ[x,y,z]

    A = unit_matrix_factorization(x^3, [:x], [:y])
    B = unit_matrix_factorization(y^3, [:y], [:z])

    @test A^2 == (x^3 - y^3) * I
    @test B^2 == (y^3 - z^3) * I

    @test (A ⨶ B)^2 == (x^3 - z^3) * I

    @test fuse(A ⨶ B, :y)^2 == (x^3 - z^3) * I

end
