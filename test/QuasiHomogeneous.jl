using Test

import PolynomialRings: namingscheme, @namingscheme, @ring!
import LandauGinzburgCategories.QuasiHomogeneous: Gradings, quasidegree, find_quasihomogeneous_degrees
import LandauGinzburgCategories.QuasiHomogeneous: centralcharge, generic_quasihomogeneous_polynomial
import LandauGinzburgCategories.QuasiHomogeneous: generic_quasihomogeneous_array

@testset "QuasiHomogeneous" begin
    @testset "Gradings" begin
        @ring! Int[x,y]
        gr = Gradings(x=3, y=4)
        @test namingscheme(gr) == @namingscheme((x,y))
        @test quasidegree(x^2 * y^4, gr) == 22

        @test find_quasihomogeneous_degrees(x^4 + y^2, x, y) == Gradings(x=1, y=2)
        @test find_quasihomogeneous_degrees(x^4 + x*y^2, x, y) == Gradings(x=2, y=3)

        @test centralcharge(x^4 + y^2, x, y) == 1//2
        @test centralcharge(x^4 + x*y^2, x, y) == 3//4
    end

    @testset "Generic polynomials" begin
        @ring! Int[c[]][x,y]

        @test generic_quasihomogeneous_polynomial(4, Gradings(x=1, y=2), reset(c)) ==
                c[3]*x^4 + c[2]*x^2*y + c[1]*y^2
        @test generic_quasihomogeneous_polynomial(3, Gradings(x=1, y=1), reset(c)) ==
                c[4]*x^3 + c[3]*x^2*y + c[2]*x*y^2 + c[1]*y^3

        @test generic_quasihomogeneous_array([1 2; 3 4], Gradings(x=1, y=2), reset(c)) == [
            c[1]*x              c[5]*x^2 + c[4]*y
            c[3]*x^3 + c[2]*x*y c[8]*x^4 + c[7]*x^2*y + c[6]*y^2
        ]
    end
end
