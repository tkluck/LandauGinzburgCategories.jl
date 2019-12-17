"""
Orbifold equivalences from

> Rachel Newton, Ana Ros Camacho
> "Strangely dual orbifold equivalence I"
> arXiv preprint https://arxiv.org/abs/1509.08069v1
"""

function _orbifold_equivalence_def(::Type{ThreeVariables.Q₁₀}, ::Type{ThreeVariables.E₁₄{:v1}}, left_vars, right_vars)
    R = @ring! ℚ[a,b,c]/(
        # f1 ; other component is f2 = 4+3a^4+8a^3*b+8a^2*b^2-4a^3*c-8a^2*b*c
        -4+3a^4+8a^3*b+8a^2*b^2-4a^3*c-8*a^2*b*c,
        # g
        a^2*(a^4-8a^2*b^2-16a*b^3-8b^4+8a^2*b*c+24a*b^2*c+16b^3*c-2a^2*c^2-8a*b*c^2-8b^2*c^2),
    )
    S = @ring! R[x,y,z,u,v,w]

    κ1 = a^3//2 + a^2*b + a*b^2 - a^2*c//2 - a*b*c
    κ2 = 1 + 3a^4//4 + 3a^3*b + 4a^2*b^2 + 2a*b^3  - a^3*c - 3*a^2*b*c - 2a*b^2*c

    Q = ARC_shaped_equivalence(
        d15 = κ1*u^3+a*u*x+z,
        d16 = v^2+v*y+y^2,
        d17 = κ2*u^4//2+w-a*(-a-2b)*u^2*x//2+x^2+b*u*z,
        d25 = y-v,
        d26 = (-b-b^2*κ1+(c-a)*κ2//2)*u^5+(-a-2b+c)*u*w+c*u*x^2+b*(-a-b+c)*u^2*z-x*z,
        d35 = (-1+(-a-2b+c)*κ1+κ2//2)*u^4-w+a*(-a-2b+2c)*u^2*x//2+x^2+(-a-b+c)*u*z,
    )

    # express in the user's variables
    Q = Q(x=left_vars[1],  y=left_vars[2],  z=left_vars[3],
          u=right_vars[1], v=right_vars[2], w=right_vars[3])

    return Q
end

function _orbifold_equivalence_def(::Type{ThreeVariables.Q₁₀}, ::Type{ThreeVariables.E₁₄{:v2}}, left_vars, right_vars)
    # compared to the desciption in the paper, we add an extra parameter ("e")
    # and an extra equation be == 1. This allows representing division by b
    # as a polynomial.
    R = @ring! ℚ[a,b,c,d,e]/(
        a^2 - 1,
        b^2 + 4*c*e - c^2 - 4c*d + b*c^2*d + 2b*c*d^2,
        -2 + 2b*c + 2c^2*e^2 - c^4//4 + 2b*d - 2c^2*d*e + c^2*d^2,
        -2e^2 + 2d*e - d^2,
        b*e - 1,
    )
    S = @ring! R[x,y,z,u,v,w]

    Q = ARC_shaped_equivalence(
        d15 = b*v^3 + c*v*x + z,
        d16 = u^2 + u*y + y^2,
        d17 = v^4 + a*w + (1//2)*(c^2 + 2c*d)*v^2*x + x^2 + d*v*z,
        d25 =  -u + y,
        d26 = -2a*v*w*e + (b + 2c*e^2 - 2c*d*e + c^2*d + 2c*d^2 - (c^2 + 2c*d)*e)*v^3*x +
              (-2e + c+ 2*d)*v*x^2 - 2v^2*z*e^2 - x*z,
        d35 = -v^4 - a*w + (c*(-2e + c + 2d) + (1//2)*(-c^2 - 2c*d))*v^2*x +
              x^2 + (-2e + d)*v*z,
    )

    # express in the user's variables
    Q = Q(x=left_vars[1],  y=left_vars[2],  z=left_vars[3],
          u=right_vars[2], v=right_vars[1], w=right_vars[3])

    return Q
end
