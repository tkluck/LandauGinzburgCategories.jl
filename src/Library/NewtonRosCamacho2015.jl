"""
Orbifold equivalences from

> Rachel Newton, Ana Ros Camacho
> "Strangely dual orbifold equivalence I"
> arXiv preprint https://arxiv.org/abs/1509.08069v1
"""

function ARC_shaped_equivalence(;substitutions...)
    @ring! ℚ[d15,d16,d17,d25,d26,d35]

    [   0    0    0    0  d15  d16  d17    0;
        0    0    0    0  d25  d26    0  d17;
        0    0    0    0  d35    0  d26 -d16;
        0    0    0    0    0  d35 -d25  d15;
      d26 -d16 -d17    0    0    0    0    0;
     -d25  d15    0 -d17    0    0    0    0;
     -d35    0  d15  d16    0    0    0    0;
        0 -d35  d25  d26    0    0    0    0](;substitutions...)
end

function _orbifold_equivalence_def(::Type{ThreeVariables.E₁₄{:v1}}, ::Type{ThreeVariables.E₁₄{:v2}}, left_vars, right_vars)
    F = @ring! ℚ[c1]/(c1^4 - 2c1^2 + 2) # other component is given by c1^4 + 2c1^2 + 2)
    @ring! F[a2,a3,a4,a5,b22,b23,b24][x,y,z,u,v,w]

    # ARC: conditions on the parameters -----------------------
    b21=-1//2*(-2*b22*c1+2*b24*c1^2-2*b23*c1^3+c1^4+2*a5*c1^4)
    a1=b21+a2*c1-b22*c1-a4*c1^2+b24*c1^2+a3*c1^3-b23*c1^3+c1^4
    d11=a1*a5+a5*b21+a3*b22+a2*b23+a4*b24-a2*a5*c1-a5*b22*c1-a4*b23*c1-a3*b24*c1+a4*a5*c1^2+a3*b23*c1^2+a5*b24*c1^2-a3*a5*c1^3-a5*b23*c1^3+a5^2*c1^4
    b25=a5
    d14=a5*b25
    d155=a5*b23+a3*b25-c1*d14
    d13=a3*b23+a5*b24+a4*b25-c1*d155
    d12=a5*b22+a4*b23+a3*b24+a2*b25-c1*d13
    d10=a5+b25
    d8=a3*b21+a4*b22+a1*b23+a2*b24-c1*d11
    d7=a3+b23-c1*d10
    d5=a4*b21+a2*b22+a1*b24-a3*b21*c1-a4*b22*c1-a1*b23*c1-a2*b24*c1+a1*a5*c1^2+a5*b21*c1^2+a3*b22*c1^2+a2*b23*c1^2+a4*b24*c1^2-a2*a5*c1^3-a5*b22*c1^3-a4*b23*c1^3-a3*b24*c1^3+a4*a5*c1^4+a3*b23*c1^4+a5*b24*c1^4-a3*a5*c1^5-a5*b23*c1^5+a5^2*c1^6
    d2=a2*b21+a1*b22-c1*d5
    d4=a4+b24-a3*c1-b23*c1+2*a5*c1^2
    d1=a2+b22-c1*d4
    d9=-a3+b23-c1
    d6=-a4+b24-c1*d9
    d3=-a2+b22-c1*d6
    # ----------------------------------------------------------

    # ARC: In this case, we used directly the structure of the ansatz for a matrix factorization, which makes us have less parameters to take care of
    Q = ARC_shaped_equivalence(
        d15=z+w+a1*u^4+a2*x*u^3+a3*x^3*u+a4*x^2*u^2+a5*x^4,
        d16=y^2+y*v+v^2,
        d17=x^3*z+d1*u^3*w+d2*u^7+d3*z*u^3+d4*x*u^2*w+d5*x*u^6+d6*x*z*u^2+d7*x^2*u*w+d8*x^2*u^5+d9*x^2*z*u+d10*x^3*w+d11*x^3*u^4+d12*x^4*u^3+d13*x^5*u^2+d14*x^7+d155*x^6*u,
        d25=y-v,
        d26=-z+w+b21*u^4+b22*x*u^3+b23*x^3*u+b24*x^2*u^2+b25*x^4,
        d35=x+c1*u,
    )
    # ----------------------------------------------------------

    # express in the user's variables
    Q = Q(x=left_vars[1],  y=left_vars[2],  z=left_vars[3],
          u=right_vars[1], v=right_vars[2], w=right_vars[3])

    return Q
end

function _orbifold_equivalence_def(::Type{ThreeVariables.U₁₂{:v3}}, ::Type{ThreeVariables.U₁₂{:v1}}, left_vars, right_vars)
    F = @ring! ℚ[a1, a2, b1, b2]/(
        -1 + a1^2*b1 - a1*b1^2,
        2a1*b1*a2 - b1^2*a2 + a1^2*b2 - 2a1*b1*b2,
        b1*a2^2 + 2a1*a2*b2 - 2b1*a2*b2 - a1*b2^2,
        -1 + a2^2*b2 - a2*b2^2,
    )
    @ring! F[c1,c2,x,y,z,u,v,w]

    Q = ARC_shaped_equivalence(
        d15 = z + a1*v + a2*w,
        d16 = -y*z + (-a1^2 + a1*b1 + a1*c1)*v^2 + (-2a1*a2 + b1*a2 + c1*a2 + a1*b2 + a1*c2)*v*w +
              + a2*(-a2+b2+c2)*w^2 + c1*v*z + c2*w*z,
        d17 = u^3 + u^2*x + u*x^2 + x^3,
        d25 = -y + b1*v + b2*w,
        d26 = -y*z + b1*c1*v^2 + (c1*b2 + b1*c2)*v*w + b2*c2*w^2 - (-a1+b1+c1)*v*y -
              (-a2+b2+c2)*w*y,
        d35 = x - u,
    )
    # express in the user's variables
    Q = Q(x=left_vars[1],  y=left_vars[2],  z=left_vars[3],
          u=right_vars[1], v=right_vars[2], w=right_vars[3])

    return Q
end

function _orbifold_equivalence_def(::Type{ThreeVariables.W₁₂{:v2}}, ::Type{ThreeVariables.W₁₂{:v1}}, left_vars, right_vars)
    F = @ring! ℚ[a,b,c]/(
        a^2-a*b,
        -a*b^3+4*a^3*c+4*a*b^2*c-5*a*b*c^2+2*a*c^3+a^2*(b^2-6*b*c+5*c^2)+1//4*(4+b^4-4*b^3*c+6*b^2*c^2-4*b*c^3+c^4),
    )
    @ring! F[x,y,z,u,v,w]

    Q1 = [w+a*u*x+1//2*(-2*a*b+b^2+2*a*c-2*b*c+c^2)*x^2+z u*w+c*w*x+(1//2*(-a+b)*(-2*a*b+b^2+2*a*c-2*b*c+c^2)+a*(b*(-2*a+b-c)+1//2*(2*a*b-b^2-2*a*c+2*b*c-c^2)))*x^3+b*x*z v^4+v^3*y+v^2*y^2+v*y^3+y^4 0
          u+(-2*a+b-c)*x -w+(-a+b)*u*x+(b*(-2*a+b-c)+1//2*(2*a*b-b^2-2*a*c+2*b*c-c^2))*x^2+z 0 v^4+v^3*y+v^2*y^2+v*y^3+y^4
          v-y 0 -w+(-a+b)*u*x+(b*(-2*a+b-c)+1//2*(2*a*b-b^2-2*a*c+2*b*c-c^2))*x^2+z -u*w-c*w*x+(-(1//2)*(-a+b)*(-2*a*b+b^2+2*a*c-2*b*c+c^2)-a*(b*(-2*a+b-c)+1//2*(2*a*b-b^2-2*a*c+2*b*c-c^2)))*x^3-b*x*z
          0 v-y -u+(2*a-b+c)*x w+a*u*x+1//2*(-2*a*b+b^2+2*a*c-2*b*c+c^2)*x^2+z
         ]

    Q2 = [w+(a-b)*u*x+(-b*(-2*a+b-c)+1//2*(-2*a*b+b^2+2*a*c-2*b*c+c^2))*x^2-z u*w+c*w*x+(1//2*(-a+b)*(-2*a*b+b^2+2*a*c-2*b*c+c^2)+a*(b*(-2*a+b-c)+1//2*(2*a*b-b^2-2*a*c+2*b*c-c^2)))*x^3+b*x*z v^4+v^3*y+v^2*y^2+v*y^3+y^4 0
          u+(-2*a+b-c)*x -w-a*u*x+1//2*(2*a*b-b^2-2*a*c+2*b*c-c^2)*x^2-z 0 v^4+v^3*y+v^2*y^2+v*y^3+y^4
          v-y 0 -w-a*u*x+1//2*(2*a*b-b^2-2*a*c+2*b*c-c^2)*x^2-z -u*w-c*w*x+(-(1//2)*(-a+b)*(-2*a*b+b^2+2*a*c-2*b*c+c^2)-a*(b*(-2*a+b-c)+1//2*(2*a*b-b^2-2*a*c+2*b*c-c^2)))*x^3-b*x*z
          0 v-y -u+(2*a-b+c)*x w+(a-b)*u*x+(-b*(-2*a+b-c)+1//2*(-2*a*b+b^2+2*a*c-2*b*c+c^2))*x^2-z
         ]

    z = zero(Q1)

    Q = [z Q1; Q2 z]

    # express in the user's variables (note different convention for (x,y) in
    # the shapes of E₁₄{:v1} and E₁₄{:v2})
    Q = Q(y=left_vars[1],  x=left_vars[2],  z=left_vars[3],
          v=right_vars[1], u=right_vars[2], w=right_vars[3])

    return Q
end
