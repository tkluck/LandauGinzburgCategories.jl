"""
Orbifold equivalences from

> ...
"""

import PolynomialRings: Ideal, map_coefficients, base_extend

function _orbifold_equivalence_def(::Type{TwoVariables.A₂A₂}, ::Type{TwoVariables.A₅}, left_vars, right_vars)
    C = Complex{BigInt}
    R = @ring! C[a[]][x,y,u,v]

    # compatibility with copy-pasted notation below
    A = a
    a(i) = A[i]
    lookup = Dict(1:4 .=> (u,v,x,y))
    x(i) = lookup[i]
    # from https://nms.kcl.ac.uk/andreas.recknagel/oeq-page/defectslistforweb.txt
    Q = zeros(R, 4, 4)
    Q[1,3]=x(1)^2-x(3)*a(2)-x(4)*a(2)
    Q[1,4]=x(1)*x(3)*a(1)-x(1)*x(4)*a(1)+x(2)
    Q[2,3]=-x(1)*x(3)*a(1)+x(1)*x(4)*a(1)+x(2)
    Q[2,4]=-64*x(4)^2*a(2)^8+16*x(3)*x(4)*a(2)^5-x(1)^4-x(1)^2*x(3)*a(2)-x(1)^2*x(4)*a(2)-4*x(3)^2*a(2)^2
    Q[3,1]=64*x(4)^2*a(2)^8-16*x(3)*x(4)*a(2)^5+x(1)^4+x(1)^2*x(3)*a(2)+x(1)^2*x(4)*a(2)+4*x(3)^2*a(2)^2
    Q[3,2]=x(1)*x(3)*a(1)-x(1)*x(4)*a(1)+x(2)
    Q[4,1]=-x(1)*x(3)*a(1)+x(1)*x(4)*a(1)+x(2)
    Q[4,2]=-x(1)^2+x(3)*a(2)+x(4)*a(2)

    # post-processing: sign difference in our convention for A₅
    Q = Q(v=im*v)

    # express in the user's variables
    Q = Q(x=left_vars[1], y=left_vars[2], u=right_vars[1], v=right_vars[2])

    # quotient out the equations resulting from imposing Q^2 == f - g
    S′ = @ring! C[a[]]
    eqns = [
        -64a[2]^9 + 1
        -64a[2]^9 + 16a[2]^6
         16a[2]^6 -  4a[2]^3
         -4a[2]^3 + 1
         64a[2]^8 -  a[1]^2 -  a[2]^2
        -16a[2]^5 + 2a[1]^2 - 2a[2]^2
          -a[1]^2 + 3a[2]^2
    ]

    S = S′/Ideal(eqns)
    Q = base_extend.(Q, Ref(S))

    return Q
end
