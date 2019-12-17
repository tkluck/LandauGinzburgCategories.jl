"""
Orbifold equivalences from (for small n)

> Orbifold completion of defect bicategories
> Nils Carqueville, Ingo Runkel
> http://arxiv.org/abs/1210.6363v4

And for general n from

> Orbifold equivalent potentials
> Nils Carqueville, Ana Ros Camacho, Ingo Runkel
> https://arxiv.org/abs/1311.3354
"""

import ...Operations: ⨷, ⨶

function _orbifold_equivalence_def(::Type{TwoVariables.A{n}}, ::Type{TwoVariables.D{m}}, left_vars, right_vars) where {n, m}
    if n == 2m - 3
        # opposite convention to agree with notation in our source
        u, v = left_vars
        x, y = right_vars
        d = m - 1
        f = div(x^d - u^2d, x - u^2)
        return [0 x - u^2; f - y^2 0] ⨶ [0 v - u*y; v + u*y 0]
    else
        return false
    end
end
