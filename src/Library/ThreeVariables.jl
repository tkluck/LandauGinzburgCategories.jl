module ThreeVariables

import ...Library: Potential, subscript, leftvars, rightvars

struct E₁₃ <: Potential{3}
    E₁₃() = E₁₃(leftvars(3)...)
    E₁₃(x, y, z) = y^3 + y*x^5 + z^2
end

struct E₁₄{variant} <: Potential{3}
    E₁₄() = E₁₄{:v1}()
    E₁₄{variant}() where variant = E₁₄{variant}(leftvars(3)...)
    E₁₄{:v1}(x, y, z) = x^4*z + y^3 + z^2
    E₁₄{:v2}(x, y, z) = x^8   + y^3 + z^2
end

struct Q₁₀ <: Potential{3}
    Q₁₀() = Q₁₀(leftvars(3)...)
    Q₁₀(x, y, z) = x^4 + y^3 + x*z^2
end

struct Q₁₁ <: Potential{3}
    Q₁₁() = Q₁₁(leftvars(3)...)
    Q₁₁(x, y, z) = y*z^3 + y^3 + x^2 * z
end

struct Q₁₂{variant} <: Potential{3}
    Q₁₂() = E₁₄{:v1}()
    Q₁₂{variant}() where variant = Q₁₂{variant}(leftvars(3)...)
    Q₁₂{:v1}(x, y, z) = x^3*z + y^3   + x*z^2
    Q₁₂{:v2}(x, y, z) = x^5   + y^3   + x*z^2
end

struct S₁₁ <: Potential{3}
    S₁₁() = S₁₁(leftvars(3)...)
    S₁₁(x, y, z) = x^2*z + y*z^2 + y^4
end

struct S₁₂ <: Potential{3}
    S₁₂() = S₁₂(leftvars(3)...)
    S₁₂(x, y, z) = x^3*y + y^2*z + x*z^2
end

struct U₁₂{variant} <: Potential{3}
    U₁₂() = U₁₂{:v1}()
    U₁₂{variant}() where variant = U₁₂{variant}(leftvars(3)...)
    U₁₂{:v1}(x, y, z) = x^4 + y^3   + z^3
    U₁₂{:v2}(x, y, z) = x^4 + y^3   + z^2*y
    U₁₂{:v3}(x, y, z) = x^4 + y^2*z + z^2*y
end

struct W₁₂{variant} <: Potential{3}
    W₁₂() = W₁₂{:v1}()
    W₁₂{variant}() where variant = W₁₂{variant}(leftvars(3)...)
    W₁₂{:v1}(x, y, z) = x^5   + y^2*z + z^2
    W₁₂{:v2}(x, y, z) = x^5   + y^4   + z^2
end

struct W₁₃{variant} <: Potential{3}
    W₁₃() = W₁₃{:v1}()
    W₁₃{variant}() where variant = W₁₃{variant}(leftvars(3)...)
    W₁₃{:v1}(x, y, z) = x^2  + y^4   + y*z^4
    W₁₃{:v2}(x, y, z) = y*x^4 + y^2*z + z^2
    W₁₃{:v3}(x, y, z) = x^4*y + y^4   + z^2
end

struct Z₁₂ <: Potential{3}
    Z₁₂() = Z₁₂(leftvars(3)...)
    Z₁₂(x, y, z) = y*x^4  + x*y^3 + z^2
end

struct Z₁₃{variant} <: Potential{3}
    Z₁₃() = Z₁₃{:v1}()
    Z₁₃{variant}() where variant = Z₁₃{variant}(leftvars(3)...)
    Z₁₃{:v1}(x, y, z) = x^6   + y^3*x + z^2
    Z₁₃{:v2}(x, y, z) = x^3*z + x*y^3 + z^2
end

end
