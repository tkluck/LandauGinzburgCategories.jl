using LandauGinzburgCategories, LandauGinzburgCategories.Library
@ring! ℚ[x,y,s,t,u,v]
f = TwoVariables.A₂A₂; f(x, y)
g = TwoVariables.A₅; g(s, t)
h = TwoVariables.D₄; h(u, v)
A = orbifold_equivalence(f, g, [x, y], [s,t])
quantum_dimensions(A, [:x, :y], [:s, :t])
B = orbifold_equivalence(g, h, [s, t], [u, v])
quantum_dimensions(B, [:s, :t], [:u, :v])
AB = fuse(A, B, s, t)
quantum_dimensions(AB, [:x, :y], [:u, :v])
