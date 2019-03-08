@doc read(open(joinpath(@__DIR__, "..", "README.md")), String)
module LandauGinzburgCategory

include("Operations.jl")
include("Library.jl")

import .Operations: ⨷, ⨶, fuse, unit_matrix_factorization
export ⨷, ⨶, fuse, unit_matrix_factorization

end # module
