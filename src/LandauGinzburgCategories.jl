@doc read(open(joinpath(@__DIR__, "..", "README.md")), String)
module LandauGinzburgCategories

include("QuasiHomogeneous.jl")
include("Operations.jl")
include("Library.jl")

import .Operations: ⨷, ⨶, fuse, fuse_abstract
import .Operations: dual, unit_matrix_factorization
import .QuasiHomogeneous: find_quasihomogeneous_degrees, quasidegree, centralcharge

export ⨷, ⨶, fuse, fuse_abstract
export dual, unit_matrix_factorization
export find_quasihomogeneous_degrees, quasidegree, centralcharge

end # module
