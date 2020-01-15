@doc read(open(joinpath(@__DIR__, "..", "README.md")), String)
module LandauGinzburgCategories

include("QuasiHomogeneous.jl")
include("MatrixUtil.jl")
include("Operations.jl")
include("OrbifoldEquivalence.jl")
include("Library.jl")

import .Operations: ⨷, ⨶, fuse, fuse_abstract
import .Operations: dual, unit_matrix_factorization
import .Operations: getpotential
import .OrbifoldEquivalence: quantum_dimension, quantum_dimensions
import .QuasiHomogeneous: find_quasihomogeneous_degrees, quasidegree, centralcharge

export ⨷, ⨶, fuse, fuse_abstract
export dual, unit_matrix_factorization
export getpotential
export quantum_dimension, quantum_dimensions
export find_quasihomogeneous_degrees, quasidegree, centralcharge

end # module
