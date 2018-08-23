"""
    Solar

Defines type: panel, includes various methods for solar flux calculation.

"""
module Solar

__precompile__()

using GeometricalPredicates
using Interpolations
using LinearAlgebra
using Statistics
import CSV

include("solartypes.jl")
include("solarcapture.jl")

end # module
