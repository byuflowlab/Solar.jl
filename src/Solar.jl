"""
    Solar

Defines type: panel, includes various methods for solar flux calculation.

"""
module Solar

using Interpolations
import CSV

include("solartypes.jl")
include("solarcapture.jl")

end # module
