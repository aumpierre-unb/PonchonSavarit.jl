@doc raw"""

`PonchonSavarit` provides a set of functions to compute 
the number of theoretical stages of a distillation column 
using the Ponchon-Savarit method.

See also: `stages`, `refmin`, `qR2S`.
"""
module PonchonSavarit

using Plots
using Test

export stages, refmin, qR2S

include("stages.jl")
include("refmin.jl")
include("qR2S.jl")
include("bissection.jl")
include("interp1.jl")
include("interp2.jl")

end