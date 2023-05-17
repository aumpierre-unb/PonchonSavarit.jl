@doc raw"""
`PonchonSavarit` provides a set of functions to compute 
the number of theoretical stages of a distillation column 
using the Ponchón-Savarit method.

Author: Alexandre Umpierre `aumpierre@gmail.com`

Maintainer's repository: `https://github.com/aumpierre-unb/Ponchón-Savarit.jl`

Citation (any version): `DOI 10.5281/zenodo.7218007`

See also: `stages`, `refmin`, `qR2S`, `qS2R`, `RS2q`.
"""
module PonchonSavarit

using Plots
using Test

export stages, refmin, qR2S, qS2R, RS2q

include("stages.jl")
include("refmin.jl")
include("qR2S.jl")
include("qS2R.jl")
include("RS2q.jl")
include("bissection.jl")
include("newtonraphson.jl")
include("interp1.jl")
include("interp2.jl")

end