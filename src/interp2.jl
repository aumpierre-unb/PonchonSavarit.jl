include("bissection.jl")

@doc raw"""
`interp2` computes the interpolation of a number.

`interp2` is an internal function of
the `PonchonSavarit` toolbox.
"""
function interp2(f, X, P, Q)
    g(x) = (Q[2] - P[2]) / (Q[1] - P[1]) * (x - P[1]) + P[2]
    h(x) = f(x) - g(x)
    return bissection(h, X[3], X[1])
end
