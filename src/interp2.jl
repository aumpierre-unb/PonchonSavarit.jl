include("bissection.jl")

function interp2(f, X, P, Q)
    g(x) = (Q[2] - P[2]) / (Q[1] - P[1]) * (x - P[1]) + P[2]
    h(x) = f(x) - g(x)
    return bissection(h, X[3], X[1])
end
