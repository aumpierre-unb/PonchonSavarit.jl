@doc raw"""
`newtonraphson` computes computes the root of
a function using the method of Newton-Raphson
from the initial guess values.

`newtonraphson` is an internal function of
the `PonchonSavarit` toolbox.
"""
function newtonraphson(f, x)
    while abs(f(x)) > 5e-4
        a = (f(x + 1e-7) - f(x)) / 1e-7
        if x - f(x) / a > 1
            x = (x + 1) / 2
        elseif x - f(x) / a < 0
            x = x / 2
        else
            x = x - f(x) / a
        end
    end
    x
end