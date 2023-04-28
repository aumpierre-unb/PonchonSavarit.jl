include("bissection.jl")
include("interp1.jl")

@doc raw"""
`qR2S(data::Matrix{Float64}, z::Vector{Float64}, q::Number, R::Number)`

`RS2q` computes the feed quality
of a distillation column
using the Ponchon-Savarit method given
a x-h-y-H matrix of the liquid and the vapor fractions
at equilibrium and their enthalpies,
the vector of the fractions of the products and the feed,
the reflux ratio at the top and
the reflux ratio at the bottom of the column.

If feed is a saturated liquid, feed quality q = 1,
feed quality is reset to q = 1 - 1e-10.

`RS2q` is a main function of
the `PonchonSavarit` toolbox for Julia.

See also: `stages`, `refmin`, `qR2S`, `qS2R`.

Examples
==========
Compute the reflux ratio at the bottom
of a distillation column for oxygen and nitrogen given
the composition of the distillate is 88 %,
the composition of the feed is 44 %,
the composition of the column's bottom product is 8 %,
the feed quality is 55 % and
the reflux ratio at the top of the column is 2:

```
data=[0.    0.420 0.    1.840; # enthalpy in kcal/mmol
      0.075 0.418 0.193 1.755;
      0.17  0.415 0.359 1.685;
      0.275 0.410 0.50  1.625;
      0.39  0.398 0.63  1.570;
      0.525 0.378 0.75  1.515;
      0.685 0.349 0.86  1.465;
      0.88  0.300 0.955 1.425;
      1.    0.263 1.    1.405];
x=[0.88;0.55;0.08];
q=RS2q(data,x,2,2.4)
```

Compute the reflux ratio at the bottom
of a distillation column for acetone and methanol given
the composition of the distillate is 88 %,
the composition of the feed is 44 %,
the composition of the column's bottom product is 8 %,
the feed is saturated liquid and
the reflux ratio at the top of the column is 2:

```
data=[2.5e-4 3.235 1.675e-3 20.720; # enthalpy in kcal/mol
      0.05   2.666 0.267    20.520;
      0.1    2.527 0.418    20.340;
      0.15   2.459 0.517    20.160;
      0.2    2.422 0.579    20.000;
      0.3    2.384 0.665    19.640;
      0.4    2.358 0.729    19.310;
      0.5    2.338 0.779    18.970;
      0.6    2.320 0.825    18.650;
      0.7    2.302 0.87     18.310;
      0.8    2.284 0.915    17.980;
      0.9    2.266 0.958    17.680;
      1      2.250 1        17.390];
x=[0.88;0.55;0.08];
q=RS2q(data,x,2,2.4)
```
"""
function RS2q(data::Matrix{Float64}, z::Vector{Float64}, R::Number, S::Number)
    xD, xF, xB = z
    if xD < xF || xB > xF
        error("Inconsistent feed and/or products compositions.")
    end
    x2h(x) = interp1(data[:, 1], data[:, 2], x)
    y2H(y) = interp1(data[:, 3], data[:, 4], y)
    x2y(x) = interp1(data[:, 1], data[:, 3], x)
    h2 = x2h(xD)
    H2 = y2H(xD)
    hdelta = (H2 - h2) * R + H2
    h3 = x2h(xB)
    H3 = y2H(xB)
    hlambda = (h3 - H3) * S + h3
    hF = (hdelta - hlambda) / (xD - xB) * (xF - xB) + hlambda
    foo(x) = (y2H(x2y(x)) - x2h(x)) / (x2y(x) - x) - (y2H(x2y(x)) - hF) / (x2y(x) - xF)
    x1 = newtonraphson(foo, xF)
    h1 = x2h(x1)
    y1 = x2y(x1)
    H1 = y2H(y1)
    (y1 - xF) / (y1 - x1)
end
