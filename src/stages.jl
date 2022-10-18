include("interp1.jl")
include("interp2.jl")
include("bissection.jl")
include("refmin.jl")

using Plots

@doc raw"""
`N=stages(y,X,q,R[,fig])`

`stages` computes the number of theoretical stages
of a distillation column
using the Ponchon-Savarit method given
a x-h-y-H matrix of the liquid and the vapor fractions
at equilibrium and their enthalpies,
the vector of the fractions of the products and the feed,
the feed quality, and
the reflux ratio at the top of the column.

If feed is a saturated liquid, feed quality q = 1,
feed quality is reset to q = 1 - 1e-10.

By default, `stages` plots a schematic diagram of the solution, fig = true.
If fig = false is given, no plot is shown.

See also: `refmin`, `qR2S`.

Examples
==========
Compute the number of theoretical stages
of a distillation column for oxygen and nitrogen
from the bottom to the top of the column given
a matrix that relates the liquid and the vapor fractions
and their enthalpies at equilibrium,
the composition of the distillate is 88 %,
the composition of the feed is 46 %,
the composition of the column's bottom product is 11 %,
the feed quality is 54 %, and
the reflux ratio at the top of the column is
70 % higher that the minimum reflux ratio:

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
x=[0.88 0.46 0.11];
q=0.56;
r=refmin(data,x,q);
R=1.70*r;
N=stages(data,x,q,R)
```

Compute the number of theoretical stages
of a distillation column for acetone and methanol
from the bottom to the top of the column given
a matrix that relates the liquid and the vapor fractions
and their enthalpies at equilibrium,
the composition of the distillate is 88 %,
the composition of the feed is 46 %,
the composition of the column's bottom product is 11 %,
the feed is a saturated liquid, and
the reflux ratio at the top of the column is
70 % higher that the minimum reflux ratio:

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
      1.     2.250 1.       17.390];
x=[0.88 0.46 0.11];
q=1;
r=refmin(data,x,q);
R=1.70*r;
N=stages(data,x,q,R)
```
"""
function stages(data, X, q, R, fig=true)
    xD = X[1]
    xF = X[2]
    xB = X[3]
    if xD < xF || xB > xF
        println("Inconsistent feed and/or products compositions.")
        return
    end
    if q == 1
        q = 1 - 1e-10
    end
    if R <= refmin(data, X, q)
        println("Minimum reflux ratio exceeded.")
        return
    end
    f(x) = interp1(data[:, 3], data[:, 1], x)
    g(x) = interp1(data[:, 1], data[:, 2], x)
    k(x) = interp1(data[:, 3], data[:, 4], x)

    foo(x) = q / (q - 1) * x - xF / (q - 1)
    bar(x) = interp1(data[:, 1], data[:, 3], x) - foo(x)
    x1 = bissection(bar, xB, xD)
    h1 = g(x1)
    y1 = interp1(data[:, 1], data[:, 3], x1)
    H1 = k(y1)
    hF = (H1 - h1) * (1 - q) + h1

    hliq = g(xD)
    Hvap = k(xD)
    hdelta = (Hvap - hliq) * R + Hvap

    hlambda = (hdelta - hF) / (xD - xF) * (xB - xF) + hF

    x2 = interp2(g, X, [xD; hdelta], [xB; hlambda])

    y = [xD]
    x = [f(y[end])]
    while x[end] > xB
        if x[end] > x2
            P = [xD; hdelta]
        else
            P = [xB; hlambda]
        end
        Q = [x[end]; g(x[end])]
        y = [y; interp2(k, X, P, Q)]
        x = [x; f(y[end])]
    end

    x = [xD; x]
    y = [y; x[end]]

    h = g.(x)
    H = k.(y)

    if fig
        p1 = plot(xlabel="x,y", ylabel="h,H",
            xlims=(0, 1),
            legend=false,
            framestyle=:box,
            grid=:true,
            minorgrid=:true)

        plot!(data[:, 1], data[:, 2],
            seriescolor=:blue,
            linestyle=:solid,
            markershape=:diamond,
            markerstrokecolor=:blue,
            markersize=3)
        plot!(data[:, 3], data[:, 4],
            seriescolor=:red,
            linestyle=:solid,
            markershape=:diamond,
            markerstrokecolor=:red,
            markersize=3)
        plot!([xD; xD; xD; xF; xB; xB; xB],
            [g(xD); k(xD); hdelta; hF; hlambda; g(xB); k(xB)],
            seriescolor=:green,
            linestyle=:solid,
            markershape=:circle,
            markerstrokecolor=:green,
            markersize=3)
        plot!([x1 xF y1], [h1 hF H1],
            color=:magenta,
            linestyle=:dash)
        plot!(reshape([x y]'[1:end-1], (2 * size(x, 1) - 1, 1)),
            reshape([h H]'[1:end-1], (2 * size(x, 1) - 1, 1)),
            color=:cyan,
            linestyle=:solid)

        p2 = plot(xlabel="x", ylabel="y",
            xlims=(0, 1), ylims=(0, 1),
            legend=false,
            framestyle=:box,
            grid=:true,
            minorgrid=:true)

        X = data[:, 1]
        Y = data[:, 3]
        plot!(X, Y,
            seriescolor=:black,
            markershape=:diamond,
            markerstrokecolor=:black,
            markersize=3)

        X = [0; 1]
        Y = X
        plot!(X, Y,
            seriestype=:line, color=:black,
            linestyle=:dash)

        plot!([xD],
            seriestype=:vline, color=:blue,
            linestyle=:dash)

        plot!([xB],
            seriestype=:vline, color=:red,
            linestyle=:dash)

        if q != 1 - 1e-10
            Y = q / (q - 1) .* X .- xF / (q - 1)
            plot!(X, Y,
                linestyle=:dash, color=:magenta)
        end

        plot!([xF],
            seriestype=:vline,
            color=:magenta,
            linestyle=:dash)

        plot!(x, y,
            seriestype=:steppost, color=:cyan,
            linestyle=:solid)
        plot!(x, y,
            seriescolor=:green,
            linestyle=:solid,
            markershape=:circle,
            markerstrokecolor=:green,
            markersize=3)

        display(plot(layout=(2, 1), p1, p2,
            size=(500, 800),
            pos=(0,500),
            margin=5Plots.mm))
    end
    return size(x, 1) - 1 - 1 + (xB - x[end-1]) / (x[end] - x[end-1])
end

