include("interp1.jl")
include("interp2.jl")
include("newtonraphson.jl")
include("refmin.jl")

using Plots

@doc raw"""
`stages(data::Union{Matrix{Float64},Function}, z::Vector{Float64}; q::Number=NaN, R::Number=NaN, S::Number=NaN, fig::Bool=true)`

`stages` computes the number of theoretical stages
of a distillation column
using the Ponchon-Savarit method
from the top to the bottom of the column given
a x-h-y-H matrix of the liquid and the vapor fractions
at equilibrium and their enthalpies,
the vector of the fractions of the products and the feed,
the feed quality, and
the reflux ratio at the top of the column.

By default, fig = true, `stages` plots a schematic diagram of the solution.
If fig = false is given, no plot is shown.

`stages` is a main function of
the `PonchonSavarit` toolbox for Julia.

See also: `refmin`, `qR2S`, `qS2R`, `RS2q`.

Examples
==========
Compute the number of theoretical stages
of a distillation column for acetone and methanol
from the bottom to the top of the column given
a matrix that relates the liquid and the vapor fractions
and their enthalpies at equilibrium,
the composition of the distillate is 93 %,
the composition of the feed is 41 %,
the composition of the bottoms is 7 %,
the feed is a saturated liquid and
the reflux ratio at the top of the column is
55 % higher that the minimum reflux ratio,
and plot a schematic diagram of the solution:

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
x=[0.93;0.41;0.07];
r,s=refmin(data,x)
N=stages(data,x,q=1,R=1.55*r)
```

Compute the number of theoretical stages
of a distillation column for oxygen and nitrogen
from the bottom to the top of the column given
a matrix that relates the liquid and the vapor fractions
and their enthalpies at equilibrium,
the composition of the distillate is 88 %,
the composition of the feed is 44 %,
the composition of the bottoms is 8 %,
the feed quality is 55 % and
the reflux ratio at the bottom of the column is
46 % higher that the minimum reflux ratio and
plot a schematic diagram of the solution:

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
x=[0.88;0.44;0.08];
r,s=refmin(data,x,q=0.55)
N=stages(data,x,q=0.55,S=1.46*s)
```

Compute the number of theoretical stages
of a distillation column for oxygen and nitrogen
from the bottom to the top of the column given
a matrix that relates the liquid and the vapor fractions
and their enthalpies at equilibrium,
the composition of the distillate is 92 %,
the composition of the feed is 59 %,
the composition of the bottoms is 9 %,
the reflux ratio at the bottom of the column is 1.9,
the reflux ratio at the top of the column is 1.5:

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
x=[0.92;0.59;0.09];
q=RS2q(data,x,1.5,1.9)
N=stages(data,x,R=1.5,S=1.9,fig=false)
```
"""
function stages(data::Union{Matrix{Float64},Function}, z::Vector{Float64}; q::Number=NaN, R::Number=NaN, S::Number=NaN, fig::Bool=true)
    xD, xF, xB = z
    if xD < xF || xB > xF
        error("Inconsistent feed and/or products compositions.")
    end
    a = isnan.([q, R, S]) .!= 1
    if sum(a) != 2
        error("""stages requires that two parameter among
        the feed quality,
        the reflux ratio at the top of the column and
        the reflux ratio at the bottom of the column
        be given alone.""")
    end
    if a == [1, 0, 1]
        R = qS2R(data, z, q, S)
    elseif a == [0, 1, 1]
        q = RS2q(data, z, R, S)
    end
    r = refmin(data, z, q=q)[1]
    if R <= r
        error("Minimum reflux ratios exceeded.")
    end
    x2y(x) = interp1(data[:, 1], data[:, 3], x)
    y2x(y) = interp1(data[:, 3], data[:, 1], y)
    x2h(x) = interp1(data[:, 1], data[:, 2], x)
    y2H(x) = interp1(data[:, 3], data[:, 4], x)
    foo(x) = q - (x2y(x) - xF) / (x2y(x) - x)
    x1 = newtonraphson(foo, xB)
    h1 = x2h(x1)
    y1 = x2y(x1)
    H1 = y2H(y1)
    hF = (H1 - h1) / (y1 - x1) * (xF - x1) + h1
    h2 = x2h(xD)
    H2 = y2H(xD)
    hdelta = (H2 - h2) * R + H2
    hlambda = (hdelta - hF) / (xD - xF) * (xB - xF) + hF
    g(x) = (H1 - h1) / (y1 - x1) * (x - x1) + h1
    xi = interp2(g, z, [xD; hdelta], [xB; hlambda])
    y = [xD]
    x = [y2x(y[end])]
    while x[end] > xB
        if x[end] > x1
            P = [xD; hdelta]
        else
            P = [xB; hlambda]
        end
        Q = [x[end]; x2h(x[end])]
        y = [y; interp2(y2H, z, P, Q)]
        x = [x; y2x(y[end])]
    end
    x = [xD; x]
    y = [y; x[end]]
    h = x2h.(x)
    H = y2H.(y)
    if !fig
        return size(x, 1) - 1 - 1 + (xB - x[end-1]) / (x[end] - x[end-1])
    end
    plot1 = plot(xlabel="x,y", ylabel="h,H",
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
        [x2h(xD); y2H(xD); hdelta; hF; hlambda; x2h(xB); y2H(xB)],
        seriescolor=:green,
        linestyle=:solid,
        markershape=:circle,
        markerstrokecolor=:green,
        markersize=3)
    plot!([x1; xF; y1], [h1; hF; H1],
        linestyle=:dash,
        color=:magenta,
        markershape=:square,
        markerstrokecolor=:magenta,
        markersize=3)
    plot!(reshape([x y]'[1:end-1], (2 * size(x, 1) - 1, 1)),
        reshape([h H]'[1:end-1], (2 * size(x, 1) - 1, 1)),
        color=:cyan,
        linestyle=:solid)
    plot2 = plot(xlabel="x", ylabel="y",
        xlims=(0, 1), ylims=(0, 1),
        legend=false,
        framestyle=:box,
        grid=:true,
        minorgrid=:true)
    plot!(data[:, 1], data[:, 3],
        seriescolor=:black,
        markershape=:diamond,
        markerstrokecolor=:black,
        markersize=3)
    plot!([0; 1], [0; 1],
        seriestype=:line, color=:black,
        linestyle=:dash)
    plot!([xD],
        seriestype=:vline, color=:blue,
        linestyle=:dash)
    plot!([xB],
        seriestype=:vline, color=:red,
        linestyle=:dash)
    if q != 1 - 1e-10
        plot!([0; 1], q / (q - 1) .* [0; 1] .- xF / (q - 1),
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
    plot!(layout=(2, 1), plot1, plot2,
        size=(600, 900),
        margin=5Plots.mm)
    display(plot!())
    size(x, 1) - 1 - 1 + (xB - x[end-1]) / (x[end] - x[end-1])
end
