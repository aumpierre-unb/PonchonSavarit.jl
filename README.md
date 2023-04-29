# PonchonSavarit.jl

[![DOI](https://zenodo.org/badge/553159585.svg)](https://zenodo.org/badge/latestdoi/553159585)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![version](https://juliahub.com/docs/PonchonSavarit/version.svg)](https://juliahub.com/ui/Packages/PonchonSavarit/WauTj)

## Installing and Loading PonchonSavarit

PonchonSavarit can be installed and loaded either
from the JuliaHub repository (last released version) or from the
[maintainer's repository](https://github.com/aumpierre-unb/PonchonSavarit.jl).

### Last Released Version

The last version of PonchonSavarit can be installed from JuliaHub repository:

```julia
using Pkg
Pkg.add("PonchonSavarit")
using PonchonSavarit
```

If PonchonSavarit is already installed, it can be updated:

```julia
using Pkg
Pkg.update("PonchonSavarit")
using PonchonSavarit
```

### Pre-Release (Under Construction) Version

The pre-release (under construction) version of PonchonSavarit
can be installed from the [maintainer's repository](https://github.com/aumpierre-unb/PonchonSavarit.jl).

```julia
using Pkg
Pkg.add(path="https://github.com/aumpierre-unb/PonchonSavarit.jl")
using PonchonSavarit
```

## Citation of PonchonSavarit

You can cite all versions (both released and pre-released), by using
[10.5281/zenodo.7218007](https://doi.org/10.5281/zenodo.7218007).
This DOI represents all versions, and will always resolve to the latest one.

## The PonchonSavarit Module for Julia

PonchonSavarit provides the following functions:

- stages
- refmin
- qR2S
- qS2R
- RS2q

### stages

stages computes the number of theoretical stages of a distillation column using the Ponchon-Savarit method given a x-h-y-H matrix of the liquid and the vapor fractions at equilibrium and their enthalpies, the vector of the fractions of the products and the feed, the feed quality, and the reflux ratio at the top of the column.

If feed is a saturated liquid, feed quality q = 1, feed quality is reset to q = 1 - 1e-10.

By default, fig = true, `stages` plots a schematic diagram of the solution. If fig = false is given, no plot is shown.

**Syntax:**

```julia
stages(data::Union{Matrix{Float64},Function},z::Vector{Float64};q::Number=NaN,R::Number=NaN,S::Number=NaN,fig::Bool=true)
```

**Examples:**

Compute the number of theoretical stages of a distillation column for oxygen and nitrogen from the bottom to the top of the column given a matrix that relates the liquid and the vapor fractions and their enthalpies at equilibrium, the composition of the distillate is 88 %, the composition of the feed is 44 %, the composition of the bottoms is 8 %, the feed quality is 55 % and the reflux ratio at the top of the column is 70 % higher that the minimum reflux ratio:

```julia
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
r=refmin(data,x,q=0.55)[1];
N=stages(data,x,q=0.55,R=1.70*r,fig=false)
```

Compute the number of theoretical stages of a distillation column for oxygen and nitrogen from the bottom to the top of the column given a matrix that relates the liquid and the vapor fractions and their enthalpies at equilibrium, the composition of the distillate is 92 %, the composition of the feed is 59 %, the composition of the bottoms is 9 %, the feed quality is 37 % and the reflux ratio at the bottom of the column is 40 % higher that the minimum reflux ratio, and plot a schematic diagram of the solution:

```julia
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
r,s=refmin(data,x,q=0.37)
N=stages(data,x,q=0.37,S=1.40*s)
```

Compute the number of theoretical stages of a distillation column for acetone and methanol from the bottom to the top of the column given a matrix that relates the liquid and the vapor fractions and their enthalpies at equilibrium, the composition of the distillate is 88 %, the composition of the feed is 44 %, the composition of the column's bottom product is 11 %, the composition of the bottoms is 8 %, the feed is a saturated liquid and the reflux ratio at the top of the column is 70 % higher that the minimum reflux ratio, and plot a schematic diagram of the solution:

```julia
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
x=[0.88;0.46;0.11];
r=refmin(data,x)[1];
N=stages(data,x,q=1,R=1.70*r)
```

### refmin

refmin computes the minimum value of the reflux ratio of a distillation column using the Ponchón-Savarit method given a x-h-y-H matrix of the liquid and the vapor fractions at equilibrium and their enthalpies, the vector of the fractions of the products and the feed, and the feed quality.

If feed is a saturated liquid, feed quality q = 1, feed quality is reset to q = 1 - 1e-10.

**Syntax:**

```julia
r=refmin(y,X,q)
```

**Examples:**

Compute the minimum value of the reflux ratio of a distillation column for oxygen and nitrogen given a matrix that relates the liquid and the vapor fractions and their enthalpies at equilibrium, the composition of the distillate is 88 %, the composition of the feed is 44 %, the composition of the bottoms is 8 % and the feed quality is 55 %:

```julia
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
r,s=refmin(data,x,q=0.55)
```

Compute the minimum value of the reflux ratio of a distillation column for acetone and methanol given a matrix that relates the liquid and the vapor fractions and their enthalpies at equilibrium, the composition of the distillate is 88 %, the composition of the feed is 44 %, the composition of the bottoms is 8 % and the feed is a saturated liquid.

```julia
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
x=[0.88;0.55;0.08];
r,s=refmin(data,x)
```

### qR2S

qR2S computes the reflux ratio at the bottom of a distillation column using the Ponchon-Savarit method given a x-h-y-H matrix of the liquid and the vapor fractions at equilibrium and their enthalpies, the vector of the fractions of the products and the feed, the feed quality, and the reflux ratio at the top of the column.

If feed is a saturated liquid, feed quality q = 1, feed quality is reset to q = 1 - 1e-10.

**Syntax:**

```julia
S=qR2S(X,q,R)
```

**Examples:**

Compute the reflux ratio at the bottom of a distillation column for oxygen and nitrogen given the composition of the distillate is 88 %, the composition of the feed is 44 %, the composition of the bottoms is 8 %, the feed quality is 55 % and the reflux ratio at the top of the column is 2:

```julia
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
S=qR2S(data,x,0.54,2)
```

Compute the reflux ratio at the bottom of a distillation column for acetone and methanol given the composition of the distillate is 88 %, the composition of the feed is 44 %, the composition of the bottoms is 8 %, the feed is saturated liquid, and the reflux ratio at the top of the column is 2:

```julia
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
x=[0.88;0.55;0.08];
S=qR2S(data,x,1,2)
```

### qS2R

qS2R computes the reflux ratio at the top of a distillation column using the Ponchon-Savarit method given a x-h-y-H matrix of the liquid and the vapor fractions at equilibrium and their enthalpies, the vector of the fractions of the products and the feed, the feed quality and the reflux ratio at the bottom of the column.

If feed is a saturated liquid, feed quality q = 1, feed quality is reset to q = 1 - 1e-10.

**Syntax:**

```julia
S=qR2S(X,q,R)
```

**Examples:**

Compute the reflux ratio at the bottom of a distillation column for oxygen and nitrogen given the composition of the distillate is 88 %, the composition of the feed is 44 %, the composition of the column's bottom product is 8 %, the feed quality is 55 % and the reflux ratio at the top of the column is 2:

```julia
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
R=qS2R(data,x,0.54,2.4)
```

Compute the reflux ratio at the bottom of a distillation column for acetone and methanol given the composition of the distillate is 88 %, the composition of the feed is 44 %, the composition of the column's bottom product is 8 %, the feed is saturated liquid and the reflux ratio at the top of the column is 2:

```julia
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
R=qS2R(data,x,1,2.4)
```

### qR2S

RS2q computes the feed quality
of a distillation column using the Ponchon-Savarit method given a x-h-y-H matrix of the liquid and the vapor fractions at equilibrium and their enthalpies, the vector of the fractions of the products and the feed, the reflux ratio at the top and the reflux ratio at the bottom of the column.

If feed is a saturated liquid, feed quality q = 1, feed quality is reset to q = 1 - 1e-10.

**Syntax:**

```julia
S=qR2S(X,q,R)
```

**Examples:**

Compute the reflux ratio at the bottom of a distillation column for oxygen and nitrogen given the composition of the distillate is 88 %, the composition of the feed is 44 %, the composition of the column's bottom product is 8 %, the feed quality is 55 % and the reflux ratio at the top of the column is 2:

```julia
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

Compute the reflux ratio at the bottom of a distillation column for acetone and methanol given the composition of the distillate is 88 %, the composition of the feed is 44 %, the composition of the column's bottom product is 8 %, the feed is saturated liquid and the reflux ratio at the top of the column is 2:

```julia
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

### See Also

[McCabeThiele.jl](https://github.com/aumpierre-unb/McCabeThiele.jl),
[Psychrometrics.jl](https://github.com/aumpierre-unb/Psychrometrics.jl),
[InternalFluidFlow.jl](https://github.com/aumpierre-unb/InternalFluidFlow.jl).

Copyright &copy; 2022 2023 Alexandre Umpierre

email: <aumpierre@gmail.com>
