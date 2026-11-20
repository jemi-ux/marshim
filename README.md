# marshim

Solving nonlinear equations and an application to hydraulics.

R package developed by **Jemima Nkama**, as part of the introductory R
course, under the supervision of **Mr Yvan Ngneunmeu**.

## Contents

- `bisection()`, `fixed_point()`, `newton()`, `secant()`: four methods
  to solve f(x) = 0.
- `colebrook_friction()`: friction factor of a pipe
  (Colebrook-White equation).
- `head_loss()`: head loss in the pipe (Darcy-Weisbach).
- `reynolds()`, `colebrook_g()`, `colebrook_dg()`: helper functions.

## Installation

```r
# from the provided .tar.gz file
install.packages("marshim_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Example

```r
library(marshim)

newton(function(x) x^2 - 2, function(x) 2 * x, 1)

# 100 m steel pipe, 10 cm diameter, water at 1 m/s
head_loss(D = 0.10, L = 100, V = 1)
```

## Documentation

- Help for each function: `?newton`, `?colebrook_friction`, ...
- User guide: `vignette("introduction", package = "marshim")`
- Presentation scripts: `inst/scripts/` folder
