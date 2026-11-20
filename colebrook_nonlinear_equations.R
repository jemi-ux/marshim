# =====================================================================
# PRESENTATION: NONLINEAR EQUATIONS  (introductory R course)
# Application: civil engineering / hydraulics
#
# Problem: compute the friction factor f of a water pipe using the
# Colebrook-White equation:
#
#   1/sqrt(f) = -2 * log10( eps/(3.7*D) + 2.51/(Re*sqrt(f)) )
#
# f is on both sides: it cannot be isolated, so a numerical method is
# needed. We set x = 1/sqrt(f) and look for x such that g(x) = 0.
# =====================================================================


# ---------------------------------------------------------------------
# STEP 1: the data (variables)
# ---------------------------------------------------------------------
D   <- 0.10      # pipe diameter (m)
eps <- 0.000045  # steel roughness (m)
Re  <- 100000    # Reynolds number (turbulent flow)
tol <- 1e-10     # requested precision


# ---------------------------------------------------------------------
# STEP 2: the function g(x) and its derivative dg(x)
# ---------------------------------------------------------------------
# g(x) = 0  <=>  Colebrook equation with x = 1/sqrt(f)
g <- function(x) {
  x + 2 * log10(eps / (3.7 * D) + 2.51 * x / Re)
}

# dg(x): derivative of g, needed for Newton's method
dg <- function(x) {
  1 + (2 / log(10)) * (2.51 / Re) / (eps / (3.7 * D) + 2.51 * x / Re)
}

# phi(x): we rewrite g(x) = 0 in the form x = phi(x) (fixed point)
phi <- function(x) {
  -2 * log10(eps / (3.7 * D) + 2.51 * x / Re)
}


# ---------------------------------------------------------------------
# STEP 3: the four methods, written as R functions
# Each function returns the root and the number of iterations.
# ---------------------------------------------------------------------

# --- Method 1: bisection ----------------------------------------------
# We cut the interval [a, b] in half until the precision is reached.
bisection <- function(f, a, b, tol) {
  k <- 0                              # iteration counter
  repeat {
    k <- k + 1
    m <- (a + b) / 2                  # midpoint of the interval
    if (abs(f(m)) < tol || (b - a) / 2 < tol) break   # precision reached
    if (f(a) * f(m) < 0) {            # the root is in [a, m]
      b <- m
    } else {                          # the root is in [m, b]
      a <- m
    }
  }
  c(root = m, iterations = k)
}

# --- Method 2: fixed point --------------------------------------------
# We repeat x <- phi(x) until x almost stops changing.
fixed_point <- function(phi, x0, tol) {
  x <- x0
  k <- 0
  repeat {
    k <- k + 1
    x_new <- phi(x)
    if (abs(x_new - x) < tol) break
    x <- x_new
  }
  c(root = x_new, iterations = k)
}

# --- Method 3: Newton-Raphson -----------------------------------------
# x_new = x - f(x) / f'(x)   (we follow the tangent of the curve)
newton <- function(f, df, x0, tol) {
  x <- x0
  k <- 0
  repeat {
    k <- k + 1
    x_new <- x - f(x) / df(x)
    if (abs(x_new - x) < tol) break
    x <- x_new
  }
  c(root = x_new, iterations = k)
}

# --- Method 4: secant -------------------------------------------------
# Like Newton, but the derivative is replaced by a slope between
# two points: no need to compute f'.
secant <- function(f, x0, x1, tol) {
  k <- 0
  repeat {
    k <- k + 1
    x2 <- x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
    if (abs(x2 - x1) < tol) break
    x0 <- x1
    x1 <- x2
  }
  c(root = x2, iterations = k)
}


# ---------------------------------------------------------------------
# STEP 4: we run the four methods
# ---------------------------------------------------------------------
r1 <- bisection(g, 3, 12, tol)        # we search between x = 3 and x = 12
r2 <- fixed_point(phi, 8, tol)        # starting point x0 = 8
r3 <- newton(g, dg, 8, tol)           # starting point x0 = 8
r4 <- secant(g, 6, 10, tol)           # two starting points

# Reference solution given directly by R
x_ref <- uniroot(g, c(3, 12), tol = 1e-12)$root


# ---------------------------------------------------------------------
# STEP 5: comparison table (data.frame)
# ---------------------------------------------------------------------
x_sol  <- unname(c(r1["root"], r2["root"], r3["root"], r4["root"]))
n_iter <- unname(c(r1["iterations"], r2["iterations"],
                   r3["iterations"], r4["iterations"]))

results <- data.frame(
  Method     = c("Bisection", "Fixed point", "Newton", "Secant"),
  x          = x_sol,
  f          = 1 / x_sol^2,             # back to f = 1/x^2
  Iterations = n_iter
)

print(results, digits = 8)
cat("\nReference solution (uniroot): x =", x_ref,
    "  f =", 1 / x_ref^2, "\n")


# ---------------------------------------------------------------------
# STEP 6: checking the Colebrook equation
# ---------------------------------------------------------------------
f_sol <- 1 / x_sol[3]^2                 # friction factor (Newton)
lhs   <- 1 / sqrt(f_sol)
rhs   <- -2 * log10(eps / (3.7 * D) + 2.51 / (Re * sqrt(f_sol)))
cat("\nCheck: 1/sqrt(f) =", lhs,
    "  right-hand side =", rhs, "\n")
cat("Difference =", abs(lhs - rhs), "(close to 0: the equation is satisfied)\n")


# ---------------------------------------------------------------------
# STEP 7: application, the head loss (Darcy-Weisbach)
# ---------------------------------------------------------------------
L    <- 100                 # pipe length (m)
V    <- 1                   # water velocity (m/s)
grav <- 9.81                # gravity (m/s2)

hf <- f_sol * (L / D) * V^2 / (2 * grav)
cat("\nHead loss over", L, "m: hf =", round(hf, 2), "m of water\n")


# ---------------------------------------------------------------------
# STEP 8: charts
# ---------------------------------------------------------------------
par(mfrow = c(1, 2))        # two charts side by side

# Chart 1: the curve g(x). The root is where it crosses the x axis.
curve(g(x), from = 3, to = 12, col = "blue", lwd = 2,
      xlab = "x = 1/sqrt(f)", ylab = "g(x)",
      main = "Colebrook equation: g(x) = 0")
abline(h = 0, lty = 2)                       # dashed horizontal axis
points(x_ref, 0, pch = 19, col = "red", cex = 1.5)

# Chart 2: number of iterations of each method
barplot(n_iter, names.arg = results$Method,
        col = c("orange", "darkgreen", "red", "purple"),
        ylab = "Number of iterations",
        main = "Which method converges fastest?")
