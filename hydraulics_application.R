# =====================================================================
# APPLICATION: HYDRAULICS (civil engineering)
# Reproduces the two charts of the presentation:
#   1) the friction factor f as a function of the Reynolds number
#   2) the head loss as a function of the pipe diameter
#
# Each point requires solving the Colebrook-White equation, which is
# nonlinear: we use R's uniroot().
# =====================================================================


# ---------------------------------------------------------------------
# STEP 1: function that computes f (Colebrook equation)
# ---------------------------------------------------------------------
# Re  : Reynolds number
# rel : relative roughness = eps / D
# We set x = 1/sqrt(f) and look for x such that g(x) = 0.
colebrook <- function(Re, rel) {
  g <- function(x) x + 2 * log10(rel / 3.7 + 2.51 * x / Re)
  x <- uniroot(g, c(3, 30), tol = 1e-12)$root   # root between 3 and 30
  1 / x^2                                        # back to f
}

# Test: steel (eps = 0.045 mm), D = 0.10 m, Re = 100000
colebrook(1e5, 0.000045 / 0.10)                  # gives about 0.02012


# ---------------------------------------------------------------------
# STEP 2: f as a function of the Reynolds number
# ---------------------------------------------------------------------
Re_vals <- c(4000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000)

f_smooth <- numeric(length(Re_vals))   # smooth pipe
f_steel  <- numeric(length(Re_vals))   # new steel   (eps/D = 0.00045)
f_rough  <- numeric(length(Re_vals))   # very rough  (eps/D = 0.005)

for (i in 1:length(Re_vals)) {
  f_smooth[i] <- colebrook(Re_vals[i], 0)
  f_steel[i]  <- colebrook(Re_vals[i], 0.000045 / 0.10)
  f_rough[i]  <- colebrook(Re_vals[i], 0.005)
}

tab_f <- data.frame(Re = Re_vals, smooth = f_smooth,
                    steel = f_steel, rough = f_rough)
print(tab_f, digits = 4)


# ---------------------------------------------------------------------
# STEP 3: head loss as a function of the diameter
# Fixed flow rate: that of the example (V = 1 m/s in the 0.10 m pipe)
# ---------------------------------------------------------------------
L    <- 100                       # length (m)
eps  <- 0.000045                  # steel roughness (m)
nu   <- 1e-6                      # kinematic viscosity of water (m2/s)
grav <- 9.81                      # gravity (m/s2)
Q    <- pi / 4 * 0.10^2 * 1       # flow rate (m3/s) = 7.85 L/s

diam <- c(0.05, 0.075, 0.10, 0.125, 0.15, 0.20)   # diameters tested (m)

velocity <- numeric(length(diam))
hf       <- numeric(length(diam))

for (i in 1:length(diam)) {
  d <- diam[i]
  velocity[i] <- 4 * Q / (pi * d^2)              # V = Q / cross-section
  Re <- velocity[i] * d / nu                     # Reynolds number
  f  <- colebrook(Re, eps / d)                   # friction factor
  hf[i] <- f * (L / d) * velocity[i]^2 / (2 * grav)   # Darcy-Weisbach
}

tab_hf <- data.frame(Diameter_cm = diam * 100,
                     Velocity_m_s = round(velocity, 2),
                     Loss_m       = round(hf, 2))
print(tab_hf)

cat("\nDiameter divided by 2 (10 cm -> 5 cm): loss multiplied by",
    round(hf[1] / hf[3]), "\n")


# ---------------------------------------------------------------------
# STEP 4: charts
# ---------------------------------------------------------------------
par(mfrow = c(1, 2))

# Chart 1: f as a function of Re (logarithmic scale on x)
plot(Re_vals, f_rough, type = "b", log = "x", col = "orange", pch = 16,
     ylim = c(0, 0.05), xlab = "Reynolds number Re",
     ylab = "Friction factor f",
     main = "Friction depending on the flow")
lines(Re_vals, f_steel,  type = "b", col = "darkcyan", pch = 17)
lines(Re_vals, f_smooth, type = "b", col = "navy",     pch = 15)
legend("topright", bty = "n", pch = c(15, 17, 16),
       col = c("navy", "darkcyan", "orange"),
       legend = c("Smooth pipe", "New steel", "Very rough"))

# Chart 2: head loss depending on the diameter
barplot(hf, names.arg = paste(diam * 100, "cm"), col = "darkcyan",
        xlab = "Pipe diameter",
        ylab = "Head loss over 100 m (m of water)",
        main = "Choosing the diameter")
