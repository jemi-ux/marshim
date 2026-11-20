f  <- function(x) x^2 - 2
df <- function(x) 2 * x

test_that("bisection finds the square root of 2", {
  r <- bisection(f, 1, 2)
  expect_equal(r[["root"]], sqrt(2), tolerance = 1e-8)
})

test_that("bisection rejects two endpoints of the same sign", {
  expect_error(bisection(function(x) x^2 + 1, 0, 1))
})

test_that("Newton finds the square root of 2 in few iterations", {
  r <- newton(f, df, 1)
  expect_equal(r[["root"]], sqrt(2), tolerance = 1e-8)
  expect_lte(r[["iterations"]], 6)
})

test_that("the secant method finds the square root of 2", {
  r <- secant(f, 1, 2)
  expect_equal(r[["root"]], sqrt(2), tolerance = 1e-8)
})

test_that("fixed point finds the solution of cos x = x", {
  r <- fixed_point(cos, 1)
  expect_equal(r[["root"]], 0.7390851332, tolerance = 1e-7)
})
