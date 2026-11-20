test_that("the friction factor of the example is about 0.02012", {
  f <- colebrook_friction(1e5, 0.000045 / 0.10)
  expect_equal(f, 0.02012, tolerance = 1e-3)
})

test_that("all methods give the same friction factor", {
  ref <- colebrook_friction(1e5, 0.000045 / 0.10, "newton")
  for (m in c("bisection", "secant", "fixed_point", "uniroot")) {
    expect_equal(colebrook_friction(1e5, 0.000045 / 0.10, m), ref,
                 tolerance = 1e-6)
  }
})

test_that("the Reynolds number is correct", {
  expect_equal(reynolds(V = 1, D = 0.10), 1e5)
})

test_that("the head loss of the example is about 1.03 m", {
  r <- head_loss(D = 0.10, L = 100, V = 1)
  expect_equal(r[["hf"]], 1.025, tolerance = 1e-2)
})

test_that("a negative Reynolds number is rejected", {
  expect_error(colebrook_friction(-5, 0.001))
})
