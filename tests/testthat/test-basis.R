test_that("make_knots incluye siempre el minimo de A (Q-17)", {
  set.seed(42)
  A <- rnorm(100, mean = 5, sd = 2)
  knots <- make_knots(A, 10)
  expect_equal(min(knots), min(A))
})

test_that("spline_basis k=0 genera indicadores correctamente", {
  x <- c(1, 2, 3, 4, 5)
  knots <- c(2, 4)
  Phi <- spline_basis(x, knots, k = 0L)
  
  expect_equal(ncol(Phi), 3) # intercepto + 2 nudos
  expect_equal(Phi[, 1], rep(1, 5))
  expect_equal(Phi[, 2], c(0, 1, 1, 1, 1)) # I(x >= 2)
  expect_equal(Phi[, 3], c(0, 0, 0, 1, 1)) # I(x >= 4)
})

test_that("spline_basis k=1 incluye termino lineal puro (Q-01)", {
  x <- c(1, 2, 3, 4, 5)
  knots <- c(2, 4)
  Phi <- spline_basis(x, knots, k = 1L)
  
  expect_equal(ncol(Phi), 4) # intercepto + x + 2 nudos
  expect_equal(Phi[, 1], rep(1, 5))
  expect_equal(Phi[, 2], x) # Termino lineal puro (borde 0)
  expect_equal(Phi[, 3], pmax(0, x - 2)) # (x - 2)_+
  expect_equal(Phi[, 4], pmax(0, x - 4)) # (x - 4)_+
})
