test_that("make_knots devuelve nudos estrictamente interiores (Q-17 resuelta)", {
  set.seed(42)
  A <- rnorm(100, mean = 5, sd = 2)
  knots <- make_knots(A, 10)
  # Los nudos deben estar estrictamente dentro del rango para evitar colinealidad
  expect_true(min(knots) > min(A))
  expect_true(max(knots) < max(A))
})

test_that("spline_basis k=0 genera indicadores correctamente", {
  x <- c(1, 2, 3, 4, 5)
  knots <- c(2, 4)
  Phi <- spline_basis(x, knots, k = 0L)
  
  expect_equal(ncol(Phi), 3) # intercepto + 2 nudos
  expect_equal(Phi[, 1], rep(1, 5))
  # Escalado: min=1, max=5. x_scaled = (x-1)/4 = 0, 0.25, 0.5, 0.75, 1
  # knots_scaled = (2-1)/4 = 0.25, (4-1)/4 = 0.75
  expect_equal(Phi[, 2], c(0, 1, 1, 1, 1)) # I(x_scaled >= 0.25)
  expect_equal(Phi[, 3], c(0, 0, 0, 1, 1)) # I(x_scaled >= 0.75)
})

test_that("spline_basis k=1 incluye termino lineal escalado (Q-01)", {
  x <- c(1, 2, 3, 4, 5)
  knots <- c(2, 4)
  Phi <- spline_basis(x, knots, k = 1L)
  
  expect_equal(ncol(Phi), 4) # intercepto + x_scaled + 2 nudos
  expect_equal(Phi[, 1], rep(1, 5))
  
  # x_scaled = (x - 1) / 4
  x_scaled <- (x - 1) / 4
  expect_equal(Phi[, 2], x_scaled) # Termino lineal puro escalado (borde 0)
  
  # (x_scaled - knots_scaled)_+
  expect_equal(Phi[, 3], pmax(0, x_scaled - 0.25)) 
  expect_equal(Phi[, 4], pmax(0, x_scaled - 0.75)) 
})
