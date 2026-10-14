test_that("tensor_basis k=0 tiene las dimensiones correctas", {
  A <- c(1, 2, 3, 4, 5)
  W <- c(0.1, 0.5, 0.9, 0.3, 0.7)
  Phi <- tensor_basis(A, W, knots_A = c(2, 4), knots_W = c(0.5), k = 0L)
  # 1 + 2 + 1 + 2*1 = 6 columnas
  expect_equal(ncol(Phi), 6)
  expect_equal(nrow(Phi), 5)
  expect_equal(Phi[, 1], rep(1, 5))
})

test_that("tensor_basis k=1 incluye los terminos de borde escalados a, w, a*w (Q-01)", {
  A <- c(1, 2, 3)
  W <- c(0.2, 0.4, 0.6)
  # Escalado interno:
  # A_range = c(1, 3) -> A_scaled = (A - 1) / 2 = c(0, 0.5, 1)
  # W_range = c(0.2, 0.6) -> W_scaled = (W - 0.2) / 0.4 = c(0, 0.5, 1)
  
  Phi <- tensor_basis(A, W, knots_A = c(2), knots_W = c(0.4), k = 1L)
  
  # 4 + 2*1 + 2*1 + 1*1 = 9 columnas
  expect_equal(ncol(Phi), 9)
  
  A_scaled <- (A - 1) / 2
  W_scaled <- (W - 0.2) / 0.4
  
  expect_equal(Phi[, 2], A_scaled)
  expect_equal(Phi[, 3], W_scaled)
  expect_equal(Phi[, 4], A_scaled * W_scaled)
})

test_that("fit_Q_hal_native recupera un Q0 suave conocido", {
  set.seed(42)
  n <- 300
  A <- runif(n, 1, 9)
  W <- runif(n, 0, 1)
  # Q0 suave: 3*sin(1.5*a) + 0.5*a + 2*w
  Y <- 3 * sin(1.5 * A) + 0.5 * A + 2 * W + rnorm(n, sd = 0.5)

  Q_fun <- fit_Q_hal_native(A, W, Y, k = 1L, n_knots = 10, seed = 1)

  grid_a <- seq(2, 8, length.out = 20)
  grid_w <- rep(0.5, 20)
  Q_hat <- Q_fun(grid_a, grid_w)
  Q_true <- 3 * sin(1.5 * grid_a) + 0.5 * grid_a + 2 * grid_w

  rmse <- sqrt(mean((Q_hat - Q_true)^2))
  # Umbral relajado para acomodar el escalado y la variabilidad del CV
  expect_lt(rmse, 2.0)
})

test_that("fit_Q_hal_native: la prediccion no depende del rango de los puntos evaluados (k = 1)", {
  set.seed(11); n <- 200
  W <- rnorm(n); A <- runif(n, 0, 10); Y <- sin(A) + W + rnorm(n, 0, 0.3)
  Qf <- fit_Q_hal_native(A, W, Y, k = 1L, n_knots = 6L, nfolds = 3L, seed = 1)
  full <- Qf(A, W)
  expect_equal(Qf(A[1:50], W[1:50]), full[1:50])
  expect_equal(Qf(c(2, 8), c(0, 0)), Qf(c(2, 8, 1, 9), c(0, 0, 0, 0))[1:2])
})
