test_that("fit_omega_kde_native devuelve una densidad valida", {
  set.seed(1)
  A <- runif(500, 1, 9)
  omega_fun <- fit_omega_kde_native(A)

  grid <- seq(1, 9, length.out = 200)
  vals <- omega_fun(grid)

  expect_true(all(vals > 0))
  # Integral sobre el rango de datos; tolerancia amplia por sesgo de borde de la KDE
  integral <- sum(vals * diff(grid)[1])
  expect_lt(abs(integral - 1), 0.15)
})

test_that("fit_omega_kde_native aproxima la uniforme (1/8) en el interior", {
  set.seed(2)
  A <- runif(500, 1, 9)
  omega_fun <- fit_omega_kde_native(A)

  # Solo el interior: la KDE tiene sesgo de borde cerca de 1 y 9
  grid <- seq(2, 8, length.out = 100)
  vals <- omega_fun(grid)
  truth <- rep(1 / 8, length(grid))

  rmse <- sqrt(mean((vals - truth)^2))
  expect_lt(rmse, 0.02)
})

test_that("fit_omega_kde_native aproxima una normal truncada a [1,9]", {
  set.seed(3)
  # Normal truncada limpia por rechazo (sin masas puntuales en los bordes)
  A <- numeric(0)
  while (length(A) < 500) {
    cand <- rnorm(600, mean = 5, sd = 2)
    A <- c(A, cand[cand >= 1 & cand <= 9])
  }
  A <- A[seq_len(500)]

  omega_fun <- fit_omega_kde_native(A)

  grid <- seq(1.5, 8.5, length.out = 100)
  vals <- omega_fun(grid)
  # Verdad: N(5, 2^2) truncada a [1,9]
  truth <- dnorm((grid - 5) / 2) / (2 * (pnorm(2) - pnorm(-2)))

  expect_gt(cor(vals, truth), 0.95)
})
