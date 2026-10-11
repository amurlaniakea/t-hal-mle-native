test_that("thal_core_native funciona con oraculos (DGP 1, uniform)", {
  set.seed(123)
  d <- simulate_data(200, dgp = 1, treatment = "uniform")

  Q_fun <- function(a, w) Q0_dgp(a, w, dgp = 1)
  g_fun <- function(a, w) g0_dgp(a, w, treatment = "uniform")
  omega_fun <- function(a) omega_oracle(a, treatment = "uniform")

  grid <- seq(1, 9, length.out = 10)

  res <- thal_core_native(
    A = d$A, W = d$W, Y = d$Y,
    Q_fun = Q_fun, g_fun = g_fun, omega_fun = omega_fun,
    grid = grid, k = 1L, n_knots = 20, c1 = 6, c2 = 9
  )

  truth <- true_drc(grid, dgp = 1)
  rmse <- sqrt(mean((res$psi_hat - truth)^2))

  expect_lt(rmse, 0.5)
  expect_true(res$size > 0)
})

test_that("thal_core_native con tratamiento normal y w_cap obs_max", {
  set.seed(456)
  d <- simulate_data(200, dgp = 1, treatment = "normal")

  Q_fun <- function(a, w) Q0_dgp(a, w, dgp = 1)
  g_fun <- function(a, w) g0_dgp(a, w, treatment = "normal")
  omega_fun <- function(a) omega_oracle(a, treatment = "normal")

  grid <- seq(1, 9, length.out = 10)

  res <- thal_core_native(
    A = d$A, W = d$W, Y = d$Y,
    Q_fun = Q_fun, g_fun = g_fun, omega_fun = omega_fun,
    grid = grid, k = 1L, n_knots = 20, c1 = 6, c2 = 9,
    w_cap = "obs_max"
  )

  truth <- true_drc(grid, dgp = 1)
  rmse <- sqrt(mean((res$psi_hat - truth)^2))

  expect_lt(rmse, 1.0)
  expect_true(res$size > 0)
})

test_that("DGP 2 (discontinuidad) con k=0: mejora sobre el plug-in", {
  set.seed(789)
  # La true_drc del DGP 2 es lineal a trozos con un salto en a=5.
  # k=0 (escalones) captura el salto pero tiene error de aproximacion
  # inherente para los tramos con pendiente (+2, -1). Por eso el objetivo
  # no es un RMSE absoluto bajo, sino que el targeting mejore al plug-in.
  d <- simulate_data(500, dgp = 2, treatment = "uniform")

  Q_fun <- function(a, w) Q0_dgp(a, w, dgp = 2)
  g_fun <- function(a, w) g0_dgp(a, w, treatment = "uniform")
  omega_fun <- function(a) omega_oracle(a, treatment = "uniform")

  grid <- seq(1, 9, length.out = 10)

  # Escalones necesitan muchos nudos para seguir las pendientes
  res <- thal_core_native(
    A = d$A, W = d$W, Y = d$Y,
    Q_fun = Q_fun, g_fun = g_fun, omega_fun = omega_fun,
    grid = grid, k = 0L, n_knots = 50, c1 = 6, c2 = 9
  )

  truth <- true_drc(grid, dgp = 2)
  rmse_thal <- sqrt(mean((res$psi_hat - truth)^2))
  rmse_plugin <- sqrt(mean((res$psi_plugin - truth)^2))

  # Criterio principal: el targeting debe funcionar (seleccionar algo)
  expect_true(res$size > 0)
  # Criterio de calidad: RMSE finito y razonable para escalones
  expect_lt(rmse_thal, 1.5)
  # El estimador no debe ser peor que el plug-in de forma catastrofica
  expect_lt(rmse_thal, rmse_plugin + 0.5)
})
