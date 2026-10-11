test_that("componente Q_W de la curva de influencia tiene media empirica exactamente 0", {
  set.seed(1)
  n <- 60; m <- 12; s <- 3
  I <- list(
    wA = rep(0, n), Qstar_A = rnorm(n),
    Phi_A = matrix(rnorm(n * s), n, s),
    Phi_proj = matrix(rnorm(m * s), m, s),
    Qstar_proj = matrix(rnorm(n * m), n, m),
    Phi_grid = matrix(rnorm(4 * s), 4, s)
  )
  ic <- influence_curve(rnorm(n), I)          # wA = 0: solo queda la componente Q_W
  expect_equal(dim(ic), c(n, 4))
  expect_lt(max(abs(colMeans(ic))), 1e-10)
})

test_that("ic_inputs es coherente con el estimador puntual (mismo Q*)", {
  set.seed(11)
  n <- 300
  d <- simulate_data(n, dgp = 1, treatment = "normal")
  Q_fun <- function(a, w) 3 * sin(1.5 * a) + 0.5 * a
  g_fun <- function(a, w) g0_dgp(a, w, "normal")
  om <- function(a) omega_oracle(a, "normal")
  grid <- seq(2, 8, length.out = 5)
  r <- thal_core_native(d$A, d$W, d$Y, Q_fun, g_fun, om, grid, k = 1L, n_knots = 50)
  I <- r$ic_inputs
  expect_false(is.null(I))

  # Q* en los datos = valor ajustado del LASSO (Q + min(omega/g, cap) * phi' eps)
  cap <- sqrt(n) * log(n / 5)
  Phi <- spline_basis(d$A, make_knots(d$A, 50), 1L, x_range = range(d$A))
  fitted <- Q_fun(d$A, d$W) + pmin(om(d$A) / g_fun(d$A, d$W), cap) * as.numeric(Phi %*% r$eps)
  expect_equal(I$Qstar_A, fitted, tolerance = 1e-8)

  # psi_hat se reconstruye exactamente a partir de los insumos de la curva de influencia
  m <- nrow(I$Phi_proj)
  alpha <- pinv_sym(crossprod(I$Phi_proj)) %*% crossprod(I$Phi_proj, colMeans(I$Qstar_proj))
  expect_equal(as.numeric(I$Phi_grid %*% alpha), r$psi_hat, tolerance = 1e-8)
})

test_that("con nuisances oraculo la media del IC es pequena frente a su desviacion", {
  set.seed(5)
  n <- 300
  d <- simulate_data(n, dgp = 1, treatment = "uniform")
  Q_fun <- function(a, w) 3 * sin(1.5 * a) + 0.5 * a
  g_fun <- function(a, w) g0_dgp(a, w, "uniform")
  om <- function(a) omega_oracle(a, "uniform")
  r <- thal_core_native(d$A, d$W, d$Y, Q_fun, g_fun, om, seq(2, 8, length.out = 5),
                        k = 1L, n_knots = 50)
  ic <- influence_curve(d$Y, r$ic_inputs)
  tstat <- colMeans(ic) * sqrt(n) / apply(ic, 2, sd)   # ~N(0,1) si el score esta resuelto
  expect_true(all(abs(tstat) < 3))
  se <- sqrt(variance_estimator(ic))
  ci <- wald_ci(r$psi_hat, se^2)
  expect_true(all(ci$lower < r$psi_hat & r$psi_hat < ci$upper))
})

test_that("select_lambda_by_size usa el exponente 1/(2k*+1) segun k", {
  df <- c(0, 10, 25, 60)
  # n = 500, k = 1: ventana [21, 31] -> indice 3; k = 0: ventana [48, 71] -> indice 4
  expect_equal(select_lambda_by_size(df, 500, k = 1L), 3L)
  expect_equal(select_lambda_by_size(df, 500, k = 0L), 4L)
})
