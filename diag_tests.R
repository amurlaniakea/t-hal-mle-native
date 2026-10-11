
pkgload::load_all(".", quiet = TRUE)

cat("=== TEST 1: DGP 1, uniform, Q misspecified ===\n")
set.seed(123)
d <- simulate_data(500, dgp = 1, treatment = "uniform")
Q_fun <- function(a, w) 3 * sin(1.5 * a) + 0.5 * a
g_fun <- function(a, w) g0_dgp(a, w, treatment = "uniform")
omega_fun <- function(a) omega_oracle(a, treatment = "uniform")
grid <- seq(2, 8, length.out = 10)

res <- thal_core_native(
  A = d$A, W = d$W, Y = d$Y,
  Q_fun = Q_fun, g_fun = g_fun, omega_fun = omega_fun,
  grid = grid, k = 1L, n_knots = 30, c1 = 6, c2 = 9
)

truth <- true_drc(grid, dgp = 1)
rmse_plugin <- sqrt(mean((res$psi_plugin - truth)^2))
rmse_thal <- sqrt(mean((res$psi_hat - truth)^2))

cat(sprintf("RMSE plugin: %.3f\n", rmse_plugin))
cat(sprintf("RMSE thal: %.3f\n", rmse_thal))
cat(sprintf("is.finite(rmse_thal): %s\n", is.finite(rmse_thal)))
cat(sprintf("size: %d\n", res$size))
cat(sprintf("psi_hat range: [%.3f, %.3f]\n", min(res$psi_hat), max(res$psi_hat)))
cat(sprintf("truth range: [%.3f, %.3f]\n", min(truth), max(truth)))

cat("\n=== TEST 2: DGP 1, normal, Q misspecified, w_cap obs_max ===\n")
set.seed(456)
d <- simulate_data(500, dgp = 1, treatment = "normal")
Q_fun <- function(a, w) 3 * sin(1.5 * a) + 0.5 * a
g_fun <- function(a, w) g0_dgp(a, w, treatment = "normal")
omega_fun <- function(a) omega_oracle(a, treatment = "normal")

res <- thal_core_native(
  A = d$A, W = d$W, Y = d$Y,
  Q_fun = Q_fun, g_fun = g_fun, omega_fun = omega_fun,
  grid = grid, k = 1L, n_knots = 30, c1 = 6, c2 = 9,
  w_cap = "obs_max"
)

truth <- true_drc(grid, dgp = 1)
rmse_thal <- sqrt(mean((res$psi_hat - truth)^2))

cat(sprintf("RMSE thal: %.3f\n", rmse_thal))
cat(sprintf("is.finite(rmse_thal): %s\n", is.finite(rmse_thal)))
cat(sprintf("size: %d\n", res$size))

cat("\n=== TEST 3: Pipeline DGP 1, uniform ===\n")
set.seed(100)
d <- simulate_data(300, dgp = 1, treatment = "uniform")
grid <- seq(2, 8, length.out = 10)

res <- thal_fit(
  A = d$A, W = d$W, Y = d$Y, grid = grid,
  k_targeting = 1L, k_Q = 1L,
  n_knots_targeting = 20, n_knots_Q = 8,
  n_bins_g = 5,
  w_cap = "obs_max",
  seed = 42
)

truth <- true_drc(grid, dgp = 1)
rmse <- sqrt(mean((res$psi_hat - truth)^2))

cat(sprintf("RMSE: %.3f\n", rmse))
cat(sprintf("is.finite(rmse): %s\n", is.finite(rmse)))
cat(sprintf("psi_hat range: [%.3f, %.3f]\n", min(res$psi_hat), max(res$psi_hat)))
cat(sprintf("truth range: [%.3f, %.3f]\n", min(truth), max(truth)))
