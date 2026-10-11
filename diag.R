
pkgload::load_all(".", quiet = TRUE)

set.seed(123)
d <- simulate_data(500, dgp = 1, treatment = "uniform")

Q_fun <- function(a, w) Q0_dgp(a, w, dgp = 1)
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

cat("=== DIAGNOSTICO DGP 1, uniform, oraculos ===\n")
cat(sprintf("RMSE plug-in (sin targeting): %.3f\n", rmse_plugin))
cat(sprintf("RMSE targeting: %.3f\n", rmse_thal))
cat(sprintf("Tamaño modelo seleccionado: %d\n", res$size))
cat(sprintf("Lambda seleccionado: %.4f\n", res$lambda))
cat("\nPrimeros 5 puntos:\n")
cat(sprintf("  a=%.2f | truth=%.3f | plugin=%.3f | thal=%.3f\n",
            grid[1:5], truth[1:5], res$psi_plugin[1:5], res$psi_hat[1:5]))
