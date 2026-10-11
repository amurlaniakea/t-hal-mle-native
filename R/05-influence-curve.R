# Copyright (C) 2026 Pedro Sordo Martinez
# SPDX-License-Identifier: AGPL-3.0-or-later

#' Calcula la curva de influencia del estimador T-HAL-MLE (Lemma 3.1 corregido)
#' @export
influence_curve <- function(A, W, Y, Q_star, g_fun, omega_fun, 
                            knots, k, selected, alpha, psi_hat, grid,
                            x_range = range(A)) {
  n <- length(A)
  m <- length(grid)
  
  Phi_A_sel <- spline_basis(A, knots, k, x_range = x_range)[, selected, drop = FALSE]
  weights_A <- omega_fun(A) / g_fun(A, W)
  residuals_A <- Y - Q_star(A, W)
  component_Q_bar <- weights_A * residuals_A * Phi_A_sel 
  
  component_Q_W <- matrix(0, nrow = n, ncol = length(selected))
  for (j in seq_along(grid)) {
    a_j <- grid[j]
    Phi_a_j <- spline_basis(a_j, knots, k, x_range = x_range)[, selected, drop = FALSE]
    Q_star_a_j_W <- Q_star(rep(a_j, n), W)
    diff_Q <- Q_star_a_j_W - psi_hat[j]
    component_Q_W <- component_Q_W + (1 / m) * (diff_Q %*% Phi_a_j)
  }
  
  ic_matrix <- matrix(0, nrow = n, ncol = m)
  for (j in seq_along(grid)) {
    a <- grid[j]
    Phi_a <- spline_basis(a, knots, k, x_range = x_range)[, selected, drop = FALSE]
    ic_matrix[, j] <- as.numeric(Phi_a %*% t(component_Q_bar + component_Q_W))
  }
  ic_matrix
}

#' @export
variance_estimator <- function(ic_matrix) {
  n <- nrow(ic_matrix)
  colMeans(ic_matrix^2) / n
}

#' @export
wald_ci <- function(psi_hat, var_hat, alpha = 0.05) {
  z <- qnorm(1 - alpha / 2)
  se <- sqrt(var_hat)
  list(lower = psi_hat - z * se, upper = psi_hat + z * se, estimate = psi_hat, se = se)
}
