# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

#' Nucleo del T-HAL-MLE (Algorithm 1, pasos 3-8)
#' @export
thal_core_native <- function(A, W, Y, Q_fun, g_fun, omega_fun, grid,
                             k = 1L, n_knots = 50L, c1 = 6, c2 = 9,
                             n_proj = 200L, nlambda = 50L, w_cap = NULL) {
  n <- length(A)
  
  if (is.null(w_cap)) {
    cap <- sqrt(n) * log(n / 5)
  } else if (identical(w_cap, "obs_max")) {
    cap <- max(omega_fun(A) / g_fun(A, W))
  } else {
    cap <- as.numeric(w_cap)
  }
  
  knots <- make_knots(A, n_knots)
  x_range_A <- range(A)
  Phi_A <- spline_basis(A, knots, k, x_range = x_range_A)
  
  wA <- pmin(omega_fun(A) / g_fun(A, W), cap)
  H <- wA * Phi_A
  
  offset <- Q_fun(A, W)
  
  fit <- lasso_path(H, Y, offset = offset, nlambda = nlambda)
  
  idx <- select_lambda_by_size(fit$df, n, k, c1, c2)
  eps <- fit$beta[, idx]
  
  marginal <- function(points) {
    m <- length(points)
    a_rep <- rep(points, each = n)
    w_rep <- rep(W, times = m)
    Qb <- matrix(Q_fun(a_rep, w_rep), n, m)
    gb <- matrix(g_fun(a_rep, w_rep), n, m)
    # CRITICO: usar x_range_A para escalar consistentemente
    corr <- as.numeric(spline_basis(points, knots, k, x_range = x_range_A) %*% eps)
    Wmat <- pmin(sweep(1 / gb, 2, omega_fun(points), "*"), cap)
    Qstar <- Qb + sweep(Wmat, 2, corr, "*")
    colMeans(Qstar)
  }
  
  psi_plugin <- marginal(grid)
  sel <- which(eps != 0)
  
  if (length(sel) == 0L) {
    warning("El paso de targeting no selecciono ninguna funcion base; devolviendo plug-in.")
    return(list(psi_hat = psi_plugin, psi_plugin = psi_plugin, grid = grid,
                eps = eps, selected = sel, size = 0L, lambda = fit$lambda[idx]))
  }
  
  proj_pts <- if (n > n_proj) {
    unique(quantile(A, probs = seq(0, 1, length.out = n_proj), names = FALSE))
  } else sort(A)
  
  # CRITICO: usar x_range_A para escalar consistentemente
  Phi_p <- spline_basis(proj_pts, knots, k, x_range = x_range_A)[, sel, drop = FALSE]
  alpha <- pinv_sym(crossprod(Phi_p)) %*% crossprod(Phi_p, marginal(proj_pts))
  psi_hat <- as.numeric(spline_basis(grid, knots, k, x_range = x_range_A)[, sel, drop = FALSE] %*% alpha)
  
  list(psi_hat = psi_hat, psi_plugin = psi_plugin, grid = grid, eps = eps,
       selected = sel, size = length(sel), lambda = fit$lambda[idx])
}
