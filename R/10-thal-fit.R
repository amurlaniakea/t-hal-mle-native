# Copyright (C) 2026 Pedro Sordo Martinez
# SPDX-License-Identifier: AGPL-3.0-or-later

#' Pipeline completo del T-HAL-MLE nativo
#' @export
thal_fit <- function(A, W, Y, grid,
                     k_targeting = 1L, k_Q = 0L,
                     n_knots_targeting = 50L, n_knots_Q = 10L,
                     n_bins_g = 5L,
                     w_cap = NULL,  # Default del paper: sqrt(n)*log(n/5)
                     c1 = 6, c2 = 9,
                     n_proj = 200L, nlambda = 50L,
                     nfolds_Q = 5L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  omega_fun <- fit_omega_kde_native(A)
  Q_fun     <- fit_Q_hal_native(A, W, Y, k = k_Q, n_knots = n_knots_Q,
                                nfolds = nfolds_Q, seed = seed)
  g_fun     <- fit_g_native(A, W, n_bins = n_bins_g)

  res <- thal_core_native(
    A = A, W = W, Y = Y,
    Q_fun = Q_fun, g_fun = g_fun, omega_fun = omega_fun,
    grid = grid,
    k = k_targeting, n_knots = n_knots_targeting,
    c1 = c1, c2 = c2,
    n_proj = n_proj, nlambda = nlambda,
    w_cap = w_cap
  )

  knots <- make_knots(A, n_knots_targeting)
  x_range_A <- range(A)
  sel <- res$selected
  eps <- res$eps

  if (length(sel) > 0) {
    Q_star <- function(a, w) {
      Phi_a <- spline_basis(a, knots, k_targeting, x_range = x_range_A)[, sel, drop = FALSE]
      corr <- as.numeric(Phi_a %*% eps[sel])
      Q_fun(a, w) + corr
    }
  } else {
    Q_star <- Q_fun
  }

  ic_mat <- influence_curve(
    A = A, W = W, Y = Y,
    Q_star = Q_star,
    g_fun = g_fun, omega_fun = omega_fun,
    knots = knots, k = k_targeting,
    selected = sel, alpha = eps[sel],
    psi_hat = res$psi_hat, grid = grid,
    x_range = x_range_A
  )

  var_hat <- variance_estimator(ic_mat)
  ci <- wald_ci(res$psi_hat, var_hat)

  list(
    psi_hat    = res$psi_hat,
    psi_plugin = res$psi_plugin,
    ci         = ci,
    ic_matrix  = ic_mat,
    grid       = grid,
    eps        = eps,
    selected   = sel,
    size       = res$size,
    lambda     = res$lambda,
    knots      = knots,
    nuisance   = list(
      omega_fun = omega_fun, Q_fun = Q_fun, g_fun = g_fun,
      k_Q = k_Q, n_knots_Q = n_knots_Q, n_bins_g = n_bins_g
    ),
    targeting  = list(
      k = k_targeting, n_knots = n_knots_targeting,
      c1 = c1, c2 = c2, w_cap = w_cap
    )
  )
}
