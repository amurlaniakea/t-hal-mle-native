# Copyright (C) 2026 Pedro Sordo Martinez
# SPDX-License-Identifier: AGPL-3.0-or-later

#' Pipeline completo del T-HAL-MLE nativo
#'
#' Ajusta los tres nuisances (omega_n, Q_n, g_n) con los estimadores nativos
#' y llama a thal_core_native. Calcula la curva de influencia y los IC de Wald.
#'
#' @param A Vector de dosis observadas.
#' @param W Vector de covariable observada.
#' @param Y Vector de outcome.
#' @param grid Puntos de la grid donde estimar psi(a).
#' @param k_targeting Orden de la base del targeting (default 1).
#' @param k_Q Orden de la base del outcome (default 0, como el paper para DGP 1/2).
#' @param n_knots_targeting Numero de nudos de la base del targeting.
#' @param n_knots_Q Numero de nudos por dimension para el outcome.
#' @param n_bins_g Numero de bins para la densidad condicional.
#' @param w_cap Cap para pesos omega/g. NULL = sqrt(n)*log(n/5); "obs_max" = max observado; numerico = fijo.
#' @param c1, c2 Constantes del selector de lambda (Seccion 5.1 del paper).
#' @param n_proj Numero maximo de puntos de dosis para la proyeccion OLS.
#' @param nlambda Numero de lambdas en el camino LASSO del targeting.
#' @param nfolds_Q Folds para el CV del outcome.
#' @param seed Semilla para los CV.
#' @return Lista con psi_hat, psi_plugin, ci (IC de Wald), ic_matrix, grid,
#'   eps (coeficientes del targeting), size, lambda, knots, selected.
#' @export
thal_fit <- function(A, W, Y, grid,
                     k_targeting = 1L, k_Q = 0L,
                     n_knots_targeting = 50L, n_knots_Q = 10L,
                     n_bins_g = 5L,
                     w_cap = "obs_max",
                     c1 = 6, c2 = 9,
                     n_proj = 200L, nlambda = 50L,
                     nfolds_Q = 5L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # 1. Ajustar nuisances nativos
  omega_fun <- fit_omega_kde_native(A)
  Q_fun     <- fit_Q_hal_native(A, W, Y, k = k_Q, n_knots = n_knots_Q,
                                nfolds = nfolds_Q, seed = seed)
  g_fun     <- fit_g_native(A, W, n_bins = n_bins_g)

  # 2. Core del targeting
  res <- thal_core_native(
    A = A, W = W, Y = Y,
    Q_fun = Q_fun, g_fun = g_fun, omega_fun = omega_fun,
    grid = grid,
    k = k_targeting, n_knots = n_knots_targeting,
    c1 = c1, c2 = c2,
    n_proj = n_proj, nlambda = nlambda,
    w_cap = w_cap
  )

  # 3. Curva de influencia y IC de Wald (D-0004)
  # Construir Q*(a, w) = Q(a, w) + phi_sel(a) %*% eps
  knots <- make_knots(A, n_knots_targeting)
  sel <- res$selected
  eps <- res$eps

  if (length(sel) > 0) {
    Q_star <- function(a, w) {
      Phi_a <- spline_basis(a, knots, k_targeting)[, sel, drop = FALSE]
      corr <- as.numeric(Phi_a %*% eps[sel])
      Q_fun(a, w) + corr
    }
  } else {
    Q_star <- Q_fun
  }

  # Calcular IC
  ic_mat <- influence_curve(
    A = A, W = W, Y = Y,
    Q_star = Q_star,
    g_fun = g_fun, omega_fun = omega_fun,
    knots = knots, k = k_targeting,
    selected = sel, alpha = eps[sel],
    psi_hat = res$psi_hat, grid = grid
  )

  var_hat <- variance_estimator(ic_mat)
  ci <- wald_ci(res$psi_hat, var_hat)

  # 4. Devolver resultado
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
    # Metadata del ajuste
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
