# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

#' Pipeline completo del T-HAL-MLE nativo
#'
#' @param A Dosis observada.
#' @param W Covariable observada (escalar).
#' @param Y Resultado continuo.
#' @param grid Dosis donde se evalua la curva.
#' @param k_targeting Orden de suavidad del targeting.
#' @param k_Q Orden de suavidad del HAL del resultado (siempre explicito).
#' @param n_knots_targeting Nudos de la base del targeting.
#' @param n_knots_Q Nudos por variable en el HAL del resultado.
#' @param n_bins_g Numero de bins de W en la densidad condicional por KDE.
#' @param w_cap Tope de pesos; ver [thal_core_native()].
#' @param c1 Constante inferior de la ventana del sieve.
#' @param c2 Constante superior de la ventana del sieve.
#' @param n_proj Dosis usadas en la proyeccion final.
#' @param nlambda Valores de lambda del camino LASSO.
#' @param nfolds_Q Folds de validacion cruzada del HAL del resultado.
#' @param n_cf Folds de cross-fitting para la varianza (1 = sin CF, 2 = K-fold CF).
#' @param g_fun Opcional: funcion `function(a, w)` con la densidad condicional de A dado W para usar en lugar de la estimada
#'   (p. ej. un oraculo en simulaciones). Por defecto se estima con [fit_g_native()].
#' @param omega_fun Opcional: funcion `function(a)` con la densidad marginal de A para usar en lugar de la estimada.
#' @param seed Semilla para reproducibilidad.
#' @return Lista con `psi_hat`, `psi_plugin`, `ci` (IC de Wald; con `n_cf = 2` usa residuos fuera de fold),
#'   `ci_insample` (IC con residuos en muestra, solo si `n_cf = 2`; mismo ajuste, para comparar pareado) y diagnosticos.
#' @export
thal_fit <- function(A, W, Y, grid,
                     k_targeting = 1L, k_Q = 0L,
                     n_knots_targeting = 50L, n_knots_Q = 10L,
                     n_bins_g = 5L,
                     w_cap = NULL,
                     c1 = 6, c2 = 9,
                     n_proj = 200L, nlambda = 50L,
                     nfolds_Q = 5L,
                     n_cf = 1L,
                     g_fun = NULL, omega_fun = NULL,
                     seed = NULL) {
  if (!n_cf %in% c(1L, 2L)) stop("n_cf debe ser 1 (sin cross-fitting) o 2")
  if (!is.null(seed)) set.seed(seed)
  n <- length(A)

  if (is.null(omega_fun)) omega_fun <- fit_omega_kde_native(A)
  Q_fun     <- fit_Q_hal_native(A, W, Y, k = k_Q, n_knots = n_knots_Q,
                                nfolds = nfolds_Q, seed = seed)
  if (is.null(g_fun)) g_fun <- fit_g_native(A, W, n_bins = n_bins_g)

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
  sel <- res$selected
  eps <- res$eps

  ci_insample <- NULL
  if (length(sel) > 0) {
    # Cross-fitting K=2 para la varianza si n_cf == 2
    if (n_cf == 2L) {
      # IC con residuos en muestra, del MISMO ajuste, para comparar pareado con el de cross-fitting
      ci_insample <- wald_ci(res$psi_hat, variance_estimator(influence_curve(Y, res$ic_inputs)))
      folds <- rep(1:2, length.out = n)
      Qstar_A_cf <- numeric(n)
      
      for (k in 1:2) {
        train_idx <- which(folds != k)
        test_idx <- which(folds == k)
        
        # Reentrenar Q solo en el fold de entrenamiento
        Q_fun_k <- fit_Q_hal_native(A[train_idx], W[train_idx], Y[train_idx], 
                                    k = k_Q, n_knots = n_knots_Q, nfolds = nfolds_Q, seed = seed)
        
        # Calcular Q* para las observaciones del fold de test
        # Q*(A_i, W_i) = Q_k(A_i, W_i) + min(ω/g, cap) * φ(A_i)' ε
        Phi_test <- spline_basis(A[test_idx], knots, k_targeting, x_range = range(A))[, sel, drop = FALSE]
        wA_test <- pmin(omega_fun(A[test_idx]) / g_fun(A[test_idx], W[test_idx]), 
                        if (is.null(w_cap)) sqrt(n) * log(n/5) else if (identical(w_cap, "obs_max")) max(omega_fun(A) / g_fun(A, W)) else as.numeric(w_cap))
        Qstar_A_cf[test_idx] <- as.numeric(Q_fun_k(A[test_idx], W[test_idx]) + wA_test * (Phi_test %*% eps[sel]))
      }
      
      # Crear ic_inputs con Qstar cross-fitted
      ic_inputs_cf <- res$ic_inputs
      ic_inputs_cf$Qstar_A <- Qstar_A_cf
      ic_mat <- influence_curve(Y, ic_inputs_cf)
    } else {
      ic_mat <- influence_curve(Y, res$ic_inputs)
    }
    
    var_hat <- variance_estimator(ic_mat)
    ci <- wald_ci(res$psi_hat, var_hat)
  } else {
    ic_mat <- matrix(NA_real_, nrow = n, ncol = length(grid))
    ci <- list(lower = rep(NA_real_, length(grid)), upper = rep(NA_real_, length(grid)),
               estimate = res$psi_hat, se = rep(NA_real_, length(grid)))
  }

  list(
    psi_hat    = res$psi_hat,
    psi_plugin = res$psi_plugin,
    ci         = ci,
    ci_insample = ci_insample,
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
