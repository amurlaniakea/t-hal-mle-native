# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

#' Nucleo del T-HAL-MLE (Algorithm 1, pasos 3-8)
#'
#' @param A Dosis observada.
#' @param W Covariable observada.
#' @param Y Resultado continuo.
#' @param Q_fun Funcion (a, w) -> estimacion inicial de E[Y|A=a,W=w].
#' @param g_fun Funcion (a, w) -> estimacion de la densidad condicional g(a|w).
#' @param omega_fun Funcion a -> densidad marginal de A.
#' @param grid Dosis donde se evalua la curva.
#' @param k Orden de suavidad del objetivo (la base del targeting).
#' @param n_knots Numero de nudos de la base del targeting.
#' @param c1 Constante inferior de la ventana del tamano del sieve.
#' @param c2 Constante superior de la ventana del tamano del sieve.
#' @param n_proj Numero maximo de dosis usadas en la proyeccion final (paso 7).
#' @param nlambda Numero de valores de lambda del camino LASSO.
#' @param w_cap Tope de los pesos omega/g: NULL (sqrt(n) log(n/5), como el paper), un numero, o "obs_max" (maximo observado; heuristica no incluida en el paper).
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
  
  # Matriz n x m de Q*(a_j, W_i) = Q(a_j, W_i) + min(omega(a_j)/g(a_j|W_i), cap) * phi(a_j)' eps
  marginal_mat <- function(points) {
    m <- length(points)
    a_rep <- rep(points, each = n)
    w_rep <- rep(W, times = m)
    Qb <- matrix(Q_fun(a_rep, w_rep), n, m)
    gb <- matrix(g_fun(a_rep, w_rep), n, m)
    # CRITICO: usar x_range_A para escalar consistentemente
    corr <- as.numeric(spline_basis(points, knots, k, x_range = x_range_A) %*% eps)
    Wmat <- pmin(sweep(1 / gb, 2, omega_fun(points), "*"), cap)
    Qb + sweep(Wmat, 2, corr, "*")
  }
  marginal <- function(points) colMeans(marginal_mat(points))

  psi_plugin <- marginal(grid)
  sel <- which(eps != 0)
  
  if (length(sel) == 0L) {
    warning("El paso de targeting no selecciono ninguna funcion base; devolviendo plug-in.")
    return(list(psi_hat = psi_plugin, psi_plugin = psi_plugin, grid = grid,
                eps = eps, selected = sel, size = 0L, lambda = fit$lambda[idx],
                ic_inputs = NULL))
  }
  
  proj_pts <- if (n > n_proj) {
    unique(quantile(A, probs = seq(0, 1, length.out = n_proj), names = FALSE))
  } else sort(A)

  # CRITICO: usar x_range_A para escalar consistentemente
  Phi_p <- spline_basis(proj_pts, knots, k, x_range = x_range_A)[, sel, drop = FALSE]
  Qstar_proj <- marginal_mat(proj_pts)
  alpha <- pinv_sym(crossprod(Phi_p)) %*% crossprod(Phi_p, colMeans(Qstar_proj))
  Phi_g <- spline_basis(grid, knots, k, x_range = x_range_A)[, sel, drop = FALSE]
  psi_hat <- as.numeric(Phi_g %*% alpha)

  # Insumos de la curva de influencia (Lemma 3.1), todos coherentes con el ajuste:
  #  - Qstar_A: valor ajustado del LASSO en los datos, Q + min(omega/g, cap) * phi' eps
  #  - Qstar_proj: Q*(a_j, W_i) en los puntos de proyeccion (medida de A del paso 7)
  Qstar_A <- as.numeric(offset + H %*% eps)
  ic_inputs <- list(
    wA = wA, Qstar_A = Qstar_A,
    Phi_A = Phi_A[, sel, drop = FALSE],
    Phi_proj = Phi_p, Qstar_proj = Qstar_proj,
    Phi_grid = Phi_g, proj_pts = proj_pts
  )

  list(psi_hat = psi_hat, psi_plugin = psi_plugin, grid = grid, eps = eps,
       selected = sel, size = length(sel), lambda = fit$lambda[idx],
       ic_inputs = ic_inputs)
}
