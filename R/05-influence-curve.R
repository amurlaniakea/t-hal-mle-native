# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

#' Curva de influencia del T-HAL-MLE (Lemma 3.1 del paper, version empirica)
#'
#' Para la base seleccionada phi (s columnas) y la medida de A del paso 7
#' (los `m` puntos de proyeccion), con Sigma = P[phi phi'] y psi*_pl(a) = mean_i Q*(a, W_i):
#'
#'   D_i(a) = phi(a)' Sigma^+ [ r_i + s_i ]
#'   r_i = min(omega/g, cap)(A_i) * phi(A_i) * (Y_i - Q*(A_i, W_i))          (componente Q-barra)
#'   s_i = (1/m) sum_j phi(a_j) * ( Q*(a_j, W_i) - psi*_pl(a_j) )            (componente Q_W)
#'
#' Q* es el MISMO Q dirigido que usa el estimador puntual (Q + omega/g * phi' eps), no una
#' version sin el peso. La media empirica de s_i es exactamente 0; la de r_i es solo
#' aproximadamente 0 porque el LASSO resuelve las ecuaciones de score de forma aproximada.
#'
#' @param Y Resultado observado (n).
#' @param ic_inputs Lista `ic_inputs` devuelta por [thal_core_native()].
#' @param tol Tolerancia relativa de la pseudo-inversa de Sigma.
#' @return Matriz n x length(grid) con D_i(a_g).
#' @export
influence_curve <- function(Y, ic_inputs, tol = 1e-8) {
  if (is.null(ic_inputs)) stop("ic_inputs es NULL: el targeting no selecciono funciones base.")
  with(ic_inputs, {
    n <- length(Y)
    m <- nrow(Phi_proj)
    Sigma <- crossprod(Phi_proj) / m
    Sigma_pinv <- pinv_sym(Sigma, tol = tol)

    psi_pl <- colMeans(Qstar_proj)                              # psi*_pl(a_j)
    R <- (wA * (Y - Qstar_A)) * Phi_A                           # n x s
    S <- (Qstar_proj - matrix(psi_pl, n, m, byrow = TRUE)) %*% Phi_proj / m  # n x s

    ((R + S) %*% Sigma_pinv) %*% t(Phi_grid)                    # n x G
  })
}

#' Componentes de la curva de influencia (Q-barra y Q_W) proyectadas a la rejilla
#'
#' Descompone `D_i(a) = R_i(a) + S_i(a)`: `R` es la componente de residuos (Q-barra) y `S` la de la variacion de W
#' (Q_W), ambas ya multiplicadas por `phi(a)' Sigma^+`. Sirve para diagnosticar cual de las dos domina la varianza
#' estimada en cada dosis.
#'
#' @inheritParams influence_curve
#' @return Lista con las matrices `R` y `S` (n x length(grid)); `R + S` coincide con [influence_curve()].
#' @export
influence_parts <- function(Y, ic_inputs, tol = 1e-8) {
  if (is.null(ic_inputs)) stop("ic_inputs es NULL: el targeting no selecciono funciones base.")
  with(ic_inputs, {
    n <- length(Y)
    m <- nrow(Phi_proj)
    Sigma_pinv <- pinv_sym(crossprod(Phi_proj) / m, tol = tol)
    psi_pl <- colMeans(Qstar_proj)
    R <- (wA * (Y - Qstar_A)) * Phi_A
    S <- (Qstar_proj - matrix(psi_pl, n, m, byrow = TRUE)) %*% Phi_proj / m
    list(R = (R %*% Sigma_pinv) %*% t(Phi_grid), S = (S %*% Sigma_pinv) %*% t(Phi_grid))
  })
}

#' Varianza de psi_hat(a): P[D(a)^2] / n
#'
#' @param ic_matrix Matriz n x G de la curva de influencia, de [influence_curve()].
#' @export
variance_estimator <- function(ic_matrix) {
  n <- nrow(ic_matrix)
  colMeans(ic_matrix^2) / n
}

#' Intervalos de Wald puntuales
#'
#' @param psi_hat Estimaciones puntuales.
#' @param var_hat Varianzas estimadas.
#' @param alpha Nivel de significacion (0.05 da IC al 95\%).
#' @export
wald_ci <- function(psi_hat, var_hat, alpha = 0.05) {
  z <- stats::qnorm(1 - alpha / 2)
  se <- sqrt(var_hat)
  list(lower = psi_hat - z * se, upper = psi_hat + z * se, estimate = psi_hat, se = se)
}
