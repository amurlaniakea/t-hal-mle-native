# Copyright (C) 2026 Pedro Sordo Martinez
# SPDX-License-Identifier: AGPL-3.0-or-later

#' Calcula la curva de influencia del estimador T-HAL-MLE
#' 
#' Segun Lemma 3.1 del paper (ec. 1):
#' D*(a) = (omega(a)/g(a|W)) * (Y - Q*(a,W)) * phi_sel(a) * alpha - psi_hat(a)
#' 
#' @param A Vector de dosis observadas
#' @param W Vector de covariables
#' @param Y Vector de outcomes
#' @param Q_star Funcion Q* ajustada (post-targeting)
#' @param g_fun Funcion de densidad condicional
#' @param omega_fun Funcion de densidad marginal
#' @param knots Nudos de la base
#' @param k Orden de la base
#' @param selected Indices de funciones seleccionadas
#' @param alpha Coeficientes de la proyeccion
#' @param psi_hat Vector de estimaciones en la grid
#' @param grid Puntos de evaluacion
#' @export
influence_curve <- function(A, W, Y, Q_star, g_fun, omega_fun, 
                            knots, k, selected, alpha, psi_hat, grid) {
  n <- length(A)
  m <- length(grid)
  
  # Para cada punto de la grid, calcular D*(a)
  ic_matrix <- matrix(0, nrow = n, ncol = m)
  
  for (j in seq_along(grid)) {
    a <- grid[j]
    
    # Evaluar phi_sel(a) para las funciones seleccionadas
    Phi_a <- spline_basis(a, knots, k)[, selected, drop = FALSE]
    
    # Componente 1: (omega(a)/g(a|W)) * (Y - Q*(a,W)) * phi_sel(a) * alpha
    # Q*(a, W_i) para cada observacion
    Q_star_aW <- Q_star(rep(a, n), W)
    g_aW <- g_fun(rep(a, n), W)
    omega_a <- omega_fun(a)
    
    # Pesos
    weights <- omega_a / g_aW
    
    # Residuos
    residuals <- Y - Q_star_aW
    
    # Contribucion de cada observacion
    contrib <- weights * residuals * as.numeric(Phi_a %*% alpha)
    
    # Componente 2: -psi_hat(a) (para que E[D*] = 0)
    ic_matrix[, j] <- contrib - psi_hat[j]
  }
  
  ic_matrix
}

#' Estimador de varianza puntual basado en la curva de influencia
#' 
#' Var(psi_hat(a)) ≈ P_n[D*(a)^2] / n
#' 
#' @param ic_matrix Matriz de curva de influencia (n x m)
#' @export
variance_estimator <- function(ic_matrix) {
  n <- nrow(ic_matrix)
  # Varianza muestral de D*(a) para cada punto de la grid
  colMeans(ic_matrix^2) / n
}

#' Intervalos de confianza de Wald puntuales
#' 
#' @param psi_hat Vector de estimaciones
#' @param var_hat Vector de varianzas estimadas
#' @param alpha Nivel de significancia (default 0.05 para IC 95%)
#' @export
wald_ci <- function(psi_hat, var_hat, alpha = 0.05) {
  z <- qnorm(1 - alpha/2)
  se <- sqrt(var_hat)
  list(
    lower = psi_hat - z * se,
    upper = psi_hat + z * se,
    estimate = psi_hat,
    se = se
  )
}
