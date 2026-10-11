# Copyright (C) 2026 Pedro Sordo Martinez
# SPDX-License-Identifier: AGPL-3.0-or-later

#' Camino de regularizacion LASSO mediante Coordinate Descent (R base puro)
#'
#' Minimiza: (1/(2n)) * ||y - offset - X beta||^2 + lambda * sum(|beta|)
#' Mantiene el vector de residuo para coste O(n*p) por barrido en vez de O(n*p^2).
#'
#' @param X Matriz de diseno (n x p).
#' @param y Vector respuesta.
#' @param offset Vector de offset (longitud n). Por defecto ceros.
#' @param nlambda Numero de lambdas si no se da `lambda`.
#' @param lambda_min_ratio Ratio lambda_min / lambda_max.
#' @param lambda Rejilla de lambda opcional (sobrescribe nlambda/lambda_min_ratio).
#' @param max_iter Iteraciones maximas de CD por lambda.
#' @param tol Tolerancia de convergencia.
#' @return Lista con lambda, beta (matriz p x n_lambda), df (grados de libertad).
#' @export
lasso_path <- function(X, y, offset = rep(0, NROW(X)), nlambda = 50L,
                       lambda_min_ratio = 1e-3, lambda = NULL,
                       max_iter = 200L, tol = 1e-6) {
  n <- NROW(X)
  p <- NCOL(X)
  r0 <- y - offset

  # lambda_max: el valor de lambda que hace todos los coeficientes cero
  lambda_max <- max(abs(crossprod(X, r0))) / n

  # Generar rejilla de lambda si no se proporciona
  if (is.null(lambda)) {
    lambda_min <- lambda_max * lambda_min_ratio
    lambda <- exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))
  } else {
    lambda <- sort(lambda, decreasing = TRUE)
  }

  # Precomputar normas al cuadrado de cada columna
  z <- colSums(X^2)

  # Inicializacion
  beta <- rep(0, p)
  r <- r0  # residuo actual: y - offset - X %*% beta

  n_lambda <- length(lambda)
  beta_path <- matrix(0, nrow = p, ncol = n_lambda)

  for (l in seq_len(n_lambda)) {
    lam <- lambda[l]
    for (iter in seq_len(max_iter)) {
      max_change <- 0
      for (j in seq_len(p)) {
        if (z[j] == 0) next
        beta_j_old <- beta[j]

        # Retirar contribucion antigua del residuo
        if (beta_j_old != 0) {
          r <- r + X[, j] * beta_j_old
        }

        # rho = X_j^T * r
        rho <- sum(X[, j] * r)

        # Soft thresholding
        if (rho > n * lam) {
          beta[j] <- (rho - n * lam) / z[j]
        } else if (rho < -n * lam) {
          beta[j] <- (rho + n * lam) / z[j]
        } else {
          beta[j] <- 0
        }

        # Anadir nueva contribucion al residuo
        if (beta[j] != 0) {
          r <- r - X[, j] * beta[j]
        }

        change <- abs(beta[j] - beta_j_old)
        if (change > max_change) max_change <- change
      }
      if (max_change < tol) break
    }
    beta_path[, l] <- beta
  }

  # Grados de libertad: numero de coeficientes no nulos
  df <- colSums(beta_path != 0)

  list(lambda = lambda, beta = beta_path, df = df)
}
