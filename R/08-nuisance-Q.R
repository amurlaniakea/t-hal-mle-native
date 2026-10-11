# Copyright (C) 2026 Pedro Sordo Martinez
# SPDX-License-Identifier: AGPL-3.0-or-later

#' Cross-validacion V-fold sobre el camino LASSO (criterio lambda.min)
#'
#' @param X Matriz de diseno (n x p).
#' @param y Vector respuesta.
#' @param offset Vector de offset (longitud n). Por defecto ceros.
#' @param nfolds Numero de folds.
#' @param nlambda Numero de lambdas en la rejilla comun.
#' @param lambda_min_ratio Ratio lambda_min / lambda_max.
#' @param seed Semilla opcional para los folds.
#' @return Lista con lambda, cv_error, lambda_min_idx, lambda_min.
#' @export
cv_lasso <- function(X, y, offset = rep(0, NROW(X)), nfolds = 5L,
                     nlambda = 30L, lambda_min_ratio = 1e-3, seed = NULL) {
  n <- NROW(X)
  if (!is.null(seed)) set.seed(seed)
  fold_id <- sample(rep_len(seq_len(nfolds), n))

  # Rejilla comun de lambda a partir de los datos completos
  r0 <- y - offset
  lambda_max <- max(abs(crossprod(X, r0))) / n
  if (lambda_max == 0) lambda_max <- 1
  lambdas <- exp(seq(log(lambda_max), log(lambda_max * lambda_min_ratio),
                     length.out = nlambda))

  cv_err <- matrix(NA_real_, nrow = nlambda, ncol = nfolds)
  for (v in seq_len(nfolds)) {
    tr <- which(fold_id != v)
    va <- which(fold_id == v)
    fit_v <- lasso_path(X[tr, , drop = FALSE], y[tr], offset = offset[tr],
                        lambda = lambdas)
    X_va <- X[va, , drop = FALSE]
    for (l in seq_along(lambdas)) {
      pred <- offset[va] + as.numeric(X_va %*% fit_v$beta[, l])
      cv_err[l, v] <- mean((y[va] - pred)^2)
    }
  }

  cv_mean <- rowMeans(cv_err)
  idx_min <- which.min(cv_mean)
  list(lambda = lambdas, cv_error = cv_mean,
       lambda_min_idx = idx_min, lambda_min = lambdas[idx_min])
}

#' Ajusta el HAL del outcome Q_n(A, W) de forma nativa
#'
#' Base bivariada tensor_basis + LASSO con lambda por CV (lambda.min).
#' Devuelve una funcion Q_fun(a, w) vectorizada por pares.
#' Nota de diseno: el intercepto es una columna mas de la base y se penaliza
#' igual que el resto (coherente con el paso de targeting).
#'
#' @param A Vector de dosis observadas.
#' @param W Vector de covariable observada.
#' @param Y Vector de outcome.
#' @param k Orden de suavidad (0 para DGP 1/2, 1 para DGP 3).
#' @param n_knots Numero de nudos por dimension.
#' @param nfolds Folds para el CV de lambda.
#' @param nlambda Numero de lambdas en el CV.
#' @param seed Semilla.
#' @return Funcion Q_fun(a, w).
#' @export
fit_Q_hal_native <- function(A, W, Y, k = 0L, n_knots = 10L, nfolds = 5L,
                             nlambda = 30L, seed = NULL) {
  knots_A <- make_knots(A, n_knots)
  knots_W <- make_knots(W, n_knots)
  Phi <- tensor_basis(A, W, knots_A, knots_W, k)

  cv <- cv_lasso(Phi, Y, nfolds = nfolds, nlambda = nlambda, seed = seed)
  fit_full <- lasso_path(Phi, Y, lambda = cv$lambda)
  beta <- fit_full$beta[, cv$lambda_min_idx]

  function(a, w) {
    Phi_new <- tensor_basis(a, w, knots_A, knots_W, k)
    as.numeric(Phi_new %*% beta)
  }
}
