# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

#' Base HAL bivariada (d=2) de orden k, segun van der Laan (2023)
#' Escala A y W a [0,1] para que los terminos de borde tengan sentido.
#'
#' k=0: 1, I(a_scaled>=ua), I(w_scaled>=uw), I(a_scaled>=ua)I(w_scaled>=uw)
#' k=1: 1, a_scaled, w_scaled, a_scaled*w_scaled, (a_scaled-ua)+, (w_scaled-uw)+,
#'      w_scaled(a_scaled-ua)+, a_scaled(w_scaled-uw)+, (a_scaled-ua)+(w_scaled-uw)+
#'
#' @export
tensor_basis <- function(A, W, knots_A, knots_W, k = 0L,
                         A_range = range(A, na.rm = TRUE),
                         W_range = range(W, na.rm = TRUE)) {
  if (length(A) != length(W)) stop("A y W deben tener la misma longitud")
  
  # Escalar A a [0,1]
  min_A <- A_range[1]; max_A <- A_range[2]
  if (max_A == min_A) {
    A_scaled <- rep(0, length(A))
    knots_A_scaled <- rep(0, length(knots_A))
  } else {
    A_scaled <- (A - min_A) / (max_A - min_A)
    knots_A_scaled <- (knots_A - min_A) / (max_A - min_A)
  }
  
  # Escalar W a [0,1]
  min_W <- W_range[1]; max_W <- W_range[2]
  if (max_W == min_W) {
    W_scaled <- rep(0, length(W))
    knots_W_scaled <- rep(0, length(knots_W))
  } else {
    W_scaled <- (W - min_W) / (max_W - min_W)
    knots_W_scaled <- (knots_W - min_W) / (max_W - min_W)
  }
  
  if (k == 0L) return(tensor_basis_k0(A_scaled, W_scaled, knots_A_scaled, knots_W_scaled))
  if (k == 1L) return(tensor_basis_k1(A_scaled, W_scaled, knots_A_scaled, knots_W_scaled))
  stop("Solo k=0 y k=1 implementados")
}

tensor_basis_k0 <- function(A, W, knots_A, knots_W) {
  n <- length(A); nA <- length(knots_A); nW <- length(knots_W)
  ind_A <- matrix(0, n, nA)
  for (j in seq_along(knots_A)) ind_A[, j] <- as.numeric(A >= knots_A[j])
  ind_W <- matrix(0, n, nW)
  for (j in seq_along(knots_W)) ind_W[, j] <- as.numeric(W >= knots_W[j])

  Phi <- matrix(0, n, 1 + nA + nW + nA * nW)
  col <- 1
  Phi[, col] <- 1; col <- col + 1
  Phi[, col:(col + nA - 1)] <- ind_A; col <- col + nA
  Phi[, col:(col + nW - 1)] <- ind_W; col <- col + nW
  for (ja in seq_len(nA)) for (jw in seq_len(nW)) {
    Phi[, col] <- ind_A[, ja] * ind_W[, jw]; col <- col + 1
  }
  Phi
}

tensor_basis_k1 <- function(A, W, knots_A, knots_W) {
  n <- length(A); nA <- length(knots_A); nW <- length(knots_W)
  sp_A <- matrix(0, n, nA)
  for (j in seq_along(knots_A)) sp_A[, j] <- pmax(0, A - knots_A[j])
  sp_W <- matrix(0, n, nW)
  for (j in seq_along(knots_W)) sp_W[, j] <- pmax(0, W - knots_W[j])

  Phi <- matrix(0, n, 4 + 2 * nA + 2 * nW + nA * nW)
  col <- 1
  Phi[, col] <- 1;     col <- col + 1
  Phi[, col] <- A;     col <- col + 1
  Phi[, col] <- W;     col <- col + 1
  Phi[, col] <- A * W; col <- col + 1
  Phi[, col:(col + nA - 1)] <- sp_A; col <- col + nA
  Phi[, col:(col + nW - 1)] <- sp_W; col <- col + nW
  for (j in seq_len(nA)) { Phi[, col] <- W * sp_A[, j]; col <- col + 1 }
  for (j in seq_len(nW)) { Phi[, col] <- A * sp_W[, j]; col <- col + 1 }
  for (ja in seq_len(nA)) for (jw in seq_len(nW)) {
    Phi[, col] <- sp_A[, ja] * sp_W[, jw]; col <- col + 1
  }
  Phi
}
