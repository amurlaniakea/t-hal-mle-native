# Copyright (C) 2026 Pedro Sordo Martinez
# SPDX-License-Identifier: AGPL-3.0-or-later

#' Base HAL bivariada (d=2) de orden k, segun van der Laan (2023)
#'
#' k=0: 1, I(a>=ua), I(w>=uw), I(a>=ua)I(w>=uw)
#' k=1: 1, a, w, a*w, (a-ua)+, (w-uw)+, w(a-ua)+, a(w-uw)+, (a-ua)+(w-uw)+
#' Los terminos de borde (orden < k) son potencias (nudo en el origen), no indicadores.
#'
#' @export
tensor_basis <- function(A, W, knots_A, knots_W, k = 0L) {
  if (length(A) != length(W)) stop("A y W deben tener la misma longitud")
  if (k == 0L) return(tensor_basis_k0(A, W, knots_A, knots_W))
  if (k == 1L) return(tensor_basis_k1(A, W, knots_A, knots_W))
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
