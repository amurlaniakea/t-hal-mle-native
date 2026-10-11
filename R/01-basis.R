# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

#' Genera nudos estrictamente interiores para evitar colinealidad con los bordes
#' @export
make_knots <- function(x, n_knots) {
  if (n_knots <= 0) return(numeric(0))
  probs <- seq(0, 1, length.out = n_knots + 2)[-c(1, n_knots + 2)]
  sort(unique(quantile(x, probs = probs, names = FALSE, type = 1L)))
}

#' Construye la base de splines de orden k segun van der Laan (2023)
#' Escala x a [0,1] para que el "borde 0" tenga sentido.
#' @export
spline_basis <- function(x, knots, k = 1L, x_range = range(x, na.rm = TRUE)) {
  min_x <- x_range[1]; max_x <- x_range[2]
  if (max_x == min_x) {
    x_scaled <- rep(0, length(x))
    knots_scaled <- rep(0, length(knots))
  } else {
    x_scaled <- (x - min_x) / (max_x - min_x)
    knots_scaled <- (knots - min_x) / (max_x - min_x)
  }
  n <- length(x_scaled)
  if (k == 0L) {
    Phi <- matrix(1, nrow = n, ncol = 1 + length(knots_scaled))
    for (j in seq_along(knots_scaled)) {
      Phi[, 1 + j] <- as.numeric(x_scaled >= knots_scaled[j])
    }
    return(Phi)
  }
  if (k == 1L) {
    Phi <- matrix(0, nrow = n, ncol = 2 + length(knots_scaled))
    Phi[, 1] <- 1
    Phi[, 2] <- x_scaled
    for (j in seq_along(knots_scaled)) {
      Phi[, 2 + j] <- pmax(0, x_scaled - knots_scaled[j])
    }
    return(Phi)
  }
  stop("Solo k=0 y k=1 implementados en v0.")
}
