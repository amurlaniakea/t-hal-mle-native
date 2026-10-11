# Copyright (C) 2026 Pedro Sordo Martinez
# SPDX-License-Identifier: AGPL-3.0-or-later

#' Genera nudos incluyendo SIEMPRE el minimo de A para evitar el hueco (Q-17)
#' @param A Vector numerico de dosis observadas.
#' @param n_knots Numero de nudos internos.
#' @export
make_knots <- function(A, n_knots) {
  u_min <- min(A)
  if (n_knots <= 1) return(u_min)
  
  # Cuantiles del rango, excluyendo el 0% para no duplicar el minimo exacto
  probs <- seq(0, 1, length.out = n_knots + 1)[-1] 
  u_rest <- quantile(A, probs = probs, names = FALSE, type = 1L)
  
  sort(unique(c(u_min, u_rest)))
}

#' Construye la base de splines de orden k segun van der Laan (2023)
#' Para k=1: 1, a (borde 0), y (a - u_j)_+
#' @param x Vector numerico.
#' @param knots Vector de nudos generado por make_knots.
#' @param k Orden de suavidad (0 o 1).
#' @export
spline_basis <- function(x, knots, k = 1L) {
  n <- length(x)
  
  if (k == 0L) {
    # Base orden 0: intercepto + I(x >= u_j)
    Phi <- matrix(1, nrow = n, ncol = 1 + length(knots))
    for (j in seq_along(knots)) {
      Phi[, 1 + j] <- as.numeric(x >= knots[j])
    }
    return(Phi)
  }
  
  if (k == 1L) {
    # Base orden 1: intercepto (1), termino lineal puro (x), y (x - u_j)_+
    # El termino lineal puro es el spline de orden 1 con nudo en el origen (borde 0)
    Phi <- matrix(0, nrow = n, ncol = 2 + length(knots))
    Phi[, 1] <- 1
    Phi[, 2] <- x 
    
    for (j in seq_along(knots)) {
      Phi[, 2 + j] <- pmax(0, x - knots[j])
    }
    return(Phi)
  }
  
  stop("Solo k=0 y k=1 implementados en v0.")
}
