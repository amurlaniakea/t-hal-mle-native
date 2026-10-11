# Copyright (C) 2026 Pedro Sordo Martinez
# SPDX-License-Identifier: AGPL-3.0-or-later

#' Pseudo-inversa de Moore-Penrose para matrices simetricas (SVD truncado por autovalores)
#' @param M Matriz simetrica.
#' @param tol Tolerancia relativa respecto al autovalor maximo.
#' @export
pinv_sym <- function(M, tol = 1e-8) {
  e <- eigen(M, symmetric = TRUE)
  vals <- e$values
  vectors <- e$vectors
  pos <- vals > tol * max(vals)
  
  if (!any(pos)) return(matrix(0, nrow = NROW(M), ncol = NCOL(M)))
  
  vectors[, pos, drop = FALSE] %*% 
    diag(1/vals[pos], nrow = sum(pos)) %*% 
    t(vectors[, pos, drop = FALSE])
}

#' Selector de lambda basado en el tamano del sieve (Seccion 5.1 del paper)
#' Busca que el numero de funciones activas caiga en [c1*n^(1/5), c2*n^(1/5)]
#' @export
select_lambda_by_size <- function(df, n, k, c1 = 6, c2 = 9) {
  target_min <- ceiling(c1 * n^(1/5))
  target_max <- floor(c2 * n^(1/5))
  
  in_window <- df >= target_min & df <= target_max
  if (any(in_window)) {
    return(which(in_window)[1])
  }
  
  # Si ningun lambda cae en la ventana, devolver el mas cercano al punto medio
  mid <- (target_min + target_max) / 2
  which.min(abs(df - mid))
}
