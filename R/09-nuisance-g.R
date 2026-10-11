# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

#' Estimador nativo de la densidad condicional g(a | w) via KDE por bins
#'
#' Alternativa nativa a haldensify (Q-04). Divide W en K bins por cuantiles
#' y ajusta una KDE 1D (Sheather-Jones) de A dentro de cada bin.
#' Resuelve Q-14: es una funcion general de (a, w) que evalua fuera de muestra.
#'
#' @param A Vector de dosis observadas.
#' @param W Vector de covariable observada.
#' @param n_bins Numero de bins para W (por defecto 5).
#' @param floor Suelo minimo de densidad (seguridad numerica para pesos omega/g).
#' @return Funcion g_fun(a, w) vectorizada por pares.
#' @export
fit_g_native <- function(A, W, n_bins = 5L, floor = 1e-4) {
  # 1. Crear breakpoints para los bins de W
  probs <- seq(0, 1, length.out = n_bins + 1)
  breaks <- quantile(W, probs = probs, names = FALSE)
  breaks[1] <- -Inf
  breaks[n_bins + 1] <- Inf
  
  bin_id <- cut(W, breaks = breaks, labels = FALSE, include.lowest = TRUE)
  
  # 2. Precomputar approxfuns de KDE para cada bin
  afuns <- lapply(seq_len(n_bins), function(b) {
    A_b <- A[bin_id == b]
    if (length(A_b) < 5) {
      # Fallback a la marginal si hay muy pocos datos en el bin
      d <- stats::density(A, bw = "SJ", n = 256)
    } else {
      d <- stats::density(A_b, bw = "SJ", n = 256)
    }
    stats::approxfun(d$x, d$y, rule = 2)
  })
  
  # 3. Funcion de evaluacion vectorizada por pares (a, w)
  function(a, w) {
    w_bin <- cut(w, breaks = breaks, labels = FALSE, include.lowest = TRUE)
    w_bin[is.na(w_bin)] <- 1L # Fallback por seguridad
    
    mapply(function(ai, bi) {
      pmax(afuns[[bi]](ai), floor)
    }, a, w_bin, USE.NAMES = FALSE)
  }
}
