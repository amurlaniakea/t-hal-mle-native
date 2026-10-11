# Copyright (C) 2026 Pedro Sordo Martinez
# SPDX-License-Identifier: AGPL-3.0-or-later

#' Estimador nativo de la densidad marginal de A (omega_n) via KDE Sheather-Jones
#'
#' Reproduce la referencia de thalmle: KDE con ancho de banda Sheather-Jones,
#' interpolacion lineal sobre una rejilla y extrapolacion constante fuera del rango.
#' Implementado en R base puro (stats::density + stats::approxfun), sin dependencias.
#'
#' @param A Vector numerico de dosis observadas.
#' @param n_grid Numero de puntos de la rejilla interna de la KDE.
#' @param floor Suelo minimo de densidad (seguridad numerica para los pesos omega/g).
#' @return Funcion omega_fun(a) vectorizada que evalua la densidad marginal.
#' @export
fit_omega_kde_native <- function(A, n_grid = 512L, floor = 1e-4) {
  d <- stats::density(A, bw = "SJ", n = n_grid)
  f <- stats::approxfun(d$x, d$y, rule = 2)
  function(a) pmax(f(a), floor)
}
