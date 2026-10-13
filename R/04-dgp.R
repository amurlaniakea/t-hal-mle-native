# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

# Cache de la densidad marginal exacta de A (caso "normal"); se calcula una sola vez.
.dgp_cache <- new.env(parent = emptyenv())

#' Componente g(w) de los DGP del paper (seccion 7.1)
#'
#' DGP 1: `5 * 1[w > 0]`. DGP 2: `4 (1[w > -1] - 1[w > 0] + 1[w > 1])`.
#' DGP 3: `3 sin(7 pi w)`.
#' @param w Covariable.
#' @param dgp Numero de DGP (1, 2 o 3).
#' @keywords internal
dgp_g_w <- function(w, dgp) {
  switch(as.character(dgp),
    "1" = 5 * (w > 0),
    "2" = 4 * ((w > -1) - (w > 0) + (w > 1)),
    "3" = 3 * sin(7 * pi * w),
    stop("dgp debe ser 1, 2 o 3"))
}

#' Simula datos con los DGP del paper (arXiv:2607.03824, seccion 7.1)
#'
#' `W ~ N(0,1)`; `Y = 3 sin(1.5 A) + 0.5 A + g(W) + (A - 5) g(W) + N(0,1)`, con `g` segun el DGP
#' (ver [dgp_g_w()]). Tratamiento `"uniform"`: `A ~ U(0, 10)`. Tratamiento `"normal"`:
#' `A | W ~ N(5 + 1.5 W, 2^2)` truncada a `[0, 10]` (muestreo por inversion, sin masa en los bordes).
#' @param n Tamano de muestra.
#' @param dgp Numero de DGP (1, 2 o 3).
#' @param treatment Tipo de tratamiento: "uniform" o "normal".
#' @return Lista con `W`, `A` e `Y`.
#' @export
simulate_data <- function(n, dgp = 1, treatment = c("uniform", "normal")) {
  treatment <- match.arg(treatment)
  W <- rnorm(n)
  if (treatment == "uniform") {
    A <- runif(n, 0, 10)
  } else {
    mu <- 5 + 1.5 * W
    A <- mu + 2 * qnorm(runif(n, pnorm((0 - mu) / 2), pnorm((10 - mu) / 2)))
  }
  gw <- dgp_g_w(W, dgp)
  Y <- 3 * sin(1.5 * A) + 0.5 * A + gw + (A - 5) * gw + rnorm(n)
  list(W = W, A = A, Y = Y)
}

#' Regresion verdadera Q0(a, w) = E[Y | A = a, W = w] del DGP
#' @param a Dosis.
#' @param w Covariable.
#' @param dgp Numero de DGP (1, 2 o 3).
#' @export
Q0_dgp <- function(a, w, dgp = 1) {
  gw <- dgp_g_w(w, dgp)
  3 * sin(1.5 * a) + 0.5 * a + gw + (a - 5) * gw
}

#' Densidad condicional verdadera g0(a | w); vale 0 fuera de [0, 10]
#' @param a Dosis.
#' @param w Covariable.
#' @param treatment "uniform" o "normal" (normal truncada a [0, 10]).
#' @export
g0_dgp <- function(a, w, treatment = c("uniform", "normal")) {
  treatment <- match.arg(treatment)
  inside <- a >= 0 & a <= 10
  if (treatment == "uniform") {
    dens <- rep(1 / 10, length(a))
  } else {
    mu <- 5 + 1.5 * w
    dens <- dnorm(a, mean = mu, sd = 2) / (pnorm(10, mean = mu, sd = 2) - pnorm(0, mean = mu, sd = 2))
  }
  dens * inside
}

#' Densidad marginal verdadera de A (oraculo exacto)
#'
#' Caso "normal": integral numerica de `g0(a | w) dN(w)` sobre una rejilla fina de `w`,
#' tabulada una vez e interpolada.
#' @param a Dosis.
#' @param treatment "uniform" o "normal".
#' @export
omega_oracle <- function(a, treatment = c("uniform", "normal")) {
  treatment <- match.arg(treatment)
  if (treatment == "uniform") return(rep(1 / 10, length(a)) * (a >= 0 & a <= 10))
  if (is.null(.dgp_cache$omega)) {
    wg <- seq(-7, 7, by = 0.01)
    ww <- dnorm(wg) * 0.01
    ag <- seq(0, 10, length.out = 1001)
    vals <- vapply(ag, function(x) sum(g0_dgp(rep(x, length(wg)), wg, "normal") * ww), numeric(1))
    .dgp_cache$omega <- stats::approxfun(ag, vals, rule = 2)
  }
  .dgp_cache$omega(a) * (a >= 0 & a <= 10)
}

#' Curva dosis-respuesta verdadera psi0(a) del DGP
#'
#' `psi0(a) = 3 sin(1.5 a) + 0.5 a + E[g(W)] (a - 4)`, con `E[g(W)] = 2.5, 2, 0` en los DGP 1, 2, 3.
#' @param a Dosis.
#' @param dgp Numero de DGP (1, 2 o 3).
#' @export
true_drc <- function(a, dgp = 1) {
  Eg <- switch(as.character(dgp),
    "1" = 5 * 0.5,
    "2" = 4 * (pnorm(1) - 0.5 + (1 - pnorm(1))),
    "3" = 0,
    stop("dgp debe ser 1, 2 o 3"))
  3 * sin(1.5 * a) + 0.5 * a + Eg * (a - 4)
}
