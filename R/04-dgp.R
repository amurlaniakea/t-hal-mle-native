# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

#' Simula datos segun los DGP exactos del paper (Tabla 1 / S6)
#' @param n Tamano de muestra
#' @param dgp Numero de DGP (1, 2 o 3)
#' @param treatment Tipo de tratamiento: "uniform" o "normal"
#' @export
simulate_data <- function(n, dgp = 1, treatment = c("uniform", "normal")) {
  treatment <- match.arg(treatment)
  W <- rnorm(n, mean = 0, sd = 1)
  if (treatment == "uniform") {
    A <- runif(n, min = 0, max = 10)
  } else {
    A <- rnorm(n, mean = 5 + 1.5 * W, sd = 2)
  }
  g_w <- function(w) 0.5 * w + 0.25 * w^2
  if (dgp == 1) {
    Q0 <- 3 * sin(1.5 * A) + 0.5 * A + g_w(W) + (A - 5) * g_w(W)
    Y <- Q0 + rnorm(n, sd = 1)
  } else if (dgp == 2) {
    Q0 <- 3 * sin(1.5 * A) + 0.5 * A + g_w(W) * (W > 0) + (A - 5) * g_w(W) * (W > 0)
    Y <- Q0 + rnorm(n, sd = 1)
  } else if (dgp == 3) {
    Q0 <- 2 * sin(A) + 0.5 * A^2 + g_w(W)
    Y <- Q0 + rnorm(n, sd = 1)
  } else {
    stop("dgp debe ser 1, 2 o 3")
  }
  list(W = W, A = A, Y = Y)
}

#' @export
Q0_dgp <- function(a, w, dgp = 1) {
  g_w <- 0.5 * w + 0.25 * w^2
  if (dgp == 1) 3 * sin(1.5 * a) + 0.5 * a + g_w + (a - 5) * g_w
  else if (dgp == 2) 3 * sin(1.5 * a) + 0.5 * a + g_w * (w > 0) + (a - 5) * g_w * (w > 0)
  else if (dgp == 3) 2 * sin(a) + 0.5 * a^2 + g_w
  else stop("dgp debe ser 1, 2 o 3")
}

#' @export
g0_dgp <- function(a, w, treatment = c("uniform", "normal")) {
  treatment <- match.arg(treatment)
  if (treatment == "uniform") rep(1 / 10, length(a))
  else dnorm(a, mean = 5 + 1.5 * w, sd = 2)
}

#' @export
omega_oracle <- function(a, treatment = c("uniform", "normal")) {
  treatment <- match.arg(treatment)
  if (treatment == "uniform") rep(1 / 10, length(a))
  else dnorm(a, mean = 5, sd = 2.5)
}

#' @export
true_drc <- function(a, dgp = 1) {
  if (dgp == 1) 3 * sin(1.5 * a) + 0.5 * a + 0.25 + (a - 5) * 0.25
  else if (dgp == 2) {
    c_val <- 0.3245
    3 * sin(1.5 * a) + 0.5 * a + c_val + (a - 5) * c_val
  }
  else if (dgp == 3) 2 * sin(a) + 0.5 * a^2 + 0.25
  else stop("dgp debe ser 1, 2 o 3")
}
