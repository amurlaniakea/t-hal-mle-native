#' Simula datos segun los DGP del paper (Tabla 1)
#' @param n Tamano de muestra
#' @param dgp Numero de DGP (1, 2 o 3)
#' @param treatment Tipo de tratamiento: "uniform" o "normal"
#' @export
simulate_data <- function(n, dgp = 1, treatment = c("uniform", "normal")) {
  treatment <- match.arg(treatment)
  
  # W ~ Uniform(0, 1) en todos los DGP
  W <- runif(n, 0, 1)
  
  # A segun treatment
  if (treatment == "uniform") {
    A <- runif(n, 1, 9)
  } else {
    # Normal centrada en 5, sd = 2, truncada a [1, 9]
    A <- pmin(pmax(rnorm(n, mean = 5, sd = 2), 1), 9)
  }
  
  # Outcome Y segun DGP
  if (dgp == 1) {
    # Q0(a, w) = 3*sin(1.5*a) + 0.5*a + 2*w
    Y <- 3 * sin(1.5 * A) + 0.5 * A + 2 * W + rnorm(n, sd = 1)
  } else if (dgp == 2) {
    # Q0(a, w) con discontinuidad en a = 5
    Y <- ifelse(A < 5, 2 * A + W, 10 - A + 2 * W) + rnorm(n, sd = 1)
  } else if (dgp == 3) {
    # Q0(a, w) suave de orden superior
    Y <- sin(A) + cos(2 * A) + 0.5 * W * A + rnorm(n, sd = 1)
  } else {
    stop("dgp debe ser 1, 2 o 3")
  }
  
  list(W = W, A = A, Y = Y)
}

#' Outcome regression verdadero Q0(a, w)
#' @export
Q0_dgp <- function(a, w, dgp = 1) {
  if (dgp == 1) {
    3 * sin(1.5 * a) + 0.5 * a + 2 * w
  } else if (dgp == 2) {
    ifelse(a < 5, 2 * a + w, 10 - a + 2 * w)
  } else if (dgp == 3) {
    sin(a) + cos(2 * a) + 0.5 * w * a
  } else {
    stop("dgp debe ser 1, 2 o 3")
  }
}

#' Densidad condicional verdadera g(a | w)
#' @export
g0_dgp <- function(a, w, treatment = c("uniform", "normal")) {
  treatment <- match.arg(treatment)
  
  if (treatment == "uniform") {
    # Uniforme en [1, 9], densidad = 1/8
    rep(1/8, length(a))
  } else {
    # Normal truncada en [1, 9], media 5, sd 2
    # g(a) = phi((a-5)/2) / (2 * (Phi(2) - Phi(-2)))
    z <- (a - 5) / 2
    phi_z <- dnorm(z)
    Phi_2 <- pnorm(2)
    Phi_minus_2 <- pnorm(-2)
    phi_z / (2 * (Phi_2 - Phi_minus_2))
  }
}

#' Densidad marginal verdadera omega(a)
#' @export
omega_oracle <- function(a, treatment = c("uniform", "normal")) {
  treatment <- match.arg(treatment)
  
  if (treatment == "uniform") {
    rep(1/8, length(a))
  } else {
    g0_dgp(a, w = NULL, treatment = "normal")
  }
}

#' Curva dosis-respuesta verdadera psi(a) = E_W[Q0(a, W)]
#' @export
true_drc <- function(a, dgp = 1) {
  # E_W[Q0(a, W)] con W ~ Uniform(0, 1)
  if (dgp == 1) {
    # E[3*sin(1.5*a) + 0.5*a + 2*W] = 3*sin(1.5*a) + 0.5*a + 2*E[W] = 3*sin(1.5*a) + 0.5*a + 1
    3 * sin(1.5 * a) + 0.5 * a + 1
  } else if (dgp == 2) {
    # E[ifelse(a < 5, 2*a + W, 10 - a + 2*W)] = ifelse(a < 5, 2*a + 0.5, 10 - a + 1)
    ifelse(a < 5, 2 * a + 0.5, 10 - a + 1)
  } else if (dgp == 3) {
    # E[sin(a) + cos(2*a) + 0.5*W*a] = sin(a) + cos(2*a) + 0.25*a
    sin(a) + cos(2 * a) + 0.25 * a
  } else {
    stop("dgp debe ser 1, 2 o 3")
  }
}
