# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

mini <- function() {
  set.seed(5)
  d <- simulate_data(120, 1, "normal")
  list(d = d, grid = seq(2, 8, length.out = 5))
}

test_that("thal_fit: inyectar g y omega estimadas reproduce el resultado por defecto", {
  m <- mini()
  r0 <- thal_fit(m$d$A, m$d$W, m$d$Y, m$grid, n_proj = 60L, w_cap = 3, seed = 1)
  r1 <- thal_fit(m$d$A, m$d$W, m$d$Y, m$grid, n_proj = 60L, w_cap = 3, seed = 1,
                 g_fun = fit_g_native(m$d$A, m$d$W, n_bins = 5L),
                 omega_fun = fit_omega_kde_native(m$d$A))
  expect_equal(r1$psi_hat, r0$psi_hat)
  expect_equal(r1$ci$se, r0$ci$se)
})

test_that("thal_fit: acepta g y omega oraculo y da estimaciones finitas", {
  m <- mini()
  r <- thal_fit(m$d$A, m$d$W, m$d$Y, m$grid, n_proj = 60L, w_cap = 3, seed = 1,
                g_fun = function(a, w) g0_dgp(a, w, "normal"),
                omega_fun = function(a) omega_oracle(a, "normal"))
  expect_true(all(is.finite(r$psi_hat)))
  expect_true(all(is.finite(r$ci$se)))
})

test_that("thal_fit: n_cf = 2 no cambia psi_hat y su IC en muestra coincide con el de n_cf = 1", {
  m <- mini()
  r1 <- thal_fit(m$d$A, m$d$W, m$d$Y, m$grid, n_proj = 60L, w_cap = 3, seed = 1, n_cf = 1L)
  r2 <- thal_fit(m$d$A, m$d$W, m$d$Y, m$grid, n_proj = 60L, w_cap = 3, seed = 1, n_cf = 2L)
  expect_equal(r2$psi_hat, r1$psi_hat)
  expect_null(r1$ci_insample)
  expect_equal(r2$ci_insample$se, r1$ci$se)
  expect_true(all(is.finite(r2$ci$se)))
  expect_error(thal_fit(m$d$A, m$d$W, m$d$Y, m$grid, n_cf = 3L), "n_cf")
})

test_that("thal_fit: se^2 = var_R + var_S + 2 cov_RS (componentes del IC)", {
  m <- mini()
  for (ncf in c(1L, 2L)) {
    r <- thal_fit(m$d$A, m$d$W, m$d$Y, m$grid, n_proj = 60L, w_cap = 3, seed = 1, n_cf = ncf)
    p <- r$ic_var_parts
    expect_equal(r$ci$se^2, p$var_R + p$var_S + 2 * p$cov_RS)
  }
})

