test_that("true_drc coincide con la integral numerica de Q0 sobre W ~ N(0,1)", {
  wg <- seq(-8, 8, by = 0.005); ww <- dnorm(wg) * 0.005
  a <- c(1, 2.5, 5, 7.5, 9)
  for (dgp in 1:3) {
    num <- vapply(a, function(x) sum(Q0_dgp(rep(x, length(wg)), wg, dgp) * ww), numeric(1))
    expect_equal(true_drc(a, dgp), num, tolerance = 5e-3)
  }
})

test_that("simulate_data respeta soporte y estructura de los DGP del paper", {
  set.seed(1)
  for (trt in c("uniform", "normal")) {
    d <- simulate_data(5000, 1, trt)
    expect_true(all(d$A >= 0 & d$A <= 10))
    expect_equal(mean(d$W), 0, tolerance = 0.1)
  }
  # con tratamiento normal hay confusion: A depende de W
  d <- simulate_data(5000, 1, "normal")
  expect_gt(cor(d$A, d$W), 0.3)
  # Y - Q0(A, W) es ruido N(0,1)
  for (dgp in 1:3) {
    d <- simulate_data(5000, dgp, "uniform")
    e <- d$Y - Q0_dgp(d$A, d$W, dgp)
    expect_equal(sd(e), 1, tolerance = 0.05)
  }
})

test_that("g0_dgp y omega_oracle son densidades que integran a 1 en [0, 10]", {
  ag <- seq(0, 10, length.out = 4001); h <- ag[2] - ag[1]
  for (w in c(-1.5, 0, 2)) expect_equal(sum(g0_dgp(ag, w, "normal")) * h, 1, tolerance = 1e-3)
  expect_equal(sum(omega_oracle(ag, "normal")) * h, 1, tolerance = 1e-3)
  expect_equal(sum(omega_oracle(ag, "uniform")) * h, 1, tolerance = 1e-3)
  expect_equal(g0_dgp(c(-1, 11), 0, "normal"), c(0, 0))
  # omega es la media sobre W de g0
  wg <- seq(-8, 8, by = 0.005); ww <- dnorm(wg) * 0.005
  num <- sum(g0_dgp(rep(9, length(wg)), wg, "normal") * ww)
  expect_equal(omega_oracle(9, "normal"), num, tolerance = 1e-3)
})
