test_that("thal_fit funciona con DGP 1, tratamiento uniforme", {
  set.seed(100)
  d <- simulate_data(300, dgp = 1, treatment = "uniform")
  grid <- seq(2, 8, length.out = 10)
  
  res <- thal_fit(
    A = d$A, W = d$W, Y = d$Y, grid = grid,
    k_targeting = 1L, k_Q = 1L,
    n_knots_targeting = 20, n_knots_Q = 8,
    n_bins_g = 5,
    w_cap = "obs_max",
    seed = 42
  )
  
  expect_true(is.list(res))
  expect_true("psi_hat" %in% names(res))
  expect_true("ci" %in% names(res))
  expect_true("ic_matrix" %in% names(res))
  expect_equal(length(res$psi_hat), length(grid))
  
  truth <- true_drc(grid, dgp = 1)
  rmse <- sqrt(mean((res$psi_hat - truth)^2))
  expect_true(is.finite(rmse))
  expect_lt(rmse, 10.0)
  
  ci_width <- mean(res$ci$upper - res$ci$lower)
  expect_gt(ci_width, 0.01)
  expect_lt(ci_width, 50)
})

test_that("thal_fit funciona con DGP 1, tratamiento normal", {
  set.seed(200)
  d <- simulate_data(300, dgp = 1, treatment = "normal")
  grid <- seq(2, 8, length.out = 10)
  
  res <- thal_fit(
    A = d$A, W = d$W, Y = d$Y, grid = grid,
    k_targeting = 1L, k_Q = 1L,
    n_knots_targeting = 20, n_knots_Q = 8,
    n_bins_g = 5,
    w_cap = "obs_max",
    seed = 42
  )
  
  truth <- true_drc(grid, dgp = 1)
  rmse <- sqrt(mean((res$psi_hat - truth)^2))
  
  expect_true(is.finite(rmse))
  expect_lt(rmse, 40.0)
})

test_that("thal_fit funciona con DGP 2 (discontinuidad), k_Q = 0", {
  set.seed(300)
  d <- simulate_data(300, dgp = 2, treatment = "uniform")
  grid <- seq(2, 8, length.out = 10)
  
  res <- thal_fit(
    A = d$A, W = d$W, Y = d$Y, grid = grid,
    k_targeting = 1L, k_Q = 0L,
    n_knots_targeting = 30, n_knots_Q = 10,
    n_bins_g = 5,
    w_cap = "obs_max",
    seed = 42
  )
  
  truth <- true_drc(grid, dgp = 2)
  rmse <- sqrt(mean((res$psi_hat - truth)^2))
  
  expect_true(is.finite(rmse))
  expect_lt(rmse, 10.0)
})
