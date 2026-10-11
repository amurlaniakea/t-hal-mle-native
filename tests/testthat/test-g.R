test_that("fit_g_native devuelve valores positivos", {
  set.seed(10)
  n <- 200
  A <- runif(n, 1, 9)
  W <- runif(n, 0, 1)
  
  g_fun <- fit_g_native(A, W, n_bins = 5)
  
  # Evaluar en una rejilla
  grid_a <- seq(2, 8, length.out = 10)
  grid_w <- seq(0.2, 0.8, length.out = 10)
  g_vals <- g_fun(rep(grid_a, each = 10), rep(grid_w, times = 10))
  
  expect_true(all(g_vals > 0))
  expect_true(all(is.finite(g_vals)))
})

test_that("fit_g_native: cuando A y W son independientes, g(a|w) ~ omega(a)", {
  set.seed(20)
  n <- 300
  A <- runif(n, 1, 9)
  W <- runif(n, 0, 1)  # Independiente de A
  
  g_fun <- fit_g_native(A, W, n_bins = 5)
  omega_fun <- fit_omega_kde_native(A)
  
  # Evaluar en puntos intermedios (evitar bordes donde KDE tiene sesgo)
  grid_a <- seq(2, 8, length.out = 20)
  grid_w <- rep(0.5, 20)
  
  g_vals <- g_fun(grid_a, grid_w)
  omega_vals <- omega_fun(grid_a)
  
  # Cuando A y W son independientes, g(a|w) debe ser cercano a omega(a)
  # Tolerancia amplia porque son estimadores diferentes
  rmse <- sqrt(mean((g_vals - omega_vals)^2))
  expect_lt(rmse, 0.15)
})

test_that("fit_g_native evalua fuera de muestra (resuelve Q-14)", {
  set.seed(30)
  n <- 200
  A <- runif(n, 1, 9)
  W <- runif(n, 0, 1)
  
  g_fun <- fit_g_native(A, W, n_bins = 5)
  
  # Evaluar con valores de W que no estaban en la muestra
  new_w <- c(0.15, 0.35, 0.55, 0.75, 0.95)
  new_a <- rep(5, length(new_w))
  
  g_vals <- g_fun(new_a, new_w)
  
  expect_true(all(g_vals > 0))
  expect_true(all(is.finite(g_vals)))
})
