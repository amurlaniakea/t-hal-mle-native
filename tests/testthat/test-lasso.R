test_that("lasso_path recupera solucion OLS cuando lambda es muy pequeno", {
  set.seed(123)
  n <- 50; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  beta_true <- c(1, -2, 0, 3, 0)
  y <- X %*% beta_true + rnorm(n, sd = 0.1)
  
  fit <- lasso_path(X, y, nlambda = 20, lambda_min_ratio = 1e-4)
  beta_hat <- fit$beta[, ncol(fit$beta)]
  
  # Comparar con OLS analitico
  beta_ols <- as.numeric(solve(crossprod(X)) %*% crossprod(X, y))
  expect_equal(beta_hat, beta_ols, tolerance = 1e-3)
})

test_that("lasso_path produce coeficientes cero para lambda grande", {
  set.seed(123)
  n <- 50; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  
  fit <- lasso_path(X, y, nlambda = 10)
  beta_first <- fit$beta[, 1]
  
  expect_true(all(abs(beta_first) < 1e-10))
})
