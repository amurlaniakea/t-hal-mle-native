# thalnative

Implementacion nativa en R base del Targeted Highly Adaptive Lasso Minimum Loss Estimation (T-HAL-MLE) para curvas dosis-respuesta continuas.

Estado: v0.0.1, 53 tests pasando, 20 funciones exportadas, 100% R base.

## Caracteristicas

- Cero dependencias externas: sin glmnet, hal9001, haldensify, ni Matrix. Solo R base + stats.
- LASSO nativo: Coordinate Descent con warm starts y mantenimiento de residuo.
- Base HAL matematicamente fiel: terminos de borde segun van der Laan (2023), nudos estrictamente interiores, escalado a [0,1] para evitar colinealidad.
- Nuisances nativos: KDE para omega, HAL bivariado con CV para Q, KDE por bins para g.
- Inferencia (experimental): curva de influencia del Lemma 3.1 (componentes Q_bar y Q_W, con Sigma^+) + IC de Wald puntuales. Con nuisances oraculo, DGP 1, n = 300 y 80 replicas la cobertura de Wald medida fue ~0.73 (nominal 0.95): aun no es fiable. Ver SPEC/decisions.md (D-0004).

## Instalacion

cd ~/repos
R CMD build t-hal-mle-native
R CMD INSTALL thalnative_0.0.1.tar.gz

## Uso rapido

library(thalnative)
set.seed(42)
d <- simulate_data(n = 500, dgp = 1, treatment = "uniform")

res <- thal_fit(
  A = d$A, W = d$W, Y = d$Y,
  grid = seq(2, 8, length.out = 20),
  k_targeting = 1L, k_Q = 1L,
  n_knots_targeting = 50, n_knots_Q = 10,
  n_bins_g = 5,
  w_cap = NULL,
  seed = 123
)

data.frame(
  a = res$grid,
  psi_hat = res$psi_hat,
  lower = res$ci$lower,
  upper = res$ci$upper,
  truth = true_drc(res$grid, dgp = 1)
)

## Funciones principales

- thal_fit(): Pipeline completo
- thal_core_native(): Nucleo del targeting
- fit_Q_hal_native(): HAL bivariado del outcome con CV
- fit_g_native(): Densidad condicional g(a|w) por KDE en bins
- fit_omega_kde_native(): Densidad marginal omega(a) por KDE
- lasso_path(): LASSO nativo
- spline_basis(): Base HAL univariada con escalado a [0,1]
- tensor_basis(): Base HAL bivariada con escalado a [0,1]^2
- influence_curve(): Curva de influencia con Q_bar y Q_W
- simulate_data(), Q0_dgp(), g0_dgp(), omega_oracle(), true_drc(): DGP del paper

## Decisiones de diseño

Ver SPEC/decisions.md para el log completo.

Decisiones clave:
- D-0001: Solo R base, sin dependencias externas
- D-0008: Truncamiento de pesos w_cap como parametro explicito (default NULL)
- D-0009: Orden de suavidad k siempre explicito
- D-0010: Repo 100% nativo, sin oraculos de validacion externos

Preguntas resueltas:
- Q-01: Base HAL segun van der Laan (2023), terminos de borde = potencias, escalado a [0,1]
- Q-04: Densidad condicional nativa (KDE por bins)
- Q-05: Sin Sigma^-1 por defecto
- Q-17: Nudos estrictamente interiores para evitar colinealidad

## Tests

cd ~/repos/t-hal-mle-native
Rscript -e 'testthat::test_local(".")'

53 tests en verde.

## Referencias

- Paper: arXiv:2607.03824v1 (T-HAL-MLE)
- Base HAL: van der Laan (2023), arXiv:2301.13354
- Implementacion de referencia: https://github.com/amurlaniakea/thalmle

## Licencia

AGPL-3.0-or-later

## Autor

Pedro Sordo Martínez, amurlaniakea@gmail.com
