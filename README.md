# thalnative

Implementacion nativa en R base del Targeted Highly Adaptive Lasso Minimum Loss Estimation (T-HAL-MLE) para curvas dosis-respuesta continuas.

Estado: v0.0.1 EXPERIMENTAL. 62 expectativas (24 bloques de test) pasando, 20 funciones exportadas, 100% R base. Prototipo de investigacion: la inferencia aun no esta validada y varias piezas difieren del paper (ver SPEC/decisions.md).

## Caracteristicas

- Cero dependencias externas: sin glmnet, hal9001, haldensify, ni Matrix. Solo R base + stats.
- LASSO nativo: Coordinate Descent con warm starts y mantenimiento de residuo.
- Base HAL segun van der Laan (2023): terminos de borde, nudos estrictamente interiores y escalado a [0,1] (sin colinealidad; rango completo verificado). No comparada con hal9001.
- Nuisances nativos: KDE para omega, HAL bivariado con CV para Q, KDE por bins para g (g_n NO es el estimador HAL de densidad del paper).
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

Preguntas abiertas:
- Q-04: Densidad condicional. Se usa KDE por bins, no el estimador HAL del paper; sin validar contra la densidad verdadera

Preguntas resueltas:
- Q-01: Base HAL segun van der Laan (2023), terminos de borde = potencias, escalado a [0,1]
- Q-05: Sin Sigma^-1 por defecto
- Q-17: Nudos estrictamente interiores para evitar colinealidad

## Evidencia: Monte Carlo piloto con los DGP del paper

n = 500, B = 50 replicas por escenario, DGP exactos de la seccion 7.1 del paper
(`experiments/mc_paper.R`; resumen en `experiments/results/piloto_n500_B50.csv`).
RMSE medio sobre 25 dosis en [1, 9]; "T-HAL mejor" = replicas de 50 con menor RMSE que el plug-in HAL.

| DGP | tratamiento | RMSE HAL | RMSE T-HAL | T-HAL mejor | cobertura Wald |
|---|---|---|---|---|---|
| 1 | uniforme | 1.17 | 0.67 | 50 | 0.81 |
| 2 | uniforme | 1.07 | 0.68 | 50 | 0.81 |
| 3 | uniforme | 1.28 | 0.77 | 47 | 0.92 |
| 1 | normal | 1.10 | 2.49 | 17 | 0.77 |
| 2 | normal | 1.06 | 2.14 | 18 | 0.72 |
| 3 | normal | 1.71 | 3.10 | 10 | 0.81 |

Lectura (con las limitaciones de B = 50 y un solo n):

- Con tratamiento uniforme, T-HAL mejora al plug-in HAL en RMSE y sesgo, en linea con el paper.
- Con tratamiento normal (con confusion), T-HAL empeora: en el DGP 1 mejora en las dosis centrales,
  pero el error se dispara en los extremos del rango de dosis (RMSE 5.1 en a = 1 y 23.0 en a = 9, frente a
  1.3 y 2.1 del plug-in HAL). El paper reporta que T-HAL tambien gana con tratamiento normal; aqui no se reproduce.
- La cobertura de Wald queda por debajo del 0.95 nominal en todos los escenarios.

## Tests

cd ~/repos/t-hal-mle-native
Rscript -e 'testthat::test_local(".")'

62 expectativas (24 bloques `test_that`) en verde al 2026-10-11.

## Referencias

- Paper: arXiv:2607.03824v1 (T-HAL-MLE)
- Base HAL: van der Laan (2023), arXiv:2301.13354
- Implementacion de referencia: https://github.com/amurlaniakea/thalmle

## Licencia

AGPL-3.0-or-later

## Autor

Pedro Sordo Martínez, amurlaniakea@gmail.com
