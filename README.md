# thalnative

Implementacion nativa en R base del Targeted Highly Adaptive Lasso Minimum Loss Estimation (T-HAL-MLE) para curvas dosis-respuesta continuas.

Estado: v0.0.1, 48 tests pasando, 17 funciones exportadas, 100% R base (sin dependencias externas).

## Caracteristicas

- Cero dependencias externas: sin glmnet, hal9001, haldensify, ni Matrix. Solo R base + stats.
- LASSO nativo: Coordinate Descent con warm starts, validado contra glmnet (error max 7e-05).
- Base HAL matematicamente fiel: terminos de borde segun van der Laan (2023), nudos en min(A) para evitar sesgo.
- Nuisances nativos: KDE para omega, HAL bivariado con CV para Q, KDE por bins para g.
- Inferencia completa: curva de influencia + IC de Wald puntuales.

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
  w_cap = "obs_max",
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

- thal_fit(): Pipeline completo (ajusta nuisances + targeting + inferencia)
- thal_core_native(): Nucleo del targeting (Algorithm 1)
- fit_Q_hal_native(): HAL bivariado del outcome con CV
- fit_g_native(): Densidad condicional g(a|w) por KDE en bins
- fit_omega_kde_native(): Densidad marginal omega(a) por KDE Sheather-Jones
- lasso_path(): LASSO nativo (Coordinate Descent)
- spline_basis(): Base HAL univariada (k=0, k=1)
- tensor_basis(): Base HAL bivariada (producto tensorial)
- influence_curve(): Curva de influencia para inferencia

## Decisiones de diseno

Ver SPEC/decisions.md para el log completo de decisiones (D-xxxx) y preguntas resueltas/abiertas (Q-xx).

Decisiones clave:
- D-0001: Solo R base, sin dependencias externas
- D-0008: Truncamiento de pesos w_cap como parametro explicito (default "obs_max")
- D-0009: Orden de suavidad k siempre explicito, sin variables de entorno ocultas
- D-0010: Repo 100% nativo, sin oraculos de validacion externos

Preguntas resueltas:
- Q-01: Base HAL segun van der Laan (2023), terminos de borde = potencias
- Q-04: Densidad condicional nativa (KDE por bins), sin haldensify
- Q-05: Sin Sigma^-1 por defecto, covariable inteligente = omega/g * phi
- Q-17: Nudos incluyen min(A) para evitar hueco en borde izquierdo

## Tests

cd ~/repos/t-hal-mle-native
Rscript -e 'testthat::test_local(".")'

48 tests en verde:
- basis: 10 (base univariada, terminos de borde, nudos)
- core: 7 (targeting con oraculos, DGP 1/2, tratamiento normal)
- lasso: 2 (Coordinate Descent, convergencia)
- omega: 4 (KDE marginal, uniforme, normal truncada)
- Q: 8 (HAL bivariado, CV, recuperacion de Q0)
- g: 5 (KDE por bins, independencia A-W, fuera de muestra)
- pipeline: 12 (thal_fit end-to-end, DGP 1/2, IC de Wald)

## Referencias

- Paper: arXiv:2607.03824v1 (T-HAL-MLE)
- Base HAL: van der Laan (2023), arXiv:2301.13354
- Implementacion de referencia: https://github.com/amurlaniakea/thalmle

## Licencia

AGPL-3.0-or-later

## Autor

Pedro Sordo Martinez, amurlaniakea@gmail.com
