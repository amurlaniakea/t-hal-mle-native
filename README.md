# thalnative

Implementacion nativa en R base del Targeted Highly Adaptive Lasso Minimum Loss Estimation (T-HAL-MLE) para curvas dosis-respuesta continuas.

Estado: v0.0.1 EXPERIMENTAL. 80 expectativas (27 bloques de test) pasando, 20 funciones exportadas, 100% R base. Prototipo de investigacion: la inferencia aun no esta validada y varias piezas difieren del paper (ver SPEC/decisions.md).

## Caracteristicas

- Cero dependencias externas: sin glmnet, hal9001, haldensify, ni Matrix. Solo R base + stats.
- LASSO nativo: Coordinate Descent con warm starts y mantenimiento de residuo.
- Base HAL segun van der Laan (2023): terminos de borde, nudos estrictamente interiores y escalado a [0,1] (sin colinealidad; rango completo verificado). No comparada con hal9001.
- Nuisances nativos: KDE para omega, HAL bivariado con CV para Q, KDE por bins para g (g_n NO es el estimador HAL de densidad del paper).
- Inferencia (experimental): curva de influencia del Lemma 3.1 (componentes Q_bar y Q_W, con Sigma^+) + IC de Wald puntuales. En el piloto con los DGP exactos del paper (n = 500, B = 50, nuisances estimadas) la cobertura de Wald fue de 0.72 a 0.92 segun el escenario (nominal 0.95): aun no es fiable. Causas candidatas, no probadas: seleccion LASSO data-adaptativa tratada como sieve fijo y residuos dentro de muestra (el paper hace cross-fitting K = 2 para la varianza). Ver SPEC/decisions.md (D-0004).

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
RMSE medio por dosis (media de los 25 RMSE por dosis, dosis en [1, 9]); "T-HAL mejor" = replicas de 50 con menor RMSE que el plug-in HAL.

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

### Topes de peso (exploratorio)

n = 500, B = 50, tratamiento normal, g y omega estimadas
(`experiments/results/wcap_exploratorio_n500_B50.csv`). RMSE agregado sobre las 25 dosis y todas las replicas (no comparable en nivel con el RMSE medio por dosis de la tabla del piloto; las diferencias frente al HAL si lo son):

| DGP | plug-in HAL | T-HAL tope por defecto | obs_max | tope 10 | tope 5 |
|---|---|---|---|---|---|
| 1 | 1.29 | 5.43 | 1.23 | 1.49 | 1.05 |
| 2 | 1.23 | 4.32 | 1.07 | 1.22 | 0.99 |

Con el tope del paper (sqrt(n) log(n/5), aprox. 103) T-HAL se dispara en las dosis extremas; con topes mas
estrictos mejora al plug-in HAL. Con g y omega oraculo el resultado es parecido, asi que la causa no es la
estimacion de g. Los topes se compararon en la misma simulacion: es evidencia exploratoria; la confirmacion con semillas nuevas esta mas abajo.

### Topes de peso (confirmacion con semilla nueva)

n = 500, B = 60, tratamiento normal, g y omega estimadas, semilla base 777001
(`experiments/results/wcap_confirmacion_n500_B60_s777001.csv`). Los topes (3, 5, 8, `obs_max` y el del paper)
se fijaron antes de correr, con semillas distintas a las del experimento exploratorio.
RMSE agregado sobre las 25 dosis y todas las replicas (no comparable en nivel con el RMSE medio por dosis de la tabla del piloto; las diferencias frente al HAL si lo son):

| DGP | plug-in HAL | T-HAL tope paper | tope 3 | tope 5 | tope 8 | obs_max |
|---|---|---|---|---|---|---|
| 1 | 1.28 | 6.97 | 0.98 | 1.19 | 1.54 | 1.32 |
| 2 | 1.24 | 5.32 | 0.92 | 1.04 | 1.22 | 1.13 |
| 3 | 1.72 | 6.53 | 1.49 | 1.67 | 1.94 | 1.72 |

Diferencia de RMSE por replica frente al plug-in HAL (negativa = T-HAL mejor), con su error estandar entre
parentesis, y replicas (de 60) en las que T-HAL gana:

| tope | DGP 1 | DGP 2 | DGP 3 |
|---|---|---|---|
| 3 | -0.36 (0.037), 55 | -0.37 (0.029), 58 | -0.25 (0.096), 40 |
| 5 | -0.22 (0.063), 45 | -0.27 (0.041), 48 | -0.11 (0.110), 37 |
| 8 | +0.01 (0.100), 42 | -0.13 (0.062), 47 | +0.12 (0.129), 31 |
| obs_max | -0.16 (0.082), 46 | -0.20 (0.052), 46 | -0.05 (0.109), 32 |

Lectura (n = 500 y tratamiento normal unicamente):

- En los tres DGP el orden es el mismo: cuanto mas bajo el tope, mejor (3 < 5 < 8 en RMSE). El tope 3 es el mas
  bajo probado, asi que el optimo no esta acotado: puede estar por debajo.
- La ganancia es clara en los DGP 1 y 2; en el DGP 3 es de unos 2.7 errores estandar con el tope 3 y el tope 5
  apenas se distingue del plug-in HAL.
- La mejora viene de las dosis centrales (RMSE ~0.43-0.58 frente a 0.80-1.51 del plug-in HAL). En los extremos la
  mejora es modesta (DGP 1 y 2) o nula: en el DGP 3 con tope 3 el RMSE en extremos (2.15) y en a = 9 (4.26) es
  peor que el del plug-in HAL (1.95 y 2.97).
- Con tratamiento uniforme omega/g = 1 y el tope no actua. Con tratamiento normal, E[omega/g] = 1, asi que un
  tope de 3 equivale a truncar a 3 veces el peso medio (el del paper, ~103, a ~100 veces).
- Aun no hay evidencia de cobertura de Wald con tope bajo, ni de n = 1000 o 2000. Por eso el tope por defecto
  del paquete sigue siendo el del paper; `w_cap` es una opcion explicita.

### Cobertura con topes bajos (semilla nueva)

n = 500, B = 60, tratamiento normal, semilla base 888001, nuisances estimadas
(`experiments/results/cobertura_topes_n500_B60_s888001.csv`; `mc_paper.R ... 888001 <tope>`).
RMSE medio por dosis, como en la tabla del piloto:

| DGP | tope | RMSE HAL | RMSE T-HAL | dif. frente a HAL (EE) | cob. Wald | cob. MC-SE | ancho Wald | ancho MC-SE |
|---|---|---|---|---|---|---|---|---|
| 1 | 3 | 1.10 | 0.81 | -0.34 (0.059) | 0.82 | 0.94 | 1.95 | 3.07 |
| 1 | 2 | 1.10 | 0.74 | -0.44 (0.042) | 0.81 | 0.94 | 1.86 | 2.81 |
| 2 | 3 | 1.12 | 0.82 | -0.39 (0.035) | 0.75 | 0.92 | 1.94 | 2.97 |
| 2 | 2 | 1.12 | 0.79 | -0.42 (0.029) | 0.73 | 0.92 | 1.82 | 2.86 |
| 3 | 3 | 1.64 | 1.15 | -0.19 (0.085) | 0.87 | 0.93 | 3.62 | 4.33 |
| 3 | 2 | 1.64 | 1.08 | -0.31 (0.079) | 0.87 | 0.93 | 3.37 | 4.03 |

Con el tope del paper (piloto, semilla distinta) la cobertura de Wald fue 0.77 / 0.72 / 0.81 y el ancho MC-SE
9.1 / 7.9 / 11.9. Lectura (n = 500 y tratamiento normal; semillas y B distintos entre piloto y esta tanda):

- Un tope bajo reduce el ancho MC-SE a un tercio y mantiene la ventaja sobre el plug-in HAL, pero la cobertura de
  Wald solo sube a 0.73-0.87: el IC de Wald sigue siendo demasiado estrecho (ancho Wald / ancho MC-SE de 0.64-0.66
  en los DGP 1 y 2 y 0.84 en el DGP 3).
- El tope no explica por si solo el defecto: con tratamiento uniforme (pesos iguales a 1) el piloto ya daba una
  cobertura de Wald de 0.81 / 0.81 / 0.92.
- Por dosis (`experiments/diag_cobertura.R`): el cociente SE medio / sd Monte Carlo es menor que 1 en casi todas las
  dosis de los DGP 1 y 2 (0.44-0.93 con tope 2) y mas bajo en a = 9 en los tres DGP (0.33-0.56). En el DGP 3 el SE es
  adecuado o conservador en las dosis 1-7 y la infracobertura se concentra en a = 8-9. En el DGP 2 el sesgo es de
  0.5-0.7 desviaciones tipicas en varias dosis (a = 1, 2, 3, 7), de modo que el sesgo tambien contribuye.
- Las causas de D-0004 (LASSO data-adaptativo tratado como fijo, residuos dentro de muestra) siguen sin
  probarse; este experimento no las distingue de la estimacion de g y omega.

## Tests

cd ~/repos/t-hal-mle-native
Rscript -e 'testthat::test_local(".")'

80 expectativas (27 bloques `test_that`) en verde al 2026-10-13.

## Referencias

- Paper: arXiv:2607.03824v1 (T-HAL-MLE)
- Base HAL: van der Laan (2023), arXiv:2301.13354
- Implementacion de referencia: https://github.com/amurlaniakea/thalmle

## Licencia

AGPL-3.0-or-later

## Autor

Pedro Sordo Martínez, amurlaniakea@gmail.com
