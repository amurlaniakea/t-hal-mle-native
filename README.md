# thalnative

Implementacion nativa en R base del Targeted Highly Adaptive Lasso Minimum Loss Estimation (T-HAL-MLE) para curvas dosis-respuesta continuas.

Estado: v0.0.1 EXPERIMENTAL. 95 expectativas (33 bloques de test) pasando, 20 funciones exportadas, 100% R base. Prototipo de investigacion: la inferencia aun no esta validada y varias piezas difieren del paper (ver SPEC/decisions.md).

## Caracteristicas

- Cero dependencias externas: sin glmnet, hal9001, haldensify, ni Matrix. Solo R base + stats.
- LASSO nativo: Coordinate Descent con warm starts y mantenimiento de residuo.
- Base HAL segun van der Laan (2023): terminos de borde, nudos estrictamente interiores y escalado a [0,1] (sin colinealidad; rango completo verificado). No comparada con hal9001.
- Nuisances nativos: KDE para omega, HAL bivariado con CV para Q, KDE por bins para g (g_n NO es el estimador HAL de densidad del paper).
- Inferencia (experimental): curva de influencia del Lemma 3.1 (componentes Q_bar y Q_W, con Sigma^+) + IC de Wald puntuales. En el piloto con los DGP exactos del paper (n = 500, B = 50, nuisances estimadas) la cobertura de Wald fue de 0.72 a 0.92 segun el escenario (nominal 0.95): aun no es fiable. Causa identificada en parte: con `k_Q = 0` y pocos nudos en el HAL de Q (`n_knots_Q = 10`, por defecto) la varianza real de `psi_hat` tiene picos en las dosis cercanas a los nudos que el IC no recoge; con `n_knots_Q = 20` o `k_Q = 1` desaparecen (ver "Hallazgo principal"). Siguen sin explicar las dosis extremas con tratamiento normal. Ver SPEC/decisions.md (D-0004, D-0012).

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

## Hallazgo principal: resolucion del HAL de Q

Con `k_Q = 0` (escalones) y `n_knots_Q = 10`, los nudos del HAL de Q (cuantiles de la muestra) son pocos y cambian de
una muestra a otra. La funcion estimada cambia de nivel en esos nudos y la varianza de `psi_hat` crece cerca de ellos en
proporcion a la pendiente al cuadrado de la verdad (Spearman 0.70-0.90 entre `var_MC` y `psi0'(a)^2` en seis corridas
uniformes de los DGP 1 y 2; ~0 en el DGP 3, donde `k_Q = 1`). El IC trata la base como fija y no recoge ese termino
(Spearman del SE^2 con la pendiente: 0.05-0.58). Dos arreglos lo eliminan:

DGP 1, tratamiento uniforme, ventana `c1 = 9, c2 = 14`, cross-fitting K = 2 (RMSE medio por dosis):

| configuracion | var_MC en a = 4 | ancho MC-SE | cob. Wald en muestra | SE/sd en muestra | RMSE HAL / T-HAL | coste por replica |
|---|---|---|---|---|---|---|
| `k_Q = 0`, 10 nudos (B = 60, s555001) | 1.84 | 2.85 | 0.797 | 0.87 | 1.13 / 0.747 | ~55 s |
| `k_Q = 0`, 20 nudos (B = 30, s888001) | 0.18 | 1.59 | 0.956 | 1.10 | 0.767 / 0.426 | ~168 s |
| `k_Q = 1`, 10 nudos (B = 60, s999001) | 0.03 | 1.85 | 0.968 | 1.52 | 1.10 / 0.490 | ~47 s |

(`experiments/results/nkq20_*`, `kq1_*` y `sieve_*`.) Con `n_knots_Q = 20` en otros escenarios (B = 30, semilla 123001,
`experiments/results/nkq20_cf_dgp123_*`):

| escenario | dif. RMSE frente a HAL (EE) | T-HAL gana en | cob. Wald en muestra -> CF | SE/sd en muestra -> CF | cob. MC-SE |
|---|---|---|---|---|---|
| DGP 2, uniforme, ventana 9-14 | -0.269 (0.017) | 30 de 30 | 0.927 -> 0.976 | 1.01 -> 1.38 | 0.940 |
| DGP 1, normal, tope 3 | -0.266 (0.035) | 27 de 30 | 0.799 -> 0.881 | 0.74 -> 1.00 | 0.932 |
| DGP 3, normal, tope 3 | +0.255 (0.092) | 10 de 30 | 0.889 -> 0.908 | 1.13 -> 1.21 | 0.936 |

Lectura (n = 500; B = 30-60; semillas distintas entre configuraciones, comparaciones no pareadas):

- Con un Q mas fino (o `k_Q = 1`) desaparecen los picos de varianza, el ancho del IC ideal se reduce hasta un 44 % y la
  cobertura de Wald con el SE en muestra sube a 0.93-0.97 en los escenarios uniformes. Tambien mejoran el HAL y T-HAL.
- Con `k_Q = 1` el SE queda sobrestimado (SE/sd ~1.5): el HAL con `k = 1` es sesgado para una verdad con escalones en w y
  los residuos incluyen ese sesgo. Con 20 nudos el SE en muestra queda calibrado (SE/sd 1.0-1.1).
- El cross-fitting de la varianza ya no hace falta con tratamiento uniforme y Q fino (sobrecorrige: SE/sd 1.4), pero con
  tratamiento normal en el DGP 1 corrige un SE corto (0.74 -> 1.00). No es una mejora universal.
- Con la media por dosis y la diferencia por replica T-HAL no gana en el DGP 3 con tratamiento normal (+0.255, 10 de 30);
  con la metrica del paper (mediana del RMSE por dosis) si queda por debajo del HAL. En a = 9 la varianza real de
  `psi_hat` (24.4) es ocho veces el SE^2 medio (2.95): las dosis extremas con pesos grandes estan sin resolver. Ver
  "Comparacion con las metricas del paper" y "DGP 3 con tratamiento normal: donde pierde T-HAL y por que".
- El valor por defecto sigue siendo `n_knots_Q = 10`: pasar a 20 cuesta ~3 veces mas de tiempo y falta confirmarlo con
  otros n y mas replicas. Si el HAL del paper usa muchos mas nudos (p. ej. los puntos observados), esta es una
  diferencia relevante de esta implementacion; no esta verificado contra el paper.

## Comparacion con las metricas del paper (seccion 7.3)

El paper (arXiv:2607.03824, secciones 7.1-7.3) usa la misma sd = 2 para A | W (truncada a [0, 10]), la misma rejilla de 25
dosis en [1, 9] y el mismo selector de tamano (6 y 9), y resume el RMSE y el sesgo PUNTUALES de cada dosis con diagramas de
caja sobre las 25 dosis. `experiments/resumen_tipo_paper.R` calcula esos resumenes (mediana y cuartiles del RMSE por dosis,
proporcion de dosis con |sesgo / MC-SE| <= 1/log(n), cobertura con el MC-SE y cobertura de Wald). n = 500; HAL-MLE,
T-HAL-MLE (Plugin) y T-HAL-MLE son `psi_hal`, `psi_plugin` y `psi_hat`.

RMSE puntual por dosis: mediana [Q1 - Q3] sobre las 25 dosis (salida completa en
`experiments/results/resumen_tipo_paper_compacto.txt`):

| escenario | HAL | T-HAL (Plugin) | T-HAL | RMSE max. HAL / T-HAL | cob. MC-SE min. (T-HAL) | cob. Wald mediana (min.) |
|---|---|---|---|---|---|---|
| DGP 1 uniforme, 10 nudos (B = 60) | 1.15 [0.71 - 1.61] | 0.84 [0.44 - 1.49] | 0.66 [0.46 - 0.82] | 2.29 / 1.01 | 0.73 | 0.88 (0.53) |
| DGP 2 uniforme, 10 nudos (B = 60) | 1.00 [0.60 - 1.40] | 0.74 [0.53 - 1.38] | 0.63 [0.48 - 0.79] | 2.30 / 0.95 | 0.85 | 0.93 (0.70) |
| DGP 3 uniforme, 10 nudos (B = 60) | 1.11 [0.86 - 1.24] | 0.60 [0.45 - 0.96] | 0.64 [0.46 - 0.94] | 1.66 / 1.45 | 0.63 | 0.98 (0.72) |
| DGP 1 normal, 10 nudos (B = 60) | 0.89 [0.62 - 1.30] | 0.96 [0.51 - 1.20] | 0.51 [0.46 - 0.91] | 2.99 / 2.19 | 0.88 | 0.87 (0.67) |
| DGP 2 normal, 10 nudos (B = 60) | 1.00 [0.65 - 1.20] | 1.01 [0.53 - 1.14] | 0.64 [0.41 - 1.08] | 2.69 / 2.08 | 0.83 | 0.85 (0.60) |
| DGP 3 normal, 10 nudos (B = 60) | 0.87 [0.64 - 1.18] | 0.71 [0.43 - 1.34] | 0.77 [0.47 - 1.34] | 2.08 / 3.71 | 0.88 | 0.92 (0.68) |
| DGP 2 uniforme, 20 nudos (B = 30) | 0.60 [0.44 - 0.94] | 0.54 [0.40 - 0.94] | 0.42 [0.35 - 0.50] | 1.35 / 0.97 | 0.87 | 1.00 (0.90) |
| DGP 1 normal, 20 nudos (B = 30) | 0.72 [0.50 - 0.97] | 0.68 [0.46 - 0.98] | 0.37 [0.32 - 0.93] | 2.32 / 1.41 | 0.87 | 0.90 (0.77) |
| DGP 3 normal, 20 nudos (B = 30) | 0.81 [0.59 - 1.18] | 0.60 [0.34 - 1.25] | 0.60 [0.45 - 1.23] | 2.76 / 4.90 | 0.77 | 0.93 (0.43) |

Lectura (n = 500, B = 30-60, semillas distintas entre filas):

- Con las metricas del paper se reproduce su ordenacion: con tratamiento uniforme T-HAL < T-HAL (Plugin) < HAL en los DGP 1 y 2;
  en el DGP 3 el HAL tiene la mediana mas alta y T-HAL y su plug-in son comparables (el paper dice lo mismo). Con
  tratamiento normal T-HAL tiene la mediana mas baja en los tres DGP, incluido el DGP 3 (0.77 frente a 0.87 con 10 nudos;
  0.60 frente a 0.81 con 20).
- Las perdidas de T-HAL en el DGP 3 normal que se describen mas abajo aparecen al usar la media por dosis o la diferencia por
  replica, que pesan las dosis extremas: el RMSE maximo por dosis de T-HAL (3.7 - 4.9, en a = 9) supera al del HAL (2.1 - 2.8).
- Con 10 nudos en el HAL de Q, el HAL y T-HAL tienen RMSE ~1.5 veces mayores que en la Figura 1 del paper (lectura a ojo del
  grafico, no una tabla: ~0.75 / 0.55 / 0.4 para HAL / Plugin / T-HAL en el DGP 1 uniforme); con 20 nudos la escala es parecida
  (media por dosis 0.77 y 0.43 en ese escenario).
- La cobertura Wald mediana de T-HAL con tratamiento uniforme es 0.88-0.98 y con tratamiento normal 0.85-0.93; el paper
  da ~0.95 y "considerable infracobertura", respectivamente. La cobertura MC-SE minima por dosis de T-HAL (0.63 - 0.88) esta
  dentro del ruido con B = 30-60 salvo las dosis ya descritas; la del HAL llega a 0.07 - 0.13 en el DGP 3, algo que el paper no
  muestra: nuestro HAL con pocos nudos es sesgado en alguna dosis.
- Limites: n = 500 y B = 30-60 (el paper usa B = 1000 y n hasta 2000); no se han comparado cifra a cifra con tablas del paper.

## DGP 3 con tratamiento normal: donde pierde T-HAL y por que

**Donde.** Dosis a dosis (n = 500, B = 60, tope 3, `n_knots_Q = 10`, g y omega estimadas; `z` = diferencia de MSE entre
T-HAL y HAL dividida por su error estandar, positivo = T-HAL peor; `experiments/diag_rmse_dosis.R`):

| dosis | 1 | 2.33 | 3.67 | 5 | 6.33 | 7.67 | 9 |
|---|---|---|---|---|---|---|---|
| RMSE HAL | 1.51 | 0.50 | 0.66 | 0.68 | 0.57 | 1.12 | 2.08 |
| RMSE T-HAL | 1.34 | 0.68 | 0.28 | 0.35 | 0.77 | 1.75 | 3.71 |
| z | -1.4 | +3.0 | -7.7 | -5.3 | +2.1 | +3.4 | +3.1 |

T-HAL gana con claridad en el centro y en el extremo inferior (el HAL tiene alli un sesgo de -1.2 a -1.4) y pierde en las
dosis altas: las tres mas altas aportan ~97 % del exceso de MSE. En a = 9 T-HAL reduce el sesgo (-0.95 -> +0.48) pero
duplica la desviacion tipica (1.87 -> 3.71). Con g y omega oraculo el resultado por dosis es casi identico.

**Perillas que no lo arreglan** (comparaciones sobre los mismos datos dentro de cada fila):

| palanca | valores | sd T-HAL en a = 9 (HAL) | otro efecto |
|---|---|---|---|
| tope de pesos (semilla 444001) | 2 / 1.5 / 1 | 3.48 / 3.25 / 3.09 (2.30) | RMSE medio por dosis 1.009 / 0.974 / 0.938 (HAL 1.034); cob. Wald CF 0.903 / 0.889 / 0.867 |
| tamano del targeting (semilla 515001) | 7 / 11 / 21 funciones | 3.00 / 3.12 / 5.00 (2.20) | RMSE en a = 5: 1.26 / 0.67 / 0.36 (HAL 0.69); cob. Wald en muestra 0.51 / 0.68 / 0.86 |

Bajar el tope apenas reduce la varianza del extremo; reducir el modelo de targeting la reduce pero devuelve el sesgo en el
centro y hunde la cobertura.

**Por que: falta de solapamiento.** Con A | W ~ N(5 + 1.5 W, sd^2) truncada a [0, 10], en a = 9 y para w < -0.5 la
probabilidad condicional es P(A > 9 | W) ~ 0.003 con sd = 2; el 37 % de la poblacion tiene alli un peso verdadero mayor
que 3. Si se ensancha A | W (semilla 515001, tope 3, B = 60, mismos W y ruido):

| sd de A \| W | P(A > 9 \| W < -0.5) | dif. RMSE por replica frente a HAL (EE) | T-HAL gana en | RMSE medio por dosis HAL / T-HAL | sd en a = 9 HAL / T-HAL | cob. Wald en muestra -> CF |
|---|---|---|---|---|---|---|
| 2 | 0.0027 | +0.395 (0.086) | 12 de 60 | 1.06 / 1.23 | 2.20 / 5.00 | 0.861 -> 0.883 |
| 3 | 0.0201 | +0.019 (0.046) | 33 de 60 | 1.00 / 0.92 | 1.52 / 2.58 | 0.916 -> 0.939 |
| 4 | 0.0410 | -0.114 (0.037) | 38 de 60 | 1.00 / 0.82 | 1.39 / 2.10 | 0.922 -> 0.941 |

Al mejorar el solapamiento, T-HAL pasa de perder a ganar, el cociente de desviaciones tipicas en a = 9 baja de 2.3 a 1.5 y
la cobertura de Wald con cross-fitting llega a ~0.94. Con tratamiento uniforme (solapamiento total) T-HAL gana en todos los
DGP. Lectura (un DGP, n = 500, B = 60):

- La perdida del DGP 3 normal con sd = 2 se explica en gran parte por la falta de solapamiento en las dosis altas, no por
  los topes ni por el tamano del targeting ni por estimar g y omega.
- Incluso con sd = 4, T-HAL sigue algo peor que el HAL en a = 9 (RMSE 2.09 frente a 1.68, z = +2.2).
- El paper confirma sd = 2 para A | W, el truncamiento a [0, 10] y la rejilla de 25 dosis; con sd = 2 los pesos verdaderos
  llegan a 29 en a = 9 (w = -1.5), aunque el paper dice que sus DGP no sufren pesos grandes (usa haldensify para g, que no se
  ha comparado aqui; con g oraculo no cambia nada).
- La comparacion global depende de la metrica: con tope <= 2 (semilla 444001) el RMSE medio por dosis da a T-HAL por
  debajo del HAL, mientras que la diferencia por replica (RMS entre dosis) lo da por encima. La metrica del paper (mediana
  del RMSE por dosis) lo da por debajo (ver "Comparacion con las metricas del paper").

## Evidencia: Monte Carlo piloto con los DGP del paper

> **Aviso (2026-10-17).** Las tablas de esta seccion usan `n_knots_Q = 10` (valor por defecto) salvo que se indique lo
> contrario. Esa resolucion del HAL de Q es insuficiente en los DGP 1 y 2 (`k_Q = 0`): infla la varianza de `psi_hat`
> en ciertas dosis y la cobertura de Wald medida (0.72-0.89) esta afectada (ver "Hallazgo principal"). Las conclusiones
> sobre T-HAL frente al HAL en los DGP 1 y 2 se mantienen; las cifras de cobertura y ancho no deben tomarse como el
> rendimiento del metodo con una resolucion adecuada.

> **Aviso (2026-10-15).** Hasta el commit `f4d86a5`, `Q_fun` escalaba la dosis con el rango de los puntos donde se
> evaluaba y no con el del entrenamiento. Solo afecta a `k_Q = 1` (DGP 3): el plug-in HAL de referencia (`psi_hal`) y
> `psi_plugin` evaluados en la rejilla salian mal en el DGP 3, asi que **las comparaciones T-HAL frente a HAL del DGP 3
> de las tablas de esta seccion (piloto, topes y cobertura) no son validas** (la comparacion corregida esta en "DGP 3 con el plug-in HAL corregido"). No afecta a los
> DGP 1 y 2 (`k_Q = 0`), ni a `psi_hat` ni a su error estandar.

n = 500, B = 50 replicas por escenario, DGP exactos de la seccion 7.1 del paper
(`experiments/mc_paper.R`; resumen en `experiments/results/piloto_n500_B50.csv`).
RMSE medio por dosis (media de los 25 RMSE por dosis, dosis en [1, 9]); "T-HAL mejor" = replicas de 50 con menor RMSE que el plug-in HAL.

| DGP | tratamiento | RMSE HAL | RMSE T-HAL | T-HAL mejor | cobertura Wald |
|---|---|---|---|---|---|
| 1 | uniforme | 1.17 | 0.67 | 50 | 0.81 |
| 2 | uniforme | 1.07 | 0.68 | 50 | 0.81 |
| 3 | uniforme | 1.28 (invalido) | 0.77 | 47 (invalido) | 0.92 |
| 1 | normal | 1.10 | 2.49 | 17 | 0.77 |
| 2 | normal | 1.06 | 2.14 | 18 | 0.72 |
| 3 | normal | 1.71 (invalido) | 3.10 | 10 (invalido) | 0.81 |

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
| 3 | 1.72 (HAL invalido) | 6.53 | 1.49 | 1.67 | 1.94 | 1.72 |

Diferencia de RMSE por replica frente al plug-in HAL (negativa = T-HAL mejor), con su error estandar entre
parentesis, y replicas (de 60) en las que T-HAL gana:

| tope | DGP 1 | DGP 2 | DGP 3 (dif. frente a HAL invalida) |
|---|---|---|---|
| 3 | -0.36 (0.037), 55 | -0.37 (0.029), 58 | -0.25 (0.096), 40 |
| 5 | -0.22 (0.063), 45 | -0.27 (0.041), 48 | -0.11 (0.110), 37 |
| 8 | +0.01 (0.100), 42 | -0.13 (0.062), 47 | +0.12 (0.129), 31 |
| obs_max | -0.16 (0.082), 46 | -0.20 (0.052), 46 | -0.05 (0.109), 32 |

Lectura (n = 500 y tratamiento normal unicamente):

- En los tres DGP el orden es el mismo: cuanto mas bajo el tope, mejor (3 < 5 < 8 en RMSE). El tope 3 es el mas
  bajo probado, asi que el optimo no esta acotado: puede estar por debajo.
- La ganancia frente al plug-in HAL es clara en los DGP 1 y 2. La del DGP 3 de esta tabla no es valida (ver aviso);
  el orden entre topes si se mantiene, porque todas las variantes comparten el mismo HAL.
- En los DGP 1 y 2 la mejora viene de las dosis centrales (RMSE ~0.43-0.49 frente a 0.80-0.83 del plug-in HAL) y en
  los extremos es modesta. Las cifras del DGP 3 frente al HAL de esta tabla no son validas (ver aviso).
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
| 3 | 3 | n.v. | 1.15 | n.v. | 0.87 | 0.93 | 3.62 | 4.33 |
| 3 | 2 | n.v. | 1.08 | n.v. | 0.87 | 0.93 | 3.37 | 4.03 |

Con el tope del paper (piloto, semilla distinta) la cobertura de Wald fue 0.77 / 0.72 / 0.81 y el ancho MC-SE
9.1 / 7.9 / 11.9. Lectura (n = 500 y tratamiento normal; semillas y B distintos entre piloto y esta tanda):

- Un tope bajo reduce el ancho MC-SE a un tercio y mantiene la ventaja sobre el plug-in HAL en los DGP 1 y 2, pero la cobertura de
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

### Cross-fitting K = 2 para la varianza (pareado)

n = 500, B = 60, tratamiento normal, tope 3, semilla base 222001
(`experiments/results/cf_pareado_n500_B60_s222001.csv`; `mc_cf.R`, resumen con `mc_cf_summary.R`). Cada replica
produce `psi_hat` y DOS errores estandar del mismo ajuste: con residuos en muestra y con cross-fitting K = 2
(`n_cf = 2`: residuos del Q inicial fuera de fold mas la correccion de targeting de la muestra completa). Es una
variante propia, no verificada contra el procedimiento exacto del paper. `psi_hat` no cambia; el coste es ~70 % mas
de tiempo por replica.

| DGP | g y omega | cob. Wald en muestra | cob. Wald CF | dif. pareada (EE) | SE/sd en muestra | SE/sd CF |
|---|---|---|---|---|---|---|
| 1 | estimadas | 0.801 | 0.861 | +0.060 (0.0074) | 0.73 | 0.88 |
| 1 | oraculo | 0.819 | 0.881 | +0.061 (0.0088) | 0.76 | 0.92 |
| 2 | estimadas | 0.763 | 0.852 | +0.089 (0.0103) | 0.74 | 1.00 |
| 2 | oraculo | 0.767 | 0.869 | +0.102 (0.0099) | 0.75 | 1.03 |
| 3 | estimadas | 0.877 | 0.897 | +0.021 (0.0040) | 1.00 | 1.09 |
| 3 | oraculo | 0.878 | 0.903 | +0.025 (0.0051) | 1.03 | 1.12 |

`SE/sd` = SE medio dividido por la desviacion tipica Monte Carlo de `psi_hat`, promediado sobre las 25 dosis. Con
g y omega estimadas, el ancho medio del IC de Wald con cross-fitting es 2.40 / 2.70 / 3.84 frente a un ancho MC-SE de
2.97 / 2.89 / 4.29. Lectura (n = 500, tratamiento normal, tope 3, un solo conjunto de semillas):

- El cross-fitting sube la cobertura de Wald de forma pareada (+0.02 a +0.10; unos 6-12 errores estandar en los
  DGP 1 y 2) y acerca el SE a la dispersion real: de 0.73-0.74 a 0.88-1.00 en los DGP 1 y 2. Es coherente con que los
  residuos en muestra del Q inicial subestimaban la varianza (causa candidata de D-0004).
- La cobertura sigue por debajo de 0.95 (0.85-0.90), incluso en el DGP 2 con SE/sd cercano a 1. El sesgo medio es
  de 0.23-0.32 desviaciones tipicas y quedan por separar la heterogeneidad por dosis y las colas no normales.
- Con g y omega oraculo (Q ajustada) el resultado casi no cambia (SE/sd +0.01 a +0.03; cobertura hasta +0.02): a este
  n, estimar g y omega no explica la subestimacion.
- El valor por defecto del paquete sigue siendo `n_cf = 1`; `n_cf = 2` es una opcion para la inferencia.

### DGP 3 con el plug-in HAL corregido

Mismas corridas del cross-fitting pareado (n = 500, B = 60, tratamiento normal, tope 3, semilla base 222001, g y
omega estimadas), ya con el paquete corregido (`f4d86a5`). RMSE medio por dosis:

| DGP | RMSE HAL | RMSE T-HAL | dif. frente a HAL (EE) | T-HAL mejor en | sesgo abs. HAL / T-HAL |
|---|---|---|---|---|---|
| 1 | 1.095 | 0.784 | -0.367 (0.036) | 56 de 60 | 0.52 / 0.18 |
| 2 | 1.081 | 0.775 | -0.347 (0.032) | 55 de 60 | 0.54 / 0.21 |
| 3 | 0.966 | 1.109 | +0.346 (0.068) | 16 de 60 | 0.61 / 0.18 |

Con g y omega oraculo el DGP 3 da 1.121, +0.347 (0.077) y 19 de 60. Lectura:

- En el DGP 3, con el HAL bien evaluado, T-HAL tiene menos sesgo (0.18 frente a 0.61) pero mas RMSE que el plug-in HAL
  y pierde en 44 de 60 replicas. El paper reporta que T-HAL gana en todos los DGP; aqui no se reproduce en el DGP 3
  con `k_Q = 1`, tope 3 y n = 500.
- En los DGP 1 y 2 el resultado se mantiene (T-HAL gana en 55-56 de 60).
- Falta medir T-HAL frente al HAL corregido con otros topes en el DGP 3 y explorar el HAL del resultado (nudos y
  orden) y n mayores antes de sacar conclusiones sobre la causa.

### Tratamiento uniforme (sin pesos) y tamano del working model

n = 500, B = 60, tratamiento uniforme (omega/g = 1, el tope no actua), g y omega estimadas, cross-fitting K = 2.
Ventana por defecto del sieve (tamano mediano 21; semilla 333001) y ventana mayor `c1 = 9, c2 = 14` (tamano mediano
32-33; semilla 555001). Las semillas y las muestras son distintas entre ambas tandas
(`experiments/results/uniforme_*_s333001.csv` y `sieve_*_s555001_c9-14.csv`).

| DGP | ventana | dif. RMSE frente a HAL (EE) | T-HAL gana en | cob. Wald en muestra -> CF | SE/sd CF | cob. MC-SE | sesgo/sd | ancho MC-SE |
|---|---|---|---|---|---|---|---|---|
| 1 | defecto | -0.614 (0.020) | 60 de 60 | 0.822 -> 0.863 | 0.97 | 0.917 | 0.46 | 2.28 |
| 1 | mayor | -0.446 (0.017) | 60 de 60 | 0.797 -> 0.859 | 1.02 | 0.952 | 0.25 | 2.85 |
| 2 | defecto | -0.541 (0.021) | 60 de 60 | 0.831 -> 0.897 | 1.05 | 0.929 | 0.34 | 2.32 |
| 2 | mayor | -0.366 (0.017) | 59 de 60 | 0.812 -> 0.887 | 1.11 | 0.947 | 0.22 | 2.89 |
| 3 | defecto | -0.283 (0.034) | 51 de 60 | 0.923 -> 0.951 | 1.49 | 0.897 | 0.48 | 2.59 |
| 3 | mayor | -0.270 (0.048) | 48 de 60 | 0.973 -> 0.983 | 1.72 | 0.949 | 0.28 | 3.05 |

Lectura (n = 500, tratamiento uniforme; semillas distintas entre ventanas):

- Sin pesos, T-HAL gana al plug-in HAL en los tres DGP, tambien en el DGP 3, como reporta el paper.
- El cross-fitting sube la cobertura de Wald de forma pareada tambien aqui (+0.01 a +0.07).
- Con la ventana por defecto la cobertura de MC-SE (con un SE ideal) es 0.90-0.93: el sesgo del working model
  (0.34-0.48 desviaciones tipicas) marca un techo. Con la ventana mayor el sesgo baja a ~0.2-0.3 y la cobertura de
  MC-SE llega a ~0.95, a costa de un IC mas ancho (+18 a +25 % de ancho MC-SE) y menos ventaja en RMSE.
- La cobertura de Wald con la ventana mayor no mejora en los DGP 1 y 2 (0.86-0.89) aunque el SE medio sea correcto
  (SE/sd 1.02-1.11) y el techo con SE ideal sea ~0.95: queda una subestimacion que no es de sesgo ni de tamano del
  sieve y que se investiga (variabilidad del SE entre replicas, dosis extremas).
- En el DGP 3 el SE esta sobrestimado (SE/sd 1.49 y 1.72): el 0.95-0.98 de cobertura de Wald viene de un IC
  demasiado ancho y no indica un buen IC.
- Un tamano de sieve mayor es un parametro (`c1`, `c2` de `thal_fit`), no el valor por defecto.

### DGP 3 con tratamiento normal y topes bajos

n = 500, B = 60, semilla base 444001 (mismos datos en los tres topes), g y omega estimadas, cross-fitting K = 2,
HAL corregido (RMSE medio por dosis 1.03). `experiments/results/dgp3_topes_bajos_*_s444001.csv`.

| tope | RMSE T-HAL | dif. frente a HAL (EE) | T-HAL gana en | cob. Wald CF | SE/sd CF | cob. MC-SE | sesgo abs. T-HAL |
|---|---|---|---|---|---|---|---|
| 2 | 1.009 | +0.102 (0.063) | 28 de 60 | 0.903 | 1.12 | 0.931 | 0.157 |
| 1.5 | 0.974 | +0.046 (0.056) | 28 de 60 | 0.889 | 1.09 | 0.930 | 0.158 |
| 1 | 0.938 | -0.0004 (0.046) | 33 de 60 | 0.867 | 0.97 | 0.936 | 0.149 |

Con el tope 3 (otra semilla) la diferencia era +0.346 (EE 0.068). Lectura:

- Bajar el tope mejora a T-HAL de forma monotona y con tope 1 empata con el plug-in HAL, pero no lo supera. Con
  tratamiento uniforme T-HAL si gana (-0.28 en el DGP 3), asi que la varianza de los pesos explica buena parte de la
  perdida pero no toda.
- T-HAL conserva menos sesgo que el HAL (0.15-0.16 frente a 0.64) y pierde por varianza.
- RMSE y cobertura van en direcciones opuestas: con el tope mas bajo el SE baja mas deprisa que la dispersion real
  y la cobertura de Wald cae de 0.903 a 0.867. Ningun tope es el mejor para las dos cosas.
- Con la media por dosis y la diferencia por replica T-HAL no supera al HAL en el DGP 3 normal con `k_Q = 1` y n = 500; con
  la mediana del RMSE por dosis (metrica del paper) si.

## Tests

cd ~/repos/t-hal-mle-native
Rscript -e 'testthat::test_local(".")'

95 expectativas (33 bloques `test_that`) en verde al 2026-10-17.

## Referencias

- Paper: arXiv:2607.03824v1 (T-HAL-MLE)
- Base HAL: van der Laan (2023), arXiv:2301.13354
- Implementacion de referencia: https://github.com/amurlaniakea/thalmle

## Licencia

AGPL-3.0-or-later

## Autor

Pedro Sordo Martínez, amurlaniakea@gmail.com
