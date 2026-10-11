# Decision log y preguntas abiertas
Copyright (c) Pedro Sordo Martinez. Licencia: AGPL-3.0-or-later.
Formato: Q-xx = pregunta (RESUELTA o abierta). D-xxxx = decision tomada.
Ultima actualizacion: 2026-10-11 (v0.0.1, 48 tests, 17 funciones exportadas).

## Decisiones

**D-0001 -- Dependencias de runtime**
Estado: tomada (confirmada).
Solo R base y stats. Sin Matrix, sin glmnet, sin hal9001, sin haldensify.
LASSO propio via Coordinate Descent con warm starts y mantenimiento de residuo.
Validado contra glmnet: diferencia maxima en coeficientes = 7.06e-05.

**D-0002 -- Oraculos de validacion (REVISADA por D-0010)**
Estado: eliminada.
Originalmente se permitia glmnet/hal9001/haldensify solo en validation/.
Ahora: no se usan en absoluto. Ver D-0010.

**D-0003 -- Alcance v0**
Estado: tomada.
Y continuo (perdida cuadratica), W y A univariados, k_targeting = 1.
Y binario y W multivariado quedan para versiones posteriores.

**D-0004 -- Inferencia**
Estado: tomada.
Inferencia primaria = curva de influencia (Lemma 3.1 / ec. (1)), IC de Wald puntuales.
Implementada en R/05-influence-curve.R: influence_curve(), variance_estimator(), wald_ci().
Bootstrap solo como diagnostico exploratorio.

**D-0005 -- Medida de integracion de A**
Estado: tomada.
Usar la medida empirica de A en Sigma, en la proyeccion final y en Q_W.
Coherente con el paso 7 del Algorithm 1.

**D-0006 -- Selector de lambda**
Estado: tomada.
Usar las desigualdades de seccion 5.1 (no las del Algorithm 1, que estan invertidas).
Implementado en select_lambda_by_size() con ventana [c1*n^(1/5), c2*n^(1/5)].

**D-0007 -- Afirmaciones tecnicas**
Estado: tomada.
Prohibido afirmar equivalencia con glmnet, validez de IC o uso clinico sin test citado.
Formato: "En el test X, con DGP Y, n = Z y R replicas, el error maximo fue E".

**D-0008 -- Truncamiento de pesos**
Estado: tomada.
w_cap es un parametro de primera clase en thal_fit() y thal_core_native().
Opciones: NULL (sqrt(n)*log(n/5) del paper), "obs_max" (maximo observado), numerico fijo.
Evidencia de thalmle: con tratamiento normal, obs_max y topes fijos (5, 10) mejoran
el RMSE; el default del paper es inestable en las colas.
Se informa siempre cual se uso.

**D-0009 -- Orden de suavidad**
Estado: tomada.
El orden k es siempre un argumento explicito, sin variable de entorno oculta.
k_Q = 0 para DGP 1/2 (discontinuidades), k_Q = 1 para DGP 3 (suave).
k_targeting = 1 por defecto.

**D-0010 -- Repo 100% nativo**
Estado: tomada.
No se usan dependencias externas ni siquiera como oraculos de validacion.
La fidelidad matematica se verifica contra van der Laan (2023) mediante tests unitarios.
Los nuisances se estiman con funciones nativas: KDE para omega, HAL propio para Q,
KDE por bins para g.

## Preguntas resueltas

**Q-01 (RESUELTA) -- Base HAL de orden k**
Fuente: van der Laan (2023), arXiv:2301.13354, Definicion 7 y ejemplo seccion 3.1.
k=0: 1, I(a>=u), I(w>=u), I(a>=u)I(w>=u).
k=1, d=1: 1, a (borde 0), (a-u)_+.
k=1, d=2: 1, a, w, a*w, (a-u)_+, (w-u)_+, w(a-u)_+, a(w-u)_+, (a-ua)_+(w-uw)_+.
Los terminos de borde son potencias (nudo en el origen), NO indicadores.
Implementado en R/01-basis.R (univariada) y R/07-basis-bivar.R (bivariada).

**Q-02 (RESUELTA) -- Nudos de la base del targeting**
50 nudos en cuantiles de A (N=51 con intercepto), como en thalmle.
make_knots() fuerza min(A) como primer nudo (ver Q-17).

**Q-03 (RESUELTA) -- HAL del outcome**
HAL bivariado nativo con base tensor_basis() y CV 5-fold (lambda.min).
k_Q y n_knots_Q son argumentos explicitos de fit_Q_hal_native().
Sin hal9001.

**Q-04 (RESUELTA) -- Densidad condicional g_n**
KDE Sheather-Jones por bins de W (n_bins = 5 por defecto).
Cada bin tiene su propia KDE 1D de A. Fallback a la marginal si < 5 datos en un bin.
Alternativa nativa a haldensify. Sin hazard por tramos ni penalizacion.
Implementado en R/09-nuisance-g.R.

**Q-05 (RESUELTA) -- Sigma inversa**
Sin Sigma^{-1} por defecto, como thalmle.
La covariable inteligente es omega/g * phi(A), no omega/g * Sigma^{-1} phi(A).
La "funcion seleccionada j" es phi_j para eps_j != 0.
Con Sigma^{-1} completa cada columna mezclaria todas las phi y la seleccion
perderia interpretabilidad.

**Q-09 (RESUELTA) -- Truncamiento de pesos**
Se trunca omega/g (no 1/g). Se aplica siempre.
Ver D-0008 para las opciones.

**Q-12 (RESUELTA) -- Criterios de convergencia del LASSO**
Coordinate Descent con warm starts, tolerancia 1e-6, max 200 iteraciones.
Mantenimiento del residuo para coste O(n*p) por barrido.
Validado contra glmnet: diferencia maxima 7.06e-05 en coeficientes.

**Q-13 (RESUELTA) -- smoothness_orders en DGP 1 y 2**
El README de thalmle confirma THAL_QSO=0 (orden 0).
En el repo nuevo es un argumento explicito (D-0009).

**Q-14 (RESUELTA) -- g_n evaluable fuera de muestra**
fit_g_native() devuelve una funcion general de (a, w) que evalua en cualquier punto.
Usa cut() para asignar bins a w nuevos. Resuelve la limitacion de thalmle.

**Q-17 (RESUELTA) -- Hueco en el borde izquierdo**
make_knots() fuerza min(A) como primer nudo.
El primer spline (A - min(A))_+ = A - min(A), que combinado con el intercepto
genera el span de la funcion lineal pura A desde el primer dato.
No hay hueco plano en la cola izquierda.

## Preguntas abiertas

**Q-06 -- Cross-fitting**
Evidencia de thalmle: cross-fitting K=2 empeora el RMSE con tratamiento normal
(de ~0.4 a ~3.6). El repo nativo usa ajuste en muestra completa por defecto.
Se puede implementar como opcion futura si hay evidencia de que ayuda en otros DGP.

**Q-07 -- Termino Q_W del gradiente**
Pendiente de implementar la componente Q_W en influence_curve().
Actualmente solo se calcula la componente Q_bar.
La version actual usa la medida empirica de A (D-0005).

**Q-08 -- CV honesto de lambda**
El CV del targeting usa Q_n ajustado en la muestra completa (no honesto).
Evidencia de thalmle: el CV honesto no mejora y complica la implementacion.

**Q-10 -- Forma del LASSO**
Usamos la forma penalizada de Lagrange (no la restriccion ||eps||_1 <= lambda).
Ambas dan el mismo camino de soluciones. El selector de lambda opera sobre
el camino completo, no asume monotonicidad de n(lambda).

**Q-11 -- Colas de omega y g_n ~ 0**
omega: KDE con floor = 1e-4. Sesgo de borde documentado en tests.
g_n: KDE por bins con floor = 1e-4.
w_cap recorta los pesos omega/g para evitar inestabilidad.

**Q-15 -- Reescalado a [0,1]**
Pendiente. Actualmente se usan valores crudos de A y W.
Los terminos de borde (a, w, a*w) dependen del origen.
No bloquea v0 pero debe resolverse antes de comparar con el paper.

**Q-16 -- Base de hal9001 vs nativa**
No se ha comparado directamente (D-0010: sin dependencias externas).
Nuestra implementacion incluye explicitamente los terminos de borde
que hal9001 podria omitir. La fidelidad se verifica contra el paper.
