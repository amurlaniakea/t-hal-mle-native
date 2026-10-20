# Decision log y preguntas abiertas
Copyright (c) Pedro Sordo Martínez. Licencia: AGPL-3.0-or-later.
Formato: Q-xx = pregunta (RESUELTA o abierta). D-xxxx = decision tomada.
Ultima actualizacion: 2026-10-11 (v0.0.1, 53 tests, 20 funciones exportadas).

## Decisiones

**D-0001 -- Dependencias de runtime**
Estado: tomada (confirmada).
Solo R base y stats. Sin Matrix, glmnet, hal9001, ni haldensify.
LASSO propio via Coordinate Descent con warm starts y mantenimiento de residuo.

**D-0002 -- Oraculos de validacion (REVISADA por D-0010)**
Estado: eliminada.
No se usan dependencias externas ni siquiera como oraculos de validacion.
La fidelidad matematica se verifica contra van der Laan (2023) mediante tests unitarios.

**D-0003 -- Alcance v0**
Estado: tomada.
Y continuo (perdida cuadratica), W y A univariados, k_targeting = 1.
Y binario y W multivariado quedan para versiones posteriores.

**D-0004 -- Inferencia**
Estado: tomada.
Inferencia primaria = curva de influencia (Lemma 3.1), IC de Wald puntuales.
Implementada con componentes Q_bar (en A_i observadas, con el MISMO Q* del estimador puntual)
y Q_W, ambas con Sigma^+ sobre la base seleccionada.
Estado: EXPERIMENTAL. Con nuisances oraculo (Q mal especificado, DGP 1 uniforme, n = 300,
80 replicas, 5 dosis) la media del IC es ~0, pero el SE estimado subestima la SD de Monte Carlo
(p. ej. 0.22 vs 0.66 en una dosis central) y la cobertura de Wald media fue 0.73 (nominal 0.95).
Causas candidatas, no probadas: seleccion LASSO data-adaptativa tratada como sieve fijo,
residuos dentro de muestra (el paper hace cross-fitting K = 2 para la varianza).
Actualizacion 2026-10-14: el 0.73 se midio con el DGP anterior (ver D-0011). Con los DGP exactos y nuisances
estimadas (piloto n = 500, B = 50) la cobertura de Wald va de 0.72 a 0.92 con el tope del paper y de 0.73 a 0.87 con
topes 2 y 3.
Actualizacion 2026-10-15: con cross-fitting K = 2 (n_cf = 2; residuos del Q inicial fuera de fold) la cobertura de
Wald sube 0.02-0.10 de forma pareada (n = 500, tope 3, normal) y el SE medio / sd Monte Carlo pasa de 0.73-0.74 a
0.88-1.00 en los DGP 1 y 2; sigue por debajo de 0.95 (0.85-0.90). Con g y omega oraculo casi no cambia. Ver README y
experiments/results/cf_pareado_n500_B60_s222001.csv.
Actualizacion 2026-10-15 (2): con tratamiento uniforme el cross-fitting tambien sube la cobertura (+0.01 a +0.07).
Con la ventana por defecto del sieve (tamano ~21) la cobertura con un SE ideal queda en 0.90-0.93 por sesgo; con
c1 = 9, c2 = 14 (tamano ~32) llega a ~0.95, pero la cobertura de Wald sigue en 0.86-0.89 en los DGP 1 y 2.
Actualizacion 2026-10-17: gran parte del defecto del SE en los DGP 1 y 2 venia de la resolucion del HAL de Q
(k_Q = 0, n_knots_Q = 10): picos de varianza real cerca de los nudos que el IC no recoge. Con n_knots_Q = 20 o
k_Q = 1 desaparecen y la cobertura de Wald con el SE en muestra llega a 0.93-0.97 (uniforme). Ver D-0012.

**D-0005 -- Medida de integracion de A**
Estado: tomada.
Usar la medida empirica de A en la proyeccion final.

**D-0006 -- Selector de lambda**
Estado: exponente por k implementado; seleccion por CV (lambda_cv recortado a la ventana) NO implementada:
se toma el primer lambda de la ventana.
El exponente depende del orden de suavidad k:
- k=0: exponente 1/3 (funciones escalon, menos regularidad).
- k=1: exponente 1/5 (splines lineales, mas regularidad).
Implementado en select_lambda_by_size().

**D-0007 -- Afirmaciones tecnicas**
Estado: tomada.
Prohibido afirmar equivalencia con glmnet o validez de IC sin test citado.
Formato: "En el test X, con DGP Y, n = Z, el RMSE fue E".

**D-0008 -- Truncamiento de pesos**
Estado: tomada (corregida).
w_cap es un parametro de primera clase.
Default: NULL (usa sqrt(n)*log(n/5) del paper).
Opcion "obs_max" disponible explicitamente para casos con colas pesadas.
Evidencia Monte Carlo (n = 500, tratamiento normal): con el tope del paper (~103) T-HAL tiene RMSE 4-7 frente a 1.2-1.3
del plug-in HAL en los DGP 1 y 2; con tope 3 gana (-0.34 a -0.39 de RMSE por dosis) y en el DGP 3 no supera al HAL
(+0.25 a +0.35; ~0 con tope 1). Bajar el tope mejora el RMSE pero reduce la cobertura de Wald. El default sigue
siendo el del paper; ver README.

**D-0009 -- Orden de suavidad**
Estado: tomada.
El orden k es siempre un argumento explicito, sin variable de entorno oculta.

**D-0010 -- Repo 100% nativo**
Estado: tomada.
No se usan dependencias externas. Los nuisances se estiman con funciones nativas:
KDE para omega, HAL bivariado propio para Q, KDE por bins para g.

## Preguntas resueltas

**Q-01 (RESUELTA) -- Base HAL de orden k**
Fuente: van der Laan (2023), arXiv:2301.13354.
Los terminos de borde son potencias (nudo en el origen), NO indicadores.
CORRECCION: Se implementa escalado interno de A y W a [0,1] para que el "borde 0"
tenga sentido geometrico y evitar colinealidad con el intercepto o nudos extremos.

**Q-02 (RESUELTA) -- Nudos de la base del targeting**
make_knots() genera nudos estrictamente interiores (excluyendo 0 y 1 en la escala [0,1])
para evitar colinealidad con los terminos de borde.

**Q-03 (RESUELTA) -- HAL del outcome**
HAL bivariado nativo con base tensor_basis() (escalada a [0,1]^2) y CV 5-fold (lambda.min).

**Q-04 (ABIERTA) -- Densidad condicional g_n**
OJO: el paper usa un estimador HAL de la densidad (haldensify). Esto es una alternativa no validada.
KDE Sheather-Jones por bins de W (n_bins = 5 por defecto).
Cada bin tiene su propia KDE 1D de A. Fallback a la marginal si < 5 datos en un bin.
Es una funcion general de (a, w) evaluable fuera de muestra (resuelve Q-14).

**Q-05 (RESUELTA) -- Sigma inversa**
Sin Sigma^{-1} por defecto. La covariable inteligente es omega/g * phi(A).

**Q-06 (RESUELTA) -- Cross-fitting**
El dato de thalmle (CF K=2 con tratamiento normal: RMSE ~3.69 frente a ~0.99 sin CF) corresponde a otra variante.
Aqui `n_cf = 2` es cross-fitting SOLO de la varianza: no cambia psi_hat y usa residuos del Q inicial fuera de fold.
Implementado y medido (ver D-0004 y README): su efecto sobre la cobertura depende de la resolucion del HAL de Q.

**Q-07 (RESUELTA) -- Termino Q_W del gradiente**
Implementado en influence_curve(): promedio sobre los puntos de proyeccion (la medida de A del paso 7)
de phi(a_j) * (Q*(a_j, W_i) - psi*_pl(a_j)), con psi*_pl la marginal plug-in dirigida.
Su media empirica es exactamente 0 (test en test-influence.R).

**Q-09 (RESUELTA) -- Truncamiento de pesos**
Se trunca omega/g. Se aplica siempre. Ver D-0008.

**Q-10 (RESUELTA) -- Forma del LASSO**
Usamos la forma penalizada de Lagrange. El selector de lambda opera sobre el camino
completo, no asume monotonicidad de n(lambda).

**Q-11 (RESUELTA) -- Colas de omega y g_n ~ 0**
omega: KDE con floor = 1e-4.
g_n: KDE por bins con floor = 1e-4.
w_cap recorta los pesos omega/g para evitar inestabilidad.

**Q-12 (RESUELTA) -- Criterios de convergencia del LASSO**
Coordinate Descent con warm starts, tolerancia 1e-6, max 200 iteraciones.
Mantenimiento del residuo para coste O(n*p) por barrido.

**Q-13 (RESUELTA) -- smoothness_orders en DGP 1 y 2**
El orden k es un argumento explicito (D-0009). k=0 para discontinuidades, k=1 para suaves.

**Q-14 (RESUELTA) -- g_n evaluable fuera de muestra**
fit_g_native() devuelve una funcion general de (a, w) que evalua en cualquier punto
usando cut() para asignar bins a w nuevos.

**Q-15 (RESUELTA) -- Reescalado a [0,1]**
RESUELTO: spline_basis() y tensor_basis() ahora escalan internamente las variables
a [0,1] usando el rango de los datos de entrenamiento, garantizando consistencia
entre el ajuste y la evaluacion en la grid.

**Q-16 (RESUELTA) -- Base de hal9001 vs nativa**
Nuestra implementacion incluye explicitamente los terminos de borde correctos
y evita la colinealidad mediante nudos interiores y escalado.

**Q-17 (RESUELTA) -- Hueco en el borde izquierdo**
RESUELTO: make_knots() ya NO fuerza min(A) como nudo, ya que eso generaba colinealidad
con el termino lineal (x - min(x))_+ = x - min(x). Ahora los nudos son estrictamente
interiores, y el termino lineal puro (escalado) cubre el borde 0 correctamente.

## Preguntas abiertas

- **Q-08 -- CV honesto de lambda**: El CV del targeting usa Q_n ajustado en la muestra completa.
- **Mejora de g_n**: KDE por bins es una aproximacion. Futuras versiones podrian explorar
  estimadores de densidad condicional mas sofisticados en R base.

---

## D-0011 -- DGP del paquete = DGP exactos del paper
Estado: tomada.
`simulate_data()`, `Q0_dgp()`, `g0_dgp()`, `omega_oracle()` y `true_drc()` reproducen la seccion 7.1 del
paper: W ~ N(0,1); g(w) = 5*1[w>0] (DGP 1), 4(1[w>-1]-1[w>0]+1[w>1]) (DGP 2), 3 sin(7 pi w) (DGP 3);
A ~ U(0,10) o A|W ~ N(5+1.5W, 2^2) truncada a [0,10]; `omega_oracle` es la marginal exacta (integral numerica).
Antes el paquete usaba g(w) propias (0.5w + 0.25w^2); la cobertura ~0.73 citada en D-0004 se midio con ese DGP anterior.
Con el DGP exacto y tratamiento normal, T-HAL se dispara en las dosis extremas incluso con g y omega oraculo
(ver experiments/mc_wcap.R); el test de test-core.R solo exige un resultado finito y acotado.

---

## D-0012 -- Resolucion del HAL de Q (n_knots_Q, k_Q)
Estado: abierta (valor por defecto sin cambiar: n_knots_Q = 10).
Con k_Q = 0 y 10 nudos por variable (cuantiles de la muestra) la varianza real de psi_hat tiene picos en las dosis
cercanas a los nudos, proporcionales a la pendiente al cuadrado de la verdad, que el IC no recoge. Se eliminan con
n_knots_Q = 20 (coste ~3 veces mayor) o con k_Q = 1 (sin coste, pero el SE sale ~1.5 veces sobrestimado). Pendiente:
confirmar con otros n y mas replicas, decidir el default (p. ej. escalar n_knots_Q con n o acelerar el LASSO) y
comprobar el HAL del paper. Los DGP 3 con tratamiento normal siguen sin reproducir la ventaja de T-HAL (ver README).

---

## D-0013 -- Solapamiento (A | W) y DGP 3 con tratamiento normal
Estado: abierta.
Con sd(A | W) = 2 T-HAL pierde frente al HAL en el DGP 3 normal (+0.25 a +0.40 de RMSE por replica) y casi todo el
exceso viene de las dosis altas (a >= 7.67), donde P(A > 9 | W < -0.5) ~ 0.003. Ni el tope de pesos (2, 1.5, 1) ni el
tamano del modelo de targeting ni los oraculos de g y omega lo arreglan; ensanchar A | W (sd = 3, 4) lo invierte
(+0.02, -0.11). Pendiente: comprobar con el paper la sd o varianza de A | W (D-0011), el truncamiento, la rejilla de
evaluacion y la agregacion del RMSE; decidir si la evaluacion se restringe al soporte (declarado de antemano) y como
se reporta la cobertura en las dosis sin solapamiento.
Actualizacion 2026-10-18: el paper confirma sd = 2, el truncamiento y la rejilla. Con su metrica (mediana del RMSE por dosis
sobre las 25 dosis) T-HAL queda por debajo del HAL en el DGP 3 normal; la perdida aparece con metricas que pesan los
extremos (media por dosis, RMS por replica). Ver D-0014.

---

## D-0014 -- Metrica de comparacion con el paper
Estado: tomada.
Las comparaciones con el paper se hacen con sus metricas (seccion 7.3): RMSE y sesgo PUNTUALES por dosis, resumidos con la
mediana y los cuartiles sobre las 25 dosis, proporcion de dosis con |sesgo / MC-SE| <= 1/log(n), cobertura con el MC-SE y
cobertura de Wald (experiments/resumen_tipo_paper.R). La media de los RMSE por dosis y la diferencia de RMSE por replica
(mc_summary.R) se mantienen como analisis de sensibilidad porque pesan las dosis extremas, donde T-HAL tiene mucha varianza.

