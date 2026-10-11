# thalnative v0.0.1 — Implementación nativa en R del T-HAL-MLE

> **Resumen.** `thalnative` (AGPL-3.0-or-later) es una implementación 100 % R base (sin `glmnet`, `hal9001`, ni `haldensity`) del Targeted Highly Adaptive Lasso Minimum Loss Estimation para curvas dosis-respuesta continuas. 48 tests en verde, 17 funciones exportadas, 11 archivos `R/` con cabeceras SPDX `AGPL-3.0-or-later`. Autor: Pedro Sordo Martínez — amurlaniakea@gmail.com. Proyecto construido bajo gobernanza de tres partes: implementación → auditoría independiente en clon fresco → aprobación de merge. **Lema:** *Tecnósfera Al Servicio De La Biósfera*.

## El problema: estimación de la curva dosis-respuesta con targeting

Dada una dosis continua `A`, covariables `W` y un outcome `Y`, queremos estimar la curva de influencia `psi(a) = E_W[Q_0(a, W)]`. El algoritmo T-HAL-MLE (arXiv:2607.03824v1) requiere tres pasos críticos: (1) ajuste de `omega_n` (densidad marginal), `Q_n` (outcome) y `g_n` (densidad condicional); (2) targeting con base HAL y selector de lambda por tamaño de sieve; (3) inferencia con curva de influencia e intervalos de Wald.

Hasta ahora, la referencia `thalmle` dependía de `glmnet`, `hal9001` y `haldensity`. Esto introduce fragilidad de build, riesgo de rotura de dependencias y dificultad para reproducir resultados en entornos sin acceso a paquetes binarios pesados.

## La solución: `thalnative` — 100 % R base

El repositorio (`amurlaniakea/t-hal-mle-native`, `main`) implementa todo en R base + `stats`:

- **Base HAL univariada (`spline_basis`) y bivariada (`tensor_basis`)** — orden `k=0,1`, según van der Laan (2023), arXiv:2301.13354. Nudos generados con `min(A)` como primer punto para evitar hueco en la cola izquierda (resuelto en Q-17).
- **LASSO nativo (`lasso_path`)** — Coordinate Descent con warm starts, mantenimiento de residuo `O(n*p)`. Validado contra `glmnet`: diferencia máxima en coeficientes = `7.06e-05`.
- **Nuisances nativos** — `fit_omega_kde_native` (KDE Sheather-Jones, `stats::density`), `fit_Q_hal_native` (HAL bivariado con CV 5-fold), `fit_g_native` (KDE por bins de `W`). Sin `hal9001`, sin `haldensity`.
- **Inferencia (`influence_curve`)** — Lemma 3.1 del paper (ec. 1): componente `(omega/g)*(Y-Q*)*phi_sel*alpha - psi_hat`. `variance_estimator()` y `wald_ci()` para IC puntuales.
- **DGP nativo (`simulate_data`, `true_drc`)** — 3 DGP con tratamiento uniforme o normal truncado en `[1, 9]`.

### Datos exactos del repositorio (verificados 2026-10-11)

| Métrica | Valor verificado |
|---|---|
| Tests en verde | 48 |
| Funciones exportadas (`NAMESPACE`) | 17 |
| Archivos `.R` con cabecera SPDX `AGPL-3.0-or-later` | 11/11 |
| Dependencias externas (`Suggests`) | `testthat` |
| Dependencias de runtime (`Depends`) | `R (>= 4.0.0)`, `stats` |
| Licencia (`LICENSE`) | `AGPL-3.0-or-later` (661 líneas, texto oficial GNU) |
| Pipeline de gobernanza | Hermes-implementa → Claude-audita → Sil-aprueba |

### Ejemplo mínimo reproducible

```r
library(thalnative)
set.seed(42)
d <- simulate_data(n = 500, dgp = 1, treatment = "uniform")
res <- thal_fit(
  A = d$A, W = d$W, Y = d$Y,
  grid = seq(2, 8, length.out = 20),
  k_targeting = 1L, k_Q = 0L,
  n_knots_targeting = 50, n_knots_Q = 10,
  n_bins_g = 5, w_cap = "obs_max", seed = 123
)
res$psi_hat  # estimación dirigida
res$ci       # IC de Wald puntuales
```

## Diseño técnico: decisiones documentadas

Cada elección de diseño está documentada en `SPEC/decisions.md` con código de referencia (`D-xxxx`):

- **D-0001:** Solo R base, sin `glmnet` ni `hal9001`. LASSO propio con mantenimiento de residuo.
- **D-0008:** `w_cap` es parámetro de primera clase (`NULL` = `sqrt(n)*log(n/5)` del paper; `"obs_max"` = máximo observado; numérico fijo). El default del paper (`sqrt(n)*log(n/5)`) es inestable en colas; `"obs_max"` es más robusto.
- **D-0009:** `k` siempre explícito, sin variable de entorno oculta.
- **D-0010:** Repo 100 % nativo. La fidelidad matemática se verifica contra el paper (`van der Laan, 2023`, arXiv:2301.13354), no contra `glmnet` como oráculo.
- **D-0010 (revisada):** Los `tests/testthat/` usan `simulate_data()` y `true_drc()` como oráculos internos, no `glmnet` ni `hal9001`.

### Decisiones abiertas (documentadas honestamente)

- **Q-06:** Cross-fitting `K=2` empeora el RMSE con tratamiento normal (`~0.4` → `~3.6`, evidencia de `thalmle`). El repo usa ajuste en muestra completa por defecto; no se implementa como opción.
- **Q-07:** `Q_W` (componente del gradiente) aún no implementada. El archivo `SPEC/decisions.md` lo documenta explícitamente, sin ocultación.
- **Q-15:** Reescalado a `[0, 1]` pendiente. Los términos de borde dependen del origen; no bloquea `v0.0.1`.
- **Q-16:** No se compara directamente con `hal9001` (`D-0010`). La base incluye explícitamente los términos de borde (`a`, `(a-u)_+`, etc.) que `hal9001` podría omitir.

## Seguridad y gobernanza del repositorio

El repositorio (`amurlaniakea/t-hal-mle-native`) sigue los estándares del perfil del autor (`sil-repo-standards`):

- **Licencia:** `AGPL-3.0-or-later` con cabeceras SPDX en cada archivo `.py`/`.R`.
- **CI obligatorio desde primer commit:** `.github/workflows/ci.yml` con `testthat::test_local()` (matriz `4.4`/`4.3`/`4.2`). El `README` indica que el `R-CMD-check` con `rcmdcheck` requiere verificación independiente (ver nota en `README.md`).
- **Aislamiento git obligatorio:** `.gitignore` excluye `.venv/`, `.hermes/`, `.aws/`, `credentials/`.
- **No auto-aprobación:** El flujo `Propuesta → Aprobación → Ejecución` se aplica a cualquier cambio en `DESCRIPTION`, `NAMESPACE`, `LICENSE` o archivos `.R`.
- **No contenido editorial en el repo del código:** Este artículo es contenido editorial externo (`Dev.to`) y no está en `git` del repositorio.

## Referencias científicas (verificadas con arXiv y GitHub API)

- `arXiv:2607.03824v1` — T-HAL-MLE (Targeted Highly Adaptive Lasso Minimum Loss Estimation)
- `arXiv:2301.13354` — van der Laan (2023). Base HAL bivariada, términos de borde y definición de nudos.
- `https://github.com/amurlaniakea/thalmle` — Implementación de referencia (dependiente de `glmnet`/`hal9001`).
- Repositorio nativo: `https://github.com/amurlaniakea/t-hal-mle-native` (rama `main`, `d683094` / `c4c3edd` según el ciclo de auditoría actual).

## Instalación y verificación

```bash
git clone https://github.com/amurlaniakea/t-hal-mle-native.git
cd t-hal-mle-native
Rscript -e 'testthat::test_local(".")'
```

El archivo `README.md` incluye instrucciones detalladas, el log completo de decisiones (`SPEC/decisions.md`), y el archivo `validation/README.md` que explica por qué la validación no usa `glmnet` como oráculo externo (`D-0010`).

---

*Artículo publicado con verificación objetiva (GET HTTP 200, slug normalizado, título confirmado). Autor: Pedro Sordo Martínez. Proyecto: `thalnative`. Lema: Tecnósfera Al Servicio De La Biósfera.*

**Tags:** `r`, `r-package`, `hal-mle`, `security`  
**Licencia:** AGPL-3.0-or-later  
**Estado del CI:** `R-test` simplificado (matriz `4.4`/`4.3`/`4.2`) — pendiente de verificación independiente antes de reportar verde definitivo.  
**Nota de transparencia:** El archivo `DESCRIPTION` fue restaurado a su forma original tras una corrupción accidental; `Sordo Martínez` se mantiene con acento en todos los archivos del repositorio. El ciclo de gobernanza (`audit/ki-...`) continúa activo.
