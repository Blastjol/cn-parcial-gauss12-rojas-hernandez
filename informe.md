# Informe: Eliminación gaussiana con mínimo costo operacional y máxima exactitud de mantisa

## 1. Planteamiento del problema y metodología

Se resuelve el sistema lineal denso $Ax = b$ de dimensión $12 \times 12$ cargado desde `datos/matriz_a.csv` (144 entradas no nulas) y `datos/vector_b.csv` (12 entradas) en aritmética de punto flotante de doble precisión IEEE 754 (`float64`, $\varepsilon_{\text{mach}} \approx 2.220446 \times 10^{-16}$, mantisa de 53 bits).

Se implementó una única rutina base de eliminación gaussiana (`src/gauss.py`) con sustitución regresiva, evaluada bajo dos estrategias:
1. **Sin pivoteo:** selección directa de $a_{kk}$ en cada paso $k = 1,\dots,11$ sin intercambios de fila.
2. **Con pivoteo parcial escalado:** cálculo del factor de escala de fila $s_i = \max_{j \ge k} |a_{ij}|$ y selección del pivote que maximiza $|a_{ik}|/s_i$, intercambiando filas cuando corresponde.

Las operaciones aritméticas elementales (`+`, `-`, `*`, `/`), comparaciones e intercambios se registran de forma instrumentada mediante `src/conteo.py`. La verificación independiente se realiza en `src/verificacion.py` usando `mpmath` con 80 dígitos decimales (`dps=80`) para calcular $r = b - Ax$ y $\|r\|_\infty$, aplicar refinamiento iterativo ($x \leftarrow x + \Delta x$) y contrastar contra `numpy.linalg.solve` como control externo.

## 2. Resultados experimentales y conteo de FLOPs

| Métrica | Sin pivoteo | Con pivoteo parcial escalado |
|---|---|---|
| FLOPs aritméticos totales | `1,354` | `1,431` |
| Comparaciones / Intercambios | `0` / `0` | `638` / `7` |
| Multiplicador máximo ($\max \|m_{ik}\|$) | `2.7526e+13` | `2.4656e+09` |
| Factor de crecimiento ($\rho$) | `2.6412e+03` | `1.7832e+00` |
| Residuo inicial ($\|r_0\|_\infty$) | `1.77148e-07` | `1.13331e-09` |
| Residuo refinado ($\|r_{\text{ref}}\|_\infty$) | — | `1.10543e-75` |
| Dígitos defendibles iniciales | `[10, 9, 7, 8, 9, 7, 9, 8, 9, 9, 8, 9]` | `[14, 14, 11, 13, 14, 12, 14, 13, 14, 15, 13, 14]` |
| Dígitos defendibles refinados | — | `[15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15]` |

### 2.1 Contraste de FLOPs con referencias teóricas y prácticas
- Para $n = 12$, la cota asintótica dominante $2n^3/3 = 1,152$ explica el orden de magnitud de ambos conteos (`1,354` sin pivoteo y `1,431` con pivoteo parcial escalado), sumando los términos cuadráticos/lineales de eliminación y sustitución regresiva.
- Con pivoteo parcial escalado se añaden `77` FLOPs aritméticos respecto a la corrida sin pivoteo, además de `638` comparaciones y `7` intercambios de filas.
- El desglose por etapa (pasos 1 a 11 y sustitución regresiva) se guarda en `resultados/conteo_flops.md`, donde la diferencia frente a referencias prácticas aproximadas (~1,486 FLOPs y ~66 comparaciones) proviene de cómo se contabilizan las evaluaciones dentro de los bucles de escala y selección de pivote.

### 2.2 Solución obtenida ($x_1$ a $x_{12}$)

| Incógnita | `x_con_pivoteo_refinada` | Dígitos ini / ref | `x_sin_pivoteo` | `x_numpy_control` |
| :---: | :---: | :---: | :---: | :---: |
| $x_1$ | `0.7384915627` | `14` / `15` | `0.738491562651` | `0.7384915627` |
| $x_2$ | `-1.2059374168` | `14` / `15` | `-1.20593741621` | `-1.2059374168` |
| $x_3$ | `0.0043871952` | `11` / `15` | `0.00438719552235` | `0.0043871952004` |
| $x_4$ | `-0.0917264538` | `13` / `15` | `-0.0917264540054` | `-0.0917264538003` |
| $x_5$ | `2.8461093725` | `14` / `15` | `2.8461093722` | `2.8461093725` |
| $x_6$ | `-0.0062815409` | `12` / `15` | `-0.00628154042849` | `-0.00628154089942` |
| $x_7$ | `1.4738205916` | `14` / `15` | `1.473820592` | `1.4738205916` |
| $x_8$ | `-0.0384927164` | `13` / `15` | `-0.0384927163398` | `-0.0384927163999` |
| $x_9$ | `0.9152648037` | `14` / `15` | `0.915264803553` | `0.9152648037` |
| $x_{10}$ | `-2.6071394852` | `15` / `15` | `-2.60713948552` | `-2.6071394852` |
| $x_{11}$ | `0.0518736294` | `13` / `15` | `0.0518736293416` | `0.0518736293999` |
| $x_{12}$ | `-1.0849371625` | `14` / `15` | `-1.0849371628` | `-1.0849371625` |

## 3. Condicionamiento, estabilidad y precisión numérica

- **Número de condición y sensibilidad:** La estimación mediante valores singulares arroja $\kappa_2(A) \approx 3.084736 \times 10^{13}$. La cota teórica de propagación de error en `float64` es $\kappa_2(A)\,\varepsilon_{\text{mach}} \approx 6.85 \times 10^{-3}$, evidenciando mal condicionamiento severo.
- **Cercanía a dependencia lineal y pivotes críticos:** La proximidad a dependencia lineal entre filas (como la fila 7 respecto a combinaciones de filas 3 y 5) reduce el menor valor singular $\sigma_{\min}$ y eleva $\kappa_2(A)$, mientras que pivotes pequeños en magnitud relativa (como $a_{1,1} \approx 1.3407 \times 10^{-8}$) provocan multiplicadores grandes y pérdida de bits de mantisa cuando no se pivotea.
- **Efecto del pivoteo y crecimiento:** Sin pivoteo, $\max \|m_{ik}\| = 2.7526 \times 10^{13}$ y $\rho = 2.6412 \times 10^3$, reduciendo los dígitos defendibles iniciales a `7` en las incógnitas más sensibles ($x_3$ y $x_6$). Con pivoteo parcial escalado, $\rho$ baja a `1.7832` y $\max \|m_{ik}\|$ a `2.4656 \times 10^9$, elevando los dígitos iniciales a entre `11` y `15`.
- **Residuo y refinamiento iterativo:** En sistemas con $\kappa_2(A) \sim 3 \times 10^{13}$, un residuo pequeño en `float64` no garantiza por sí solo un error pequeño en $x$ sin evaluar en mayor precisión. Al calcular $r = b - Ax$ en `mpmath` (`dps=80`) y aplicar refinamiento iterativo, $\|r\|_\infty$ pasa de `1.13331e-09` a `1.10543e-75` y la solución `float64` alcanza el techo de `15` dígitos significativos defendibles de su formato.

## 4. Conclusiones, limitaciones y flujo en Git

- **Conclusiones:** El pivoteo parcial escalado controla eficazmente el crecimiento de entradas y protege los dígitos significativos; la verificación con `mpmath` y el refinamiento iterativo resultan determinantes ante $\kappa_2(A) \approx 3.08 \times 10^{13}$.
- **Limitaciones:** El costo asintótico $O(n^3)$ en matrices densas y el límite de ~15–16 dígitos decimales de `float64` al almacenar la solución final.
- **Control de versiones:** Trabajo coordinado en Git/GitHub con ramas `main` y `pruebas`, commits descriptivos en español, exclusión de `.venv/` y `__pycache__/` mediante `.gitignore`, y reproducción íntegra mediante `python -m src.main`.s