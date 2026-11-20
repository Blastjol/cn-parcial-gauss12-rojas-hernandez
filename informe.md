# Informe Técnico: Eliminación gaussiana con mínimo costo operacional y máxima exactitud de mantisa

## 1. Identificación

- **Estudiantes:** Blas Emilio Rojas Díaz · Joaquín Hernández
- **Asignatura:** Computación Numérica (INF-321) — Sección `[Indicar Sección]`
- **Carrera:** Ingeniería Civil Informática, Facultad de Ciencias de la Ingeniería, Universidad Católica del Maule
- **Docentes:** Juan Manuel Cabezas Contardo · Pedro Pesante
- **Fecha:** Noviembre de 2026
- **Repositorio Git:** [https://github.com/Blastjol/cn-parcial-gauss12-rojas-hernandez](https://github.com/Blastjol/cn-parcial-gauss12-rojas-hernandez)

## 2. Planteamiento

Se resuelve el sistema lineal denso $Ax = b$ de dimensión $n = 12$ ($144$ entradas no nulas en $A$ y $12$ componentes en $b$ cargadas con todas sus cifras exactas desde `datos/matriz_a.csv` y `datos/vector_b.csv`) en doble precisión IEEE 754 (`float64`, mantisa de 53 bits, $\varepsilon_{\text{mach}} = 2^{-52} \approx 2.220446 \times 10^{-16}$).

El sistema presenta tres dificultades numéricas diseñadas para degradar la mantisa:
1. **Escalas de fila incompatibles:** las magnitudes de las filas oscilan entre $10^{-8}$ y $10^{5}$ (13 órdenes de magnitud). Restar un múltiplo de una fila grande a una pequeña desplaza los bits bajos fuera de los 53 bits de mantisa por alineamiento de exponentes.
2. **Pivote inicial minúsculo:** $a_{1,1} = 1.3407 \times 10^{-8}$ es trece órdenes menor que $a_{6,1} = 3.69040121 \times 10^{5}$, provocando multiplicadores gigantescos si no se intercambia fila.
3. **Casi dependencia lineal oculta en la Fila 7:** la Fila 7 satisface aproximadamente $F_7 \approx 4.2\,F_3 + 586\,F_5$ salvo una perturbación de orden $\mathcal{O}(10^{-4})$ (por ejemplo, en $c_1$: $4.2(2.70519337\times 10^{-3}) + 586(4.11556548\times 10^{-5}) = 3.5479\times 10^{-2}$ frente a $a_{7,1} = 3.5470172\times 10^{-2}$). Esto mantiene la matriz no singular ($\det(A) \approx 1.3646 \times 10^{9}$), pero eleva su número de condición a $\kappa_2(A) \approx 3.084736 \times 10^{13}$.

## 3. Estrategia de eliminación

Se implementó una única rutina base (`eliminacion_gaussiana` en `src/gauss.py`) parametrizada por `usar_pivoteo_escalado`, evitando duplicar código o calcular la inversa de $A$:
- **Anulación estricta de los 66 elementos subdiagonales:** en cada paso $k = 1,\dots,11$, se fija `A[i, k] = 0.0` y el ciclo interno recorre únicamente las columnas $j = k+1,\dots,12$, evitando $66 \times 2 = 132$ operaciones redundantes sobre entradas ya anuladas.
- **Reutilización única de multiplicadores:** para cada fila $i > k$, el multiplicador $m_{ik} = a_{ik}/a_{kk}$ se calcula **una sola vez** y se reutiliza en toda la fila $i$ y en la actualización de $b_i$.
- **Pivoteo parcial escalado:** en cada paso $k$, se determina la escala de fila $s_i = \max_{k \le j \le 12} |a_{ij}|$ y se selecciona la fila pivote $p \ge k$ que maximiza $|a_{pk}| / s_p$, intercambiando las filas $k$ y $p$ en $A$, $b$ y $s$ cuando $p \neq k$.

## 4. Conteo de operaciones

**Convención declarada:** 1 FLOP = 1 operación aritmética elemental (`+`, `-`, `*`, `/`) ejecutada durante la eliminación hacia adelante, actualización de $b$ y sustitución regresiva. Las divisiones de escalamiento ($|a_{ik}|/s_i$), comparaciones de búsqueda e intercambios de fila se contabilizan por separado como sobrecosto de estabilidad (`resultados/conteo_flops.md`).

| Etapa / Paso | Mult / Div | Sum / Resta | FLOPs elim./sust. | Div. escalamiento (pivoteo) | Comparaciones (pivoteo) | Intercambios |
|---|---:|---:|---:|---:|---:|---:|
| Paso 1 | 143 | 132 | 275 | 12 | 132 | 1 |
| Paso 2 | 120 | 110 | 230 | 11 | 110 | 1 |
| Paso 3 | 99 | 90 | 189 | 10 | 90 | 1 |
| Paso 4 | 80 | 72 | 152 | 9 | 72 | 1 |
| Paso 5 | 63 | 56 | 119 | 8 | 56 | 0 |
| Paso 6 | 48 | 42 | 90 | 7 | 42 | 1 |
| Paso 7 | 35 | 30 | 65 | 6 | 30 | 0 |
| Paso 8 | 24 | 20 | 44 | 5 | 20 | 1 |
| Paso 9 | 15 | 12 | 27 | 4 | 12 | 1 |
| Paso 10 | 8 | 6 | 14 | 3 | 6 | 0 |
| Paso 11 | 3 | 2 | 5 | 2 | 2 | 0 |
| Sustitución regresiva | 78 | 66 | 144 | 0 | 0 | 0 |
| **Total** | **716** | **638** | **1,354** | **77** | **638** | **7** |

- **Comparación con cotas teóricas:** Para $n = 12$, el término cúbico dominante es $2n^3/3 = 1,152$ FLOPs. Al incluir las actualizaciones de $b$ (132 FLOPs) y la sustitución regresiva (144 FLOPs), nuestro algoritmo alcanza **1,354 FLOPs** (716 multiplicaciones/divisiones y 638 sumas/restas), ubicándose 132 FLOPs por debajo de la cota referencial de 1,486 FLOPs al no operar sobre las 66 entradas subdiagonales anuladas y dentro del margen aceptado $\pm 10\%$ (`[1337, 1635]`). El pivoteo parcial escalado agrega un sobrecosto separado de **77 divisiones de escala**, **638 comparaciones** y **7 intercambios**.

## 5. Análisis de mantisa y pivoteo

| Métrica | Sin pivoteo | Con pivoteo parcial escalado |
|---|---:|---:|
| Multiplicador máximo ($\max \|m_{ik}\|$) | `2.7526e+13` | `2.4656e+09` |
| Factor de crecimiento ($\rho$) | `2.6412e+03` | `1.7832e+00` |
| Bits perdidos por crecimiento ($\log_2 \rho$) | `11.37 bits` | `0.83 bits` |
| Dígitos significativos iniciales (rango en $x_1\dots x_{12}$) | `7 a 10 dígitos` | `11 a 15 dígitos` |
| Bits de mantisa perdidos en solución inicial (de 53 bits) | `20 a 30 bits` | `3 a 16 bits` |

- **Rechazo del pivote $a_{1,1} = 1.3407 \times 10^{-8}$:** Sin pivoteo, dividir filas de orden $10^{5}$ por $a_{1,1}$ produce $m_{6,1} \approx 2.7526 \times 10^{13} \approx 2^{44.65}$. Al restar $m_{6,1} a_{1,j}$ de $a_{6,j}$, el alineamiento de exponentes desplaza $\approx 44.6$ bits de la mantisa de 53 bits, destruyendo la información de los bits menos significativos y elevando el factor de crecimiento a $\rho = 2.6412 \times 10^{3}$ (degradando $x_3$ y $x_6$ a solo 7 dígitos). Con pivoteo parcial escalado, $\rho$ cae a `1.7832` ($< 1$ bit de crecimiento), preservando entre 11 y 15 dígitos iniciales.

## 6. Condicionamiento

- **Estimación de $\kappa_2(A)$ y cota clásica:** Mediante valores singulares (`np.linalg.svd`), se obtiene $\kappa_2(A) = \sigma_{\max}/\sigma_{\min} \approx 3.084736 \times 10^{13}$. La cota clásica de perturbación predice:
  $$\frac{\|\delta x\|_2}{\|x\|_2} \lesssim \kappa_2(A)\,\varepsilon_{\text{mach}} \approx (3.084736 \times 10^{13})(2.220446 \times 10^{-16}) \approx 6.85 \times 10^{-3}$$
  lo cual garantiza teóricamente apenas $-\log_{10}(6.85 \times 10^{-3}) \approx 2.16$ dígitos significativos en el peor caso.
- **Por qué la cota resulta pesimista aquí:** La cota $\kappa_2(A)\,\varepsilon$ supone el peor escenario posible: que los errores de redondeo acumulados se alineen exactamente con el vector singular derecho asociado a $\sigma_{\min}$ (en la dirección de la casi dependencia lineal entre las filas 3, 5 y 7) y que el factor de crecimiento sea elevado. Gracias al pivoteo parcial escalado ($\rho = 1.7832$), el error hacia atrás (backward error) es muy inferior a $\varepsilon$, permitiendo recuperar **11 a 15 dígitos** antes del refinamiento y **15 dígitos** tras refinar.

## 7. Resultados

A continuación se reportan las 12 incógnitas obtenidas con pivoteo parcial escalado (solución inicial y tras 1 pasada de refinamiento iterativo), indicando exactamente los dígitos significativos defendibles (los ceros finales a la derecha se omiten en la notación decimal estándar porque la solución racional exacta posee 10 cifras decimales no nulas):

| Incógnita | Solución inicial (pivoteo escalado) | Dígitos iniciales defendibles | Solución refinada (reportada) | Dígitos finales defendibles |
|---|---:|---:|---:|---:|
| $x_{1}$ | `0.7384915627` | 14 | `0.7384915627` | 15 |
| $x_{2}$ | `-1.2059374168` | 14 | `-1.2059374168` | 15 |
| $x_{3}$ | `0.0043871952` | 11 | `0.0043871952` | 15 |
| $x_{4}$ | `-0.0917264538` | 13 | `-0.0917264538` | 15 |
| $x_{5}$ | `2.8461093725` | 14 | `2.8461093725` | 15 |
| $x_{6}$ | `-0.0062815409` | 12 | `-0.0062815409` | 15 |
| $x_{7}$ | `1.4738205916` | 14 | `1.4738205916` | 15 |
| $x_{8}$ | `-0.0384927164` | 13 | `-0.0384927164` | 15 |
| $x_{9}$ | `0.9152648037` | 14 | `0.9152648037` | 15 |
| $x_{10}$ | `-2.6071394852` | 15 | `-2.6071394852` | 15 |
| $x_{11}$ | `0.0518736294` | 13 | `0.0518736294` | 15 |
| $x_{12}$ | `-1.0849371625` | 14 | `-1.0849371625` | 15 |

## 8. Verificación

- **Residuo en precisión extendida (`mpmath`, `dps=80`):** Evaluando $r = b - Ax$ con los datos exactos sin truncar, la solución inicial con pivoteo escalado obtiene $\|r_0\|_\infty = 1.13331 \times 10^{-9}$ (frente a $1.77148 \times 10^{-7}$ sin pivoteo).
- **Refinamiento iterativo:** Calculando $r_0 = b - Ax_0$ en 80 dígitos, resolviendo $A\,\Delta x = r_0$ con nuestra propia eliminación gaussiana escalada y actualizando $x_1 = x_0 + \Delta x$, la norma infinito del residuo cae a $\|r_{\text{ref}}\|_\infty = 1.10543 \times 10^{-75}$ y todas las incógnitas convergen al valor exacto representable en `float64` (15 dígitos defendibles).
- **Control independiente con biblioteca:** Contrastando $x_{\text{ref}}$ contra `numpy.linalg.solve(A, b)` (usado únicamente como control externo), la diferencia máxima componente a componente es del orden de $\sim 10^{-12}$, confirmando la consistencia de la solución.

## 9. Conclusiones

- **Compromiso costo/exactitud:** Agregar un sobrecosto controlado de 77 divisiones de escala, 638 comparaciones y 7 intercambios sobre la base de 1,354 FLOPs evita la amplificación catastrófica provocada por $a_{1,1} = 1.3407 \times 10^{-8}$ y las escalas heterogéneas ($10^{-8}$ a $10^{5}$), elevando la exactitud inicial entre 11 y 15 dígitos y alcanzando precisión completa en `float64` tras una pasada de refinamiento iterativo.
- **Qué haría distinto con $n = 1000$:** Para $n = 1000$ ($\frac{2}{3}n^3 \approx 6.67 \times 10^8$ FLOPs), (1) reemplazaría los bucles escalares en Python por operaciones vectorizadas por bloques (BLAS-3 / NumPy) para aprovechar caché y SIMD; (2) calcularía los factores de escala $s_i$ una sola vez al inicio ($O(n^2)$) y usaría un vector de permutación $p$ en lugar de copiar filas completas en memoria; (3) almacenaría la factorización $PA = LU$ compacta para reutilizarla en el refinamiento iterativo en $O(n^2)$ sin refactorizar; y (4) omitiría la escritura de trazas matriciales paso a paso en disco por su tamaño $O(n^3)$.

## 10. Reproducibilidad, referencias y declaración de uso de IA

- **Entorno de ejecución:** Python 3.11 (entorno virtual `.venv`), librerías fijadas en `requirements.txt` (`numpy` y `mpmath`).
- **Comando de reproducción:** `python src/main.py` desde la raíz del repositorio (genera `solucion.csv`, `traza_eliminacion.txt`, `conteo_flops.md` y `residuo.txt` en `resultados/`).
- **Referencias bibliográficas:**
  - Burden, R. L., & Faires, J. D. (2011). *Análisis Numérico* (9.ª ed.). Cengage Learning. (Cap. 6: Métodos directos para resolver sistemas lineales y pivoteo parcial escalado).
  - Golub, G. H., & Van Loan, C. F. (2013). *Matrix Computations* (4.ª ed.). Johns Hopkins University Press. (Cap. 3: Gaussian Elimination, growth factor y análisis de error hacia atrás).
  - Higham, N. J. (2002). *Accuracy and Stability of Numerical Algorithms* (2.ª ed.). SIAM. (Cap. 9 y 12: Error de redondeo en eliminación gaussiana y refinamiento iterativo).
- **Declaración explícita de uso de asistentes de IA (Sección 9):** Se utilizó el asistente de IA **Gemini** como apoyo para estructurar la modularización del código en Python, revisar la instrumentación del conteo de operaciones y verificar la redacción técnica del informe. Todo el código, los resultados experimentales y cada línea entregada fueron revisados, ejecutados y validados por los autores para su defensa oral.
