# El pivote — Eliminación gaussiana con mínimo costo operacional y máxima exactitud de mantisa

## 1. Identificación

| Campo | Detalle |
|---|---|
| Estudiantes | Blas Emilio Rojas Díaz · Joaquín Alonso Hernández Castro |
| Asignatura | Computación Numérica (INF-321), Ingeniería Civil Informática, Universidad Católica del Maule |
| Docentes | Juan Manuel Cabezas Contardo · Pedro Pesante |
| Fecha | Talca, noviembre de 2026 |
| Repositorio | https://github.com/Blastjol/cn-parcial-gauss12-rojas-hernandez |

## 2. Planteamiento

Se resuelve $Ax=b$ con $A\in\mathbb{R}^{12\times12}$ densa (144 de 144 elementos no nulos, 9 cifras significativas por elemento) y $b$ entregado con 20 dígitos. Ambos se leen **sin truncar** desde `datos/` (`float(str)` redondea al `double` más cercano; la verificación usa los mismos textos en `mpmath`). Como control de integridad de los datos, $\det A = 1.3646\times10^{9}$ y $\kappa_2(A)\approx3.08\times10^{13}$ coinciden con los valores del enunciado.

Tres rasgos hacen costoso el problema en `float64` ($u=2^{-53}\approx1.11\times10^{-16}$, 53 bits de mantisa):

1. **Escalas de fila incompatibles.** La norma de las filas va de $2.1\times10^{-4}$ (F1) a $2.3\times10^{6}$ (F6). Restar a una fila pequeña un múltiplo de una grande alinea exponentes y descarta los bits bajos del operando pequeño.
2. **Pivote inicial minúsculo.** $a_{11}=1.3407\times10^{-8}$; $a_{21}/a_{11}\approx8.9\times10^{11}$ (doce órdenes) y $a_{61}/a_{11}\approx2.75\times10^{13}$.
3. **Casi dependencia lineal de la fila 7.** Un ajuste de mínimos cuadrados relativo da $F_7\approx4.1732\,F_3+582.69\,F_5$, con desviación relativa máxima por columna de $9\times10^{-3}$. El enunciado la señala como causa de $\kappa_2(A)\approx3\times10^{13}$ sin que $\det A$ se anule; la sección 6 mide cuánto de ese $\kappa_2$ proviene de la escala de las filas.

## 3. Estrategia de eliminación

Una única rutina (`eliminacion_gaussiana`, `src/gauss.py`) con tres políticas: `ninguno`, `parcial` y `escalado`. El orden es el clásico por columnas, $k=1,\dots,11$, y obtiene $PA=LU$.

- **Pivoteo parcial escalado.** Las escalas $s_i=\max_j|a_{ij}^{(0)}|$ se calculan **una vez** ($n(n-1)=132$ comparaciones). En el paso $k$ se elige la fila $p\ge k$ que maximiza $|a_{pk}|/s_p$ y se intercambian filas en $A$, $L$, $b$ y $s$.
- **Un multiplicador por fila.** $m_{ik}=a_{ik}/a_{kk}$ se calcula con **una** división, se guarda en $L$ y se reutiliza en toda la fila $i$ y en $b_i$. Recalcularlo por columna costaría ≈1750 flops.
- **No se toca la subdiagonal.** El barrido va de $j=k+1$ a $n$; el elemento $a_{ik}$ se fija en 0 sin operar. Así se anulan los 66 elementos bajo la diagonal sin gastar flops en ellos.
- **Una sola factorización, nunca la inversa.** $L$, $U$ y $P$ se reutilizan en el refinamiento, que cuesta 276 flops por pasada frente a 1354 de refactorizar.

**Por qué no excede el mínimo.** Con $m=n-k$ filas por debajo del pivote, el paso $k$ gasta $m$ divisiones + $m^2$ mult + $m^2$ restas (A) + $m$ mult + $m$ restas (b). Sumando $m=1..11$: $66+1012+132=1210$; la sustitución regresiva agrega $66+66+12=144$. Total $=1354$, dominado por $2n^3/3=1152$.

## 4. Conteo de operaciones

**Convención:** 1 flop = 1 operación aritmética elemental (`+ − × ÷`) en eliminación, actualización de $b$ y sustitución regresiva. Las divisiones de escalamiento, las comparaciones y los intercambios son **sobrecosto de estabilidad** y se reportan aparte. El conteo sale de contadores en `src/gauss.py` (impresos en `resultados/conteo_flops.md`); la forma cerrada solo se usa como verificación cruzada.

| Etapa | Mult/Div | Sum/Resta | FLOPs | Div. escala | Comparaciones | Intercambios |
|---|---:|---:|---:|---:|---:|---:|
| Previo (escalas $s_i$) | 0 | 0 | 0 | 0 | 132 | 0 |
| Paso 1 | 143 | 132 | 275 | 12 | 11 | 1 |
| Paso 2 | 120 | 110 | 230 | 11 | 10 | 0 |
| Paso 3 | 99 | 90 | 189 | 10 | 9 | 1 |
| Paso 4 | 80 | 72 | 152 | 9 | 8 | 1 |
| Paso 5 | 63 | 56 | 119 | 8 | 7 | 1 |
| Paso 6 | 48 | 42 | 90 | 7 | 6 | 1 |
| Paso 7 | 35 | 30 | 65 | 6 | 5 | 1 |
| Paso 8 | 24 | 20 | 44 | 5 | 4 | 1 |
| Paso 9 | 15 | 12 | 27 | 4 | 3 | 1 |
| Paso 10 | 8 | 6 | 14 | 3 | 2 | 0 |
| Paso 11 | 3 | 2 | 5 | 2 | 1 | 0 |
| Sustitución regresiva | 78 | 66 | 144 | 0 | 0 | 0 |
| **Total** | **716** | **638** | **1354** | **77** | **198** | **8** |

| Estrategia | FLOPs | vs. cota 1486 | Sobrecosto de estabilidad |
|---|---:|---:|---|
| Sin pivoteo | 1354 | −8.9 % | ninguno |
| Parcial | 1354 | −8.9 % | 66 comparaciones + 9 intercambios |
| **Parcial escalado (elegida)** | **1354** | **−8.9 %** | 77 div + 198 comp + 8 intercambios |

Los tres quedan dentro del ±10 % ([1337, 1635]). La diferencia con 1486 se debe a que aquí el barrido arranca en $k+1$ (sin tocar la columna anulada). El pivoteo no agrega un solo flop aritmético; su costo son ≈283 operaciones de otro tipo, $\mathcal{O}(n^2)$.

## 5. Análisis de mantisa y pivoteo

| Métrica | Sin pivoteo | Parcial | Parcial escalado |
|---|---:|---:|---:|
| Multiplicador máximo | $2.7526\times10^{13}$ ($2^{44.65}$) | $0.620$ | $2.4656\times10^{9}$ ($2^{31.20}$) |
| Factor de crecimiento $\rho$ | $2.6412\times10^{3}$ (11.37 bits) | $1.000$ | $1.7832$ (0.83 bits) |
| Bits del operando absorbidos (máx. en una resta) | 18.11 | 7.26 | 6.29 |
| Error relativo final $\|x_0-x\|_\infty/\|x\|_\infty$ | $2.08\times10^{-10}$ | $5.52\times10^{-13}$ | $2.36\times10^{-13}$ |
| Bits de mantisa perdidos, $\log_2(e/u)$ | 20.8 | 12.3 | 11.1 |
| Dígitos correctos (rango en $x_1..x_{12}$) | 6 a 9 | 9 a 12 | 9 a 12 |
| $\|r_0\|_\infty$ | $1.77\times10^{-7}$ | $5.4\times10^{-10}$ | $2.29\times10^{-10}$ |

**Por qué se rechaza $a_{11}=1.3407\times10^{-8}$.** Dividir por él da $m_{61}=2.75\times10^{13}=2^{44.65}$. Al hacer $a_{6j}-m_{61}a_{1j}$, el término $m_{61}a_{1j}$ domina por completo a $a_{6j}$: en el peor caso el operando original queda 18.1 bits por debajo del producto, de modo que esos bits de $a_{6j}$ se pierden de forma irrecuperable al redondear a 53 bits. Los elementos de la matriz crecen hasta $2.6\times10^{3}$ veces el máximo inicial (11.4 bits), y esa pérdida se propaga a todas las incógnitas: el error final es $2\times10^{-10}$ y $x_3,x_6$ quedan con solo 6 dígitos. Con pivoteo en el paso 1 se intercambian F1 y F8 (razón escalada $|a_{81}|/s_8=1$), y el crecimiento se mantiene en $\rho=1.78$.

**Lectura de los multiplicadores.** Que el pivoteo escalado deje un multiplicador de $2.5\times10^{9}$ no es un defecto: solo aparece cuando la fila $i$ es varios órdenes mayor que la fila pivote, y entonces $m\,a_{kj}$ es del tamaño de la propia fila $i$. Por eso los bits absorbidos son solo 6.29 y $\rho=1.78$. El multiplicador absoluto no mide el daño; $\rho$ y los bits absorbidos sí.

**Parcial vs. escalado, con honestidad.** En esta matriz el pivoteo parcial simple ya es casi tan bueno ($5.5\times10^{-13}$ vs. $2.4\times10^{-13}$): el escalado mejora el error inicial en un factor ≈2.3 a cambio de 77 divisiones y 132 comparaciones extra. Como ambos terminan en 15 dígitos tras refinar, la ventaja del escalado es pequeña pero real, y su costo es despreciable frente a los 1354 flops.

## 6. Condicionamiento

| Medida | Valor |
|---|---:|
| $\kappa_2(A)=\sigma_{\max}/\sigma_{\min}$ (SVD, solo diagnóstico) | $3.08\times10^{13}$ |
| $\kappa_1(A)$ / $\kappa_\infty(A)$ | $4.47\times10^{13}$ / $5.37\times10^{13}$ |
| $\kappa_2(D^{-1}A)$, filas equilibradas ($s_i=\max_j|a_{ij}|$) | $8.13\times10^{3}$ |

**Cota clásica.** $\|\delta x\|/\|x\|\lesssim\kappa_2\,\varepsilon=(3.08\times10^{13})(2.22\times10^{-16})\approx6.8\times10^{-3}$, que solo garantiza $-\log_{10}(6.8\times10^{-3})\approx2.2$ dígitos. Sin embargo, con pivoteo se obtienen 9–12 dígitos y con refinamiento 15.

**Por qué la cota es pesimista aquí.**
1. **Casi todo el $\kappa_2$ viene de la escala de las filas, no de la geometría.** Equilibrar las filas lo baja a $8.1\times10^{3}$ (casi diez órdenes de magnitud). El pivoteo escalado se comporta como eliminar sobre $D^{-1}A$, cuya cota $8.1\times10^{3}\cdot\varepsilon\approx1.8\times10^{-12}$ (≈11.7 dígitos) es consistente con el error observado $2.4\times10^{-13}$.
2. La cota supone un error de retroceso normwise de tamaño $\varepsilon\|A\|$, mientras que con pivoteo el error es componentwise, del orden de $\varepsilon|A|$; las filas pequeñas no se contaminan con la norma de las grandes.
3. La cota usa el peor caso de alineación del error con el vector singular de $\sigma_{\min}$; el error real rara vez lo alcanza.

La casi dependencia entre F3, F5 y F7 es detectable (sección 2), pero es el residuo, no la cota, quien decide el número de dígitos. Los datos ($A$ con 9 cifras exactas, $b$ con 20 dígitos) se tratan como exactos, así que no hay error de entrada que $\kappa$ pueda amplificar.

## 7. Resultados

Solución con pivoteo parcial escalado + 1 pasada de refinamiento. Dígitos correctos $d_i$: mayor $d$ con $|e_{\text{rel}}|\le0.5\times10^{-d}$ frente a la referencia de 80 dígitos, con tope 15 (límite de `float64`). Se reportan exactamente $d_i$ cifras.

| Incógnita | Valor reportado | Dígitos afirmados | Valor inicial (antes de refinar) | Dígitos iniciales |
|---|---:|---:|---:|---:|
| $x_1$ | `0.738491562700000` | 15 | `0.738491562700` | 12 |
| $x_2$ | `-1.20593741680000` | 15 | `-1.2059374168` | 11 |
| $x_3$ | `0.00438719520000000` | 15 | `0.00438719520` | 9 |
| $x_4$ | `-0.0917264538000000` | 15 | `-0.091726453800` | 11 |
| $x_5$ | `2.84610937250000` | 15 | `2.84610937250` | 12 |
| $x_6$ | `-0.00628154090000000` | 15 | `-0.00628154090` | 9 |
| $x_7$ | `1.47382059160000` | 15 | `1.47382059160` | 12 |
| $x_8$ | `-0.0384927164000000` | 15 | `-0.038492716400` | 11 |
| $x_9$ | `0.915264803700000` | 15 | `0.915264803700` | 12 |
| $x_{10}$ | `-2.60713948520000` | 15 | `-2.60713948520` | 12 |
| $x_{11}$ | `0.0518736294000000` | 15 | `0.051873629400` | 11 |
| $x_{12}$ | `-1.08493716250000` | 15 | `-1.08493716250` | 12 |

La solución exacta tiene a lo sumo 10 decimales; por eso, pasados esos decimales, los dígitos reportados son ceros, y son correctos. Los 12 valores superan el umbral de 8 dígitos. Los valores con 17 cifras están en `resultados/solucion.csv`.

## 8. Verificación

- **Residuo en precisión extendida** (`mpmath`, 80 dígitos, datos exactos): $\|r_0\|_\infty=2.29\times10^{-10}$ con pivoteo escalado ($1.77\times10^{-7}$ sin pivoteo). El vector completo está en `resultados/residuo.txt`.
- **Refinamiento iterativo** (residuo extendido, corrección con $P,L,U$ ya calculados, 276 flops por pasada):

| Variante | Pasada 1 | Pasada 2 | Pasada 3 | Pasada 4 |
|---|---:|---:|---:|---:|
| Iterado en `float64` (error relativo) | $6.0\times10^{-17}$ | $6.0\times10^{-17}$ | $6.0\times10^{-17}$ | — |
| Iterado en precisión extendida (error relativo) | $6.2\times10^{-26}$ | $1.6\times10^{-38}$ | $4.2\times10^{-51}$ | $1.1\times10^{-63}$ |

El error relativo cae de $2.4\times10^{-13}$ a $6.0\times10^{-17}$ en una pasada, y luego no baja más porque $x$ está guardado en `float64` (media unidad de redondeo es $1.1\times10^{-16}$). Con el iterado en precisión extendida se ve que cada pasada gana ≈12.6 dígitos, que es la calidad de la factorización. **El residuo no mejora** en `float64` ($2.29\times10^{-10}\to2.60\times10^{-10}$): su piso es $u\,|A||x|\approx10^{-10}$. Lo que mejora es el error en $x$, no el residuo, y esa diferencia es la lección de un sistema mal condicionado.
- **Control con biblioteca:** `numpy.linalg.solve` (solo control) da un error relativo de $5.2\times10^{-13}$ y 9–12 dígitos, equivalente a la corrida propia sin refinar; la solución refinada difiere de la de numpy en ≤ $1.5\times10^{-12}$.

## 9. Conclusiones

- **Compromiso elegido.** La eliminación base (1354 flops, −8.9 % de la cota) es la misma para las tres políticas. Elegir el pivoteo parcial escalado agrega solo 77 divisiones, 198 comparaciones y 8 intercambios, evita el crecimiento de $2.6\times10^{3}$ del caso sin pivoteo, y deja 9–12 dígitos. El refinamiento sube a 15 dígitos por solo 276 flops, sin refactorizar. No hay un compromiso real entre costo y exactitud: la estabilidad cuesta $\mathcal{O}(n^2)$, y el mínimo $\mathcal{O}(n^3)$ no se toca.
- **Qué dice la comparación.** La cota $\kappa\varepsilon$ prometía 2 dígitos y se obtuvieron 15; no sirve como estimación, solo como garantía. El diagnóstico útil fue equilibrar las filas ($\kappa_2$: $3\times10^{13}\to8\times10^{3}$).
- **Con $n=1000$** ($\tfrac23n^3\approx6.7\times10^{8}$ flops): (1) LU por bloques con BLAS-3 en lugar de bucles escalares; (2) pivoteo parcial simple, o equilibrar filas una vez en $\mathcal{O}(n^2)$ (el escalado dinámico cuesta divisiones extra por paso), con un vector de permutación en lugar de copiar filas; (3) factorización in situ (L y U en la misma matriz) y sin trazas por paso, que serían $\mathcal{O}(n^3)$ en disco; (4) estimar $\kappa_1$ con el estimador de Hager–Higham sobre los factores en $\mathcal{O}(n^2)$ en lugar de una SVD $\mathcal{O}(n^3)$; (5) refinamiento con factorización en baja precisión y residuo en alta.

## 10. Reproducibilidad

- **Entorno:** Python 3.12.3 (se exige ≥ 3.10), `numpy==2.2.6`, `mpmath==1.3.0` (`requirements.txt`). Probado en un entorno virtual limpio.
- **Comando:** `python src/main.py` desde la raíz; regenera `resultados/conteo_flops.md`, `traza_eliminacion.txt`, `solucion.csv` y `residuo.txt` sin intervención manual.
- **Variaciones esperables:** $\kappa_2$ y el error de `numpy.linalg.solve` dependen de la versión de LAPACK en el 4.º–5.º dígito; por eso se citan con 3 cifras. Flops, multiplicadores, $\rho$ y los errores del código propio son deterministas.
- **Fuentes:** Burden y Faires, *Análisis Numérico* (9.ª ed., 2011), cap. 6; Golub y Van Loan, *Matrix Computations* (4.ª ed., 2013), cap. 3; Higham, *Accuracy and Stability of Numerical Algorithms* (2.ª ed., 2002), caps. 9 y 12.
- **Uso de asistentes de IA.** Se usó **Gemini** para la estructura modular inicial del código y un primer borrador del informe. Se usó **Claude** (Anthropic) para revisar el repositorio contra la rúbrica, corregir cifras del informe que no coincidían con la ejecución, ampliar `gauss.py`, `conteo.py`, `verificacion.py` y `main.py`, y redactar nuevamente `informe.md` y `readme.md`. Todas las cifras de este informe provienen de la ejecución de `python src/main.py`; los autores deben poder explicar cada línea en la defensa oral.
