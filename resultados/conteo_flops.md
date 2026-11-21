# Conteo instrumentado de operaciones (FLOPs)

Generado por `python src/main.py`. Todas las cifras provienen de contadores
incrementados dentro de `src/gauss.py` en cada operación ejecutada.

## Convención de conteo

- **1 flop = 1 operación aritmética elemental** (`+`, `-`, `*`, `/`) en la eliminación, la actualización de `b` y la sustitución regresiva.
- Cada multiplicador $m_{ik} = a_{ik}/a_{kk}$ se calcula **una sola vez** (1 división) y se reutiliza en toda la fila $i$ y en $b_i$.
- El barrido del paso $k$ recorre las columnas $j = k+1,\dots,n$: **no se opera sobre la subdiagonal**, que se fija en 0 sin costo (66 elementos).
- Sobrecosto de estabilidad (divisiones $|a_{ik}|/s_i$, comparaciones de búsqueda, cálculo único de escalas $s_i$ e intercambios de filas) **se contabiliza aparte** y no se suma a los flops aritméticos.
- Los diagnósticos (factor de crecimiento, bits absorbidos, traza) no se cuentan: no forman parte del algoritmo.
- Cota de referencia del enunciado: ≈ 1486 flops (eliminación densa mínima, con `b` y sustitución regresiva); rango ±10 % = [1337, 1635]. Término cúbico dominante $2n^3/3$ = 1152.

## Estrategia: Sin pivoteo

| Etapa | Mult/Div | Sum/Resta | FLOPs aritméticos | Div. escalamiento | Comparaciones | Intercambios |
|---|---:|---:|---:|---:|---:|---:|
| Paso 1 | 143 | 132 | 275 | 0 | 0 | 0 |
| Paso 2 | 120 | 110 | 230 | 0 | 0 | 0 |
| Paso 3 | 99 | 90 | 189 | 0 | 0 | 0 |
| Paso 4 | 80 | 72 | 152 | 0 | 0 | 0 |
| Paso 5 | 63 | 56 | 119 | 0 | 0 | 0 |
| Paso 6 | 48 | 42 | 90 | 0 | 0 | 0 |
| Paso 7 | 35 | 30 | 65 | 0 | 0 | 0 |
| Paso 8 | 24 | 20 | 44 | 0 | 0 | 0 |
| Paso 9 | 15 | 12 | 27 | 0 | 0 | 0 |
| Paso 10 | 8 | 6 | 14 | 0 | 0 | 0 |
| Paso 11 | 3 | 2 | 5 | 0 | 0 | 0 |
| Sustitución regresiva | 78 | 66 | 144 | 0 | 0 | 0 |
| **Total** | **716** | **638** | **1354** | **0** | **0** | **0** |

## Estrategia: Con pivoteo parcial

| Etapa | Mult/Div | Sum/Resta | FLOPs aritméticos | Div. escalamiento | Comparaciones | Intercambios |
|---|---:|---:|---:|---:|---:|---:|
| Paso 1 | 143 | 132 | 275 | 0 | 11 | 1 |
| Paso 2 | 120 | 110 | 230 | 0 | 10 | 1 |
| Paso 3 | 99 | 90 | 189 | 0 | 9 | 1 |
| Paso 4 | 80 | 72 | 152 | 0 | 8 | 1 |
| Paso 5 | 63 | 56 | 119 | 0 | 7 | 1 |
| Paso 6 | 48 | 42 | 90 | 0 | 6 | 1 |
| Paso 7 | 35 | 30 | 65 | 0 | 5 | 0 |
| Paso 8 | 24 | 20 | 44 | 0 | 4 | 1 |
| Paso 9 | 15 | 12 | 27 | 0 | 3 | 1 |
| Paso 10 | 8 | 6 | 14 | 0 | 2 | 0 |
| Paso 11 | 3 | 2 | 5 | 0 | 1 | 1 |
| Sustitución regresiva | 78 | 66 | 144 | 0 | 0 | 0 |
| **Total** | **716** | **638** | **1354** | **0** | **66** | **9** |

## Estrategia: Con pivoteo parcial escalado

| Etapa | Mult/Div | Sum/Resta | FLOPs aritméticos | Div. escalamiento | Comparaciones | Intercambios |
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

## Resumen frente a la cota

| Estrategia | FLOPs aritméticos | Desviación vs 1486 | Dentro de ±10 % | Sobrecosto de estabilidad |
|---|---:|---:|:---:|---|
| Sin pivoteo | 1354 | -8.9 % | sí | ninguno |
| Con pivoteo parcial | 1354 | -8.9 % | sí | 0 div + 66 comp + 9 intercambios |
| Con pivoteo parcial escalado | 1354 | -8.9 % | sí | 77 div + 198 comp + 8 intercambios |

## Verificación cruzada con forma cerrada

Forma cerrada para $n=12$: 716 mult/div + 638 sum/resta = **1354** flops. Coincide con el contador instrumentado (la fórmula solo valida; no genera el conteo).

Cada pasada de refinamiento iterativo reutiliza $L$ y $U$ (sin refactorizar) y cuesta **276** flops ($\mathcal{O}(n^2)$), frente a 1354 de una nueva eliminación.

## Alternativas descartadas (referencia del enunciado)

- Multiplicadores recalculados por columna: ≈ 1750 flops (exceso penalizado).
- Vía matriz inversa: ≈ 7000 flops o más (descalifica C2). No se calcula la inversa en ningún punto.
