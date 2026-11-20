# Conteo instrumentado de operaciones (FLOPs)

## Convención de conteo
- **1 flop = 1 operación aritmética elemental** (`+`, `-`, `*`, `/`) en eliminación, actualización de `b` y sustitución regresiva.
- Cada multiplicador $m_{ik} = a_{ik}/a_{kk}$ se calcula **una sola vez por fila** y el barrido va desde $j = k+1$ hasta $n$ (sin operar sobre los 66 elementos subdiagonales anulados).
- El sobrecosto del pivoteo parcial escalado (**divisiones de escalamiento** $|a_{ik}|/s_i$, **comparaciones** de búsqueda e **intercambios** de filas) se contabiliza por separado.
- Referencia asintótica dominante $2n^3/3$ para $n=12$: `1152.00` FLOPs; cota de referencia del enunciado (incluyendo `b` y sustitución regresiva): `~1486` FLOPs (rango aceptado $\pm 10\%$: `[1337, 1635]`).

## Estrategia: Sin pivoteo

| Etapa / Paso | Mult / Div | Sum / Resta | FLOPs eliminación/sust. | Div. escalamiento | Comparaciones | Intercambios |
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

## Estrategia: Con pivoteo parcial escalado

| Etapa / Paso | Mult / Div | Sum / Resta | FLOPs eliminación/sust. | Div. escalamiento | Comparaciones | Intercambios |
|---|---:|---:|---:|---:|---:|---:|
| Paso 1 | 143 | 132 | 275 | 12 | 143 | 1 |
| Paso 2 | 120 | 110 | 230 | 11 | 120 | 0 |
| Paso 3 | 99 | 90 | 189 | 10 | 99 | 1 |
| Paso 4 | 80 | 72 | 152 | 9 | 80 | 1 |
| Paso 5 | 63 | 56 | 119 | 8 | 63 | 1 |
| Paso 6 | 48 | 42 | 90 | 7 | 48 | 1 |
| Paso 7 | 35 | 30 | 65 | 6 | 35 | 1 |
| Paso 8 | 24 | 20 | 44 | 5 | 24 | 1 |
| Paso 9 | 15 | 12 | 27 | 4 | 15 | 0 |
| Paso 10 | 8 | 6 | 14 | 3 | 8 | 0 |
| Paso 11 | 3 | 2 | 5 | 2 | 3 | 0 |
| Sustitución regresiva | 78 | 66 | 144 | 0 | 0 | 0 |
| **Total** | **716** | **638** | **1354** | **77** | **638** | **7** |
