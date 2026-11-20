# Conteo instrumentado de operaciones (FLOPs)

## Convención de conteo
- **1 flop = 1 operación aritmética elemental** (`+`, `-`, `*`, `/`) ejecutada en el código.
- Se contabilizan por separado las **multiplicaciones/divisiones** y las **sumas/restas**.
- Las **comparaciones** e **intercambios** asociados al pivoteo se registran aparte como sobrecosto de estabilidad.
- Se incluyen la actualización del vector `b` y la sustitución regresiva.
- Referencia teórica $2n^3/3$ para $n=12$: `1152.00` (referencia práctica del enunciado: ~1486 flops + ~66 comparaciones con pivoteo).

## Estrategia: Sin pivoteo

| Etapa / Paso | Mult / Div | Sum / Resta | FLOPs paso | Comparaciones | Intercambios |
|---|---:|---:|---:|---:|---:|
| Paso 1 | 143 | 132 | 275 | 0 | 0 |
| Paso 2 | 120 | 110 | 230 | 0 | 0 |
| Paso 3 | 99 | 90 | 189 | 0 | 0 |
| Paso 4 | 80 | 72 | 152 | 0 | 0 |
| Paso 5 | 63 | 56 | 119 | 0 | 0 |
| Paso 6 | 48 | 42 | 90 | 0 | 0 |
| Paso 7 | 35 | 30 | 65 | 0 | 0 |
| Paso 8 | 24 | 20 | 44 | 0 | 0 |
| Paso 9 | 15 | 12 | 27 | 0 | 0 |
| Paso 10 | 8 | 6 | 14 | 0 | 0 |
| Paso 11 | 3 | 2 | 5 | 0 | 0 |
| Sustitución regresiva | 78 | 66 | 144 | 0 | 0 |
| **Total** | **716** | **638** | **1354** | **0** | **0** |

## Estrategia: Con pivoteo parcial escalado

| Etapa / Paso | Mult / Div | Sum / Resta | FLOPs paso | Comparaciones | Intercambios |
|---|---:|---:|---:|---:|---:|
| Paso 1 | 155 | 132 | 287 | 143 | 1 |
| Paso 2 | 131 | 110 | 241 | 120 | 0 |
| Paso 3 | 109 | 90 | 199 | 99 | 1 |
| Paso 4 | 89 | 72 | 161 | 80 | 1 |
| Paso 5 | 71 | 56 | 127 | 63 | 1 |
| Paso 6 | 55 | 42 | 97 | 48 | 1 |
| Paso 7 | 41 | 30 | 71 | 35 | 1 |
| Paso 8 | 29 | 20 | 49 | 24 | 1 |
| Paso 9 | 19 | 12 | 31 | 15 | 0 |
| Paso 10 | 11 | 6 | 17 | 8 | 0 |
| Paso 11 | 5 | 2 | 7 | 3 | 0 |
| Sustitución regresiva | 78 | 66 | 144 | 0 | 0 |
| **Total** | **793** | **638** | **1431** | **638** | **7** |
