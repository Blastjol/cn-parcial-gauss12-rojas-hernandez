# Evaluacion 3 El pivote

Resuelve un sistema lineal denso de 12 incógnitas con eliminación gaussiana propia (sin pivoteo y con pivoteo parcial escalado), con conteo instrumentado de flops, residuo en precisión extendida y refinamiento iterativo.

## Autoría y asignatura

- **Estudiantes:** Blas Emilio Rojas Díaz y Joaquín Alonso Hernández Castro.
- **Asignatura:** Computación Numérica (INF-321), Ingeniería Civil Informática, Universidad Católica del Maule. Prueba parcial «El pivote», noviembre de 2026.
- **Docentes:** Juan Manuel Cabezas Contardo · Pedro Pesante.

## Requisitos

- Python **3.10 o superior** (probado con 3.12).
- Dependencias con versión fija: `numpy==2.2.6`, `mpmath==1.3.0`.

```bash
git clone https://github.com/Blastjol/cn-parcial-gauss12-rojas-hernandez.git
cd cn-parcial-gauss12-rojas-hernandez
python -m venv .venv
source .venv/bin/activate          # Windows (PowerShell): .\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

## Ejecución

```bash
python src/main.py
```

Un solo comando, sin intervención manual. Regenera en `resultados/`: `solucion.csv`, `traza_eliminacion.txt`, `conteo_flops.md` y `residuo.txt`, y escribe un resumen por consola (tarda unos segundos).

## Estructura de carpetas

```text
datos/       A (12x12) y b con todas sus cifras, sin truncar
src/         gauss.py (eliminación), conteo.py (flops), verificacion.py (residuo y refinamiento), main.py
resultados/  salidas generadas por main.py
informe.md   informe técnico y análisis de error
```

## Resultado resumido

Pivoteo parcial escalado + 1 pasada de refinamiento (10 cifras significativas; 15 dígitos correctos en cada incógnita):

| Incógnita | Valor |
|---|---:|
| $x_{1}$ | `0.7384915627` |
| $x_{2}$ | `-1.205937417` |
| $x_{3}$ | `0.004387195200` |
| $x_{4}$ | `-0.09172645380` |
| $x_{5}$ | `2.846109372` |
| $x_{6}$ | `-0.006281540900` |
| $x_{7}$ | `1.473820592` |
| $x_{8}$ | `-0.03849271640` |
| $x_{9}$ | `0.9152648037` |
| $x_{10}$ | `-2.607139485` |
| $x_{11}$ | `0.05187362940` |
| $x_{12}$ | `-1.084937162` |

| Magnitud | Valor |
|---|---:|
| Flops aritméticos (n = 12) | 1354 |
| Sobrecosto del pivoteo escalado | 77 div + 198 comp + 8 intercambios |
| ‖r‖∞ inicial / refinado (precisión extendida) | 2.29e-10 / 2.60e-10 |

## Informe

Análisis completo en [informe.md](informe.md).

## Licencia

Trabajo académico de uso educativo para INF-321. No se concede licencia de reutilización comercial.
