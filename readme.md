# Eliminación Gaussiana — Sistema 12x12

Resuelve un sistema lineal denso de $12 \times 12$ (`datos/matriz_a.csv` y `datos/vector_b.csv`) mediante eliminación gaussiana en aritmética `float64`, comparando eliminación **sin pivoteo** y **con pivoteo parcial escalado**, con conteo instrumentado de FLOPs, verificación y refinamiento iterativo en precisión extendida (`mpmath`) y contraste con `numpy.linalg.solve`.

### Autoría y asignatura

- **Integrantes:** Blas Emilio Rojas Díaz y Joaquín Alonso Hernández Castro.
- **Asignatura:** Computación Numérica (INF-321) — Ingeniería Civil Informática, Universidad Católica del Maule.
- **Docentes:** Juan Manuel Cabezas Contardo · Pedro Pesante.

### Requisitos

- **Python:** `3.10` o superior.
- **Dependencias:** `numpy`, `mpmath` (versiones fijadas en `requirements.txt`).

Crear y activar entorno virtual:

```bash
git clone [https://github.com/Blastjol/cn-parcial-gauss12-rojas-hernandez.git](https://github.com/Blastjol/cn-parcial-gauss12-rojas-hernandez.git)
cd cn-parcial-gauss12-rojas-hernandez
python -m venv .venv

# Windows (PowerShell)
.\.venv\Scripts\Activate.ps1

# Linux / macOS
source .venv/bin/activate