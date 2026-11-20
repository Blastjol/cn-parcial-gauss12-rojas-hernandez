# Eliminación Gaussiana — Sistema 12x12

Universidad Católica del Maule, Computación Numérica.

Integrantes: Blas Emilio Rojas Díaz y Joaquín Alonso Hernández Castro.

### Qué hace el proyecto

Resuelve un sistema lineal denso de $12 \times 12$ (`datos/matriz_a.csv` y `datos/vector_b.csv`) mediante eliminación gaussiana en aritmética `float64`, comparando eliminación **sin pivoteo** y **con pivoteo parcial escalado**, con conteo instrumentado de FLOPs, verificación y refinamiento iterativo en precisión extendida (`mpmath`) y contraste con `numpy.linalg.solve`. El análisis completo está en [informe.md](informe.md).

### Estructura

```text
cn-parcial-gauss12-rojas-hernandez/
├── README.md
├── informe.md
├── requirements.txt
├── datos/
│   ├── matriz_a.csv
│   └── vector_b.csv
├── resultados/
│   ├── solucion.csv
│   ├── traza_eliminacion.txt
│   ├── conteo_flops.md
│   └── residuo.txt
└── src/
    ├── __init__.py
    ├── conteo.py
    ├── gauss.py
    ├── verificacion.py
    └── main.py
```

### Scripts (`src/`)

- `main.py`: es el que se ejecuta. Carga los datos, corre ambas estrategias, verifica la solución y genera los archivos en `resultados/`.
- `gauss.py`: implementación de eliminación gaussiana y sustitución regresiva, sin pivoteo y con pivoteo parcial escalado.
- `conteo.py`: contador de operaciones aritméticas, comparaciones e intercambios, y generación del reporte de FLOPs.
- `verificacion.py`: cálculo de residuo y refinamiento iterativo en alta precisión (`mpmath`), estimación de dígitos defendibles y comparación con `numpy.linalg.solve`.

### Cómo ejecutarlo

Clonar el repositorio y crear el entorno virtual:

```bash
git clone [https://github.com/Blastjol/cn-parcial-gauss12-rojas-hernandez.git](https://github.com/Blastjol/cn-parcial-gauss12-rojas-hernandez.git)
cd cn-parcial-gauss12-rojas-hernandez
python -m venv .venv
```

Activarlo:

```bash
# Windows (PowerShell)
.\.venv\Scripts\Activate.ps1

# Linux / macOS
source .venv/bin/activate
```

Instalar dependencias y ejecutar:

```bash
pip install -r requirements.txt
python -m src.main
```

### Resultados

- `resultados/solucion.csv`: valores de las incógnitas $x_1$ a $x_{12}$, dígitos significativos defendibles, comparación entre estrategias y control con `numpy.linalg.solve`.
- `resultados/traza_eliminacion.txt`: traza paso a paso de la eliminación sin pivoteo y con pivoteo parcial escalado.
- `resultados/conteo_flops.md`: desglose de FLOPs, comparaciones e intercambios por etapa y total.
- `resultados/residuo.txt`: normas del residuo inicial y refinado, corrección $\Delta x$, número de condición $\kappa_2(A)$ y contraste con `numpy.linalg.solve`.

---