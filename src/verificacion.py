from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple
import mpmath as mp
import numpy as np

from src.gauss import eliminacion_gaussiana


@dataclass
class ResultadoVerificacion:
    residuo_inicial: List[mp.mpf]
    norma_inf_inicial: mp.mpf
    correccion_dx: np.ndarray
    solucion_refinada: np.ndarray
    residuo_refinado: List[mp.mpf]
    norma_inf_refinada: mp.mpf
    solucion_numpy: np.ndarray
    referencia_alta_precision: List[mp.mpf]
    digitos_correctos_inicial: List[int]
    digitos_correctos_refinada: List[int]
    pasadas_refinamiento: int


def _cargar_datos_mp(
    ruta_matriz_a: str | Path,
    ruta_vector_b: str | Path,
    dps: int = 80,
) -> Tuple[List[List[mp.mpf]], List[mp.mpf]]:
    mp.mp.dps = dps
    with open(ruta_matriz_a, mode="r", encoding="utf-8", newline="") as f_a:
        A_mp = [
            [mp.mpf(celda.strip()) for celda in fila if celda.strip() != ""]
            for fila in csv.reader(f_a)
            if any(celda.strip() != "" for celda in fila)
        ]
    with open(ruta_vector_b, mode="r", encoding="utf-8", newline="") as f_b:
        b_mp = [
            mp.mpf(fila[0].strip())
            for fila in csv.reader(f_b)
            if fila and fila[0].strip() != ""
        ]
    return A_mp, b_mp


def calcular_residuo_mp(
    A_mp: List[List[mp.mpf]],
    b_mp: List[mp.mpf],
    x_vec: np.ndarray | List[mp.mpf],
    dps: int = 80,
) -> Tuple[List[mp.mpf], mp.mpf]:
    """Calcula r = b - Ax y ||r||_inf en precisión extendida mpmath."""
    mp.mp.dps = dps
    n = len(b_mp)
    x_mp = [mp.mpf(str(float(val))) if not isinstance(val, mp.mpf) else val for val in x_vec]
    residuo: List[mp.mpf] = []
    max_abs = mp.mpf("0")

    for i in range(n):
        ax_i = mp.mpf("0")
        for j in range(n):
            ax_i += A_mp[i][j] * x_mp[j]
        r_i = b_mp[i] - ax_i
        residuo.append(r_i)
        abs_r_i = abs(r_i)
        if abs_r_i > max_abs:
            max_abs = abs_r_i

    return residuo, max_abs


def _estimar_referencia_mp(
    A_mp: List[List[mp.mpf]],
    b_mp: List[mp.mpf],
    x_base: np.ndarray,
    usar_pivoteo_escalado: bool = True,
    iteraciones: int = 4,
    dps: int = 80,
) -> List[mp.mpf]:
    """Construye una solución de referencia refinada en alta precisión."""
    mp.mp.dps = dps
    A_float = np.array([[float(val) for val in fila] for fila in A_mp], dtype=np.float64)
    x_actual = np.array(x_base, dtype=np.float64, copy=True)

    for _ in range(iteraciones):
        residuo_mp, _ = calcular_residuo_mp(A_mp, b_mp, x_actual, dps=dps)
        r_float = np.array([float(val) for val in residuo_mp], dtype=np.float64)
        res_corr = eliminacion_gaussiana(
            A_float,
            r_float,
            usar_pivoteo_escalado=usar_pivoteo_escalado,
            contador=None,
        )
        x_actual = x_actual + res_corr.solucion_x

    return [mp.mpf(str(val)) for val in x_actual]


def estimar_digitos_correctos(
    x_aprox: np.ndarray,
    x_referencia: List[mp.mpf],
    dps: int = 80,
) -> List[int]:
    """Estima dígitos significativos defendibles por componente frente a la referencia."""
    mp.mp.dps = dps
    digitos: List[int] = []
    for idx, ref in enumerate(x_referencia):
        aprox_mp = mp.mpf(str(float(x_aprox[idx])))
        if ref == 0:
            digitos.append(0)
            continue
        error_rel = abs(aprox_mp - ref) / abs(ref)
        if error_rel == 0:
            d = 15
        else:
            d = int(math.floor(-mp.log10(error_rel)))
        d_defendible = max(0, min(15, d))
        digitos.append(d_defendible)
    return digitos


def verificar_solucion(
    ruta_matriz_a: str | Path,
    ruta_vector_b: str | Path,
    x_inicial: np.ndarray,
    usar_pivoteo_escalado: bool = True,
    pasadas_refinamiento: int = 1,
    dps: int = 80,
) -> ResultadoVerificacion:
    """Ejecuta residuo extendido, refinamiento iterativo, control numpy y dígitos correctos."""
    if pasadas_refinamiento < 1:
        raise ValueError("Se requiere al menos una pasada de refinamiento iterativo.")

    A_mp, b_mp = _cargar_datos_mp(ruta_matriz_a, ruta_vector_b, dps=dps)
    A_float = np.array([[float(val) for val in fila] for fila in A_mp], dtype=np.float64)
    b_float = np.array([float(val) for val in b_mp], dtype=np.float64)

    res_ini, norm_ini = calcular_residuo_mp(A_mp, b_mp, x_inicial, dps=dps)

    x_actual = np.array(x_inicial, dtype=np.float64, copy=True)
    ultima_dx = np.zeros_like(x_actual)

    for _ in range(pasadas_refinamiento):
        r_actual_mp, _ = calcular_residuo_mp(A_mp, b_mp, x_actual, dps=dps)
        r_actual_float = np.array([float(val) for val in r_actual_mp], dtype=np.float64)
        corr = eliminacion_gaussiana(
            A_float,
            r_actual_float,
            usar_pivoteo_escalado=usar_pivoteo_escalado,
            contador=None,
        )
        ultima_dx = corr.solucion_x
        x_actual = x_actual + ultima_dx

    res_ref, norm_ref = calcular_residuo_mp(A_mp, b_mp, x_actual, dps=dps)

    # Control independiente mediante rutina de biblioteca
    sol_numpy = np.linalg.solve(A_float, b_float)

    # Referencia de alta precisión para defender dígitos significativos
    ref_mp = _estimar_referencia_mp(
        A_mp,
        b_mp,
        x_actual,
        usar_pivoteo_escalado=usar_pivoteo_escalado,
        iteraciones=3,
        dps=dps,
    )

    dig_ini = estimar_digitos_correctos(x_inicial, ref_mp, dps=dps)
    dig_ref = estimar_digitos_correctos(x_actual, ref_mp, dps=dps)

    return ResultadoVerificacion(
        residuo_inicial=res_ini,
        norma_inf_inicial=norm_ini,
        correccion_dx=ultima_dx,
        solucion_refinada=x_actual,
        residuo_refinado=res_ref,
        norma_inf_refinada=norm_ref,
        solucion_numpy=sol_numpy,
        referencia_alta_precision=ref_mp,
        digitos_correctos_inicial=dig_ini,
        digitos_correctos_refinada=dig_ref,
        pasadas_refinamiento=pasadas_refinamiento,
    )