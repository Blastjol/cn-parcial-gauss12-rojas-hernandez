from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

import mpmath as mp
import numpy as np

from src.conteo import ContadorOperaciones
from src.gauss import ResultadoEliminacion, resolver_con_lu

DPS = 80
MAX_DIGITOS = 15  # tope defendible en float64


@dataclass
class Iteracion:
    pasada: int
    norma_inf_residuo: mp.mpf
    error_rel_inf: float


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
    digitos_correctos_numpy: List[int]
    error_relativo_inf_inicial: float
    error_relativo_inf_refinado: float
    error_relativo_2_inicial: float
    error_relativo_2_refinado: float
    error_relativo_inf_numpy: float
    historial_float: List[Iteracion]
    historial_extendido: List[Iteracion]
    flops_por_pasada: int
    pasadas_refinamiento: int


def cargar_datos_mp(ruta_a, ruta_b, dps: int = DPS) -> Tuple[List[List[mp.mpf]], List[mp.mpf]]:
    mp.mp.dps = dps
    with open(ruta_a, encoding="utf-8", newline="") as f:
        A = [[mp.mpf(c.strip()) for c in fila if c.strip()] for fila in csv.reader(f) if any(c.strip() for c in fila)]
    with open(ruta_b, encoding="utf-8", newline="") as f:
        b = [mp.mpf(fila[0].strip()) for fila in csv.reader(f) if fila and fila[0].strip()]
    return A, b


def calcular_residuo_mp(A, b, x, dps: int = DPS) -> Tuple[List[mp.mpf], mp.mpf]:
    """r = b - A x con x convertido SIN pérdida (float -> mpf es exacto) y A, b exactos."""
    mp.mp.dps = dps
    n = len(b)
    xm = [v if isinstance(v, mp.mpf) else mp.mpf(float(v)) for v in x]
    r = []
    for i in range(n):
        s = mp.mpf(0)
        for j in range(n):
            s += A[i][j] * xm[j]
        r.append(b[i] - s)
    return r, max(abs(v) for v in r)


def referencia_mp(A, b, dps: int = DPS) -> List[mp.mpf]:
    """Eliminación propia con pivoteo escalado en dps dígitos (no usa mp.lu_solve)."""
    mp.mp.dps = dps
    n = len(b)
    M = [fila[:] for fila in A]
    rhs = b[:]
    s = [max(abs(v) for v in fila) for fila in M]
    for k in range(n - 1):
        p = max(range(k, n), key=lambda i: abs(M[i][k]) / s[i])
        if p != k:
            M[k], M[p] = M[p], M[k]
            rhs[k], rhs[p] = rhs[p], rhs[k]
            s[k], s[p] = s[p], s[k]
        for i in range(k + 1, n):
            m = M[i][k] / M[k][k]
            M[i][k] = mp.mpf(0)
            for j in range(k + 1, n):
                M[i][j] -= m * M[k][j]
            rhs[i] -= m * rhs[k]
    x = [mp.mpf(0)] * n
    for i in range(n - 1, -1, -1):
        x[i] = (rhs[i] - sum(M[i][j] * x[j] for j in range(i + 1, n))) / M[i][i]
    return x


def errores_relativos(x, ref, dps: int = DPS) -> Tuple[float, float]:
    mp.mp.dps = dps
    xm = [v if isinstance(v, mp.mpf) else mp.mpf(float(v)) for v in x]
    d = [abs(xm[i] - ref[i]) for i in range(len(ref))]
    r = [abs(v) for v in ref]
    return float(max(d) / max(r)), float(mp.sqrt(sum(v * v for v in d)) / mp.sqrt(sum(v * v for v in r)))


def digitos_correctos(x, ref, dps: int = DPS) -> List[int]:
    mp.mp.dps = dps
    out = []
    for xi, ri in zip(x, ref):
        e = abs(mp.mpf(float(xi)) - ri) / abs(ri)
        d = MAX_DIGITOS if e == 0 else int(mp.floor(-mp.log10(2 * e)))
        out.append(max(0, min(MAX_DIGITOS, d)))
    return out


def verificar_solucion(ruta_matriz_a, ruta_vector_b, resultado: ResultadoEliminacion,
                       pasadas_refinamiento: int = 1, dps: int = DPS) -> ResultadoVerificacion:
    if pasadas_refinamiento < 1:
        raise ValueError("Se requiere al menos una pasada de refinamiento.")
    mp.mp.dps = dps
    A_mp, b_mp = cargar_datos_mp(ruta_matriz_a, ruta_vector_b, dps)
    n = len(b_mp)
    A_f = np.array([[float(v) for v in f] for f in A_mp])
    b_f = np.array([float(v) for v in b_mp])
    ref = referencia_mp(A_mp, b_mp, dps)

    x0 = np.array(resultado.solucion_x, dtype=np.float64, copy=True)
    r0, n0 = calcular_residuo_mp(A_mp, b_mp, x0, dps)

    cnt = ContadorOperaciones(num_pasos=n - 1)
    # --- refinamiento iterado en float64 (3 pasadas registradas; se reporta `pasadas_refinamiento`) ---
    x = x0.copy()
    hist_f: List[Iteracion] = []
    dx = np.zeros(n)
    x_reportada = x0.copy()
    for p in range(1, max(3, pasadas_refinamiento) + 1):
        r, _ = calcular_residuo_mp(A_mp, b_mp, x, dps)
        dx = resolver_con_lu(resultado, np.array([float(v) for v in r]), contador=cnt if p == 1 else None)
        x = x + dx
        _, nr = calcular_residuo_mp(A_mp, b_mp, x, dps)
        hist_f.append(Iteracion(p, nr, errores_relativos(x, ref, dps)[0]))
        if p == pasadas_refinamiento:
            x_reportada = x.copy()
            dx_reportada = dx.copy()
    flops_pasada = cnt.totales()["total_flops"]

    # --- refinamiento con iterado en precisión extendida (diagnóstico) ---
    xe = [mp.mpf(float(v)) for v in x0]
    hist_e: List[Iteracion] = []
    for p in range(1, 5):
        r, _ = calcular_residuo_mp(A_mp, b_mp, xe, dps)
        dxe = resolver_con_lu(resultado, np.array([float(v) for v in r]))
        xe = [xe[i] + mp.mpf(float(dxe[i])) for i in range(n)]
        _, nr = calcular_residuo_mp(A_mp, b_mp, xe, dps)
        hist_e.append(Iteracion(p, nr, errores_relativos(xe, ref, dps)[0]))

    res_ref, n_ref = calcular_residuo_mp(A_mp, b_mp, x_reportada, dps)
    x_np = np.linalg.solve(A_f, b_f)  # SOLO control

    return ResultadoVerificacion(
        residuo_inicial=r0, norma_inf_inicial=n0, correccion_dx=dx_reportada,
        solucion_refinada=x_reportada, residuo_refinado=res_ref, norma_inf_refinada=n_ref,
        solucion_numpy=x_np, referencia_alta_precision=ref,
        digitos_correctos_inicial=digitos_correctos(x0, ref, dps),
        digitos_correctos_refinada=digitos_correctos(x_reportada, ref, dps),
        digitos_correctos_numpy=digitos_correctos(x_np, ref, dps),
        error_relativo_inf_inicial=errores_relativos(x0, ref, dps)[0],
        error_relativo_inf_refinado=errores_relativos(x_reportada, ref, dps)[0],
        error_relativo_2_inicial=errores_relativos(x0, ref, dps)[1],
        error_relativo_2_refinado=errores_relativos(x_reportada, ref, dps)[1],
        error_relativo_inf_numpy=errores_relativos(x_np, ref, dps)[0],
        historial_float=hist_f, historial_extendido=hist_e,
        flops_por_pasada=flops_pasada, pasadas_refinamiento=pasadas_refinamiento,
    )
