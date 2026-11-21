from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np

from src.conteo import ContadorOperaciones

PIVOTEOS = ("ninguno", "parcial", "escalado")
NOMBRES = {"ninguno": "Sin pivoteo",
           "parcial": "Con pivoteo parcial",
           "escalado": "Con pivoteo parcial escalado"}


@dataclass
class ResultadoEliminacion:
    estrategia: str
    pivoteo: str
    solucion_x: np.ndarray
    matriz_l: np.ndarray
    matriz_u: np.ndarray
    permutacion: List[int]
    vector_b_mod: np.ndarray
    multiplicadores: List[Tuple[int, int, float]]
    intercambios_realizados: List[Tuple[int, int, int]]
    multiplicador_maximo: float
    factor_crecimiento: float
    bits_absorbidos_max: float
    bits_absorbidos_por_paso: List[float]
    no_nulos_bajo_diagonal: int
    traza: str


def _fmt_matriz(M: np.ndarray) -> str:
    return "\n".join("  [ " + " ".join(f"{v:24.16e}" for v in fila) + " ]" for fila in M)


def _fmt_vector(v: np.ndarray) -> str:
    return "  [ " + " ".join(f"{x:24.16e}" for x in v) + " ]"


def _escalas_iniciales(A: np.ndarray, contador: ContadorOperaciones) -> np.ndarray:
    """s_i = max_j |a_ij^(0)|, una sola vez: n(n-1) comparaciones (paso 0 = previo)."""
    n = A.shape[0]
    s = np.zeros(n)
    for i in range(n):
        m = abs(float(A[i, 0]))
        for j in range(1, n):
            contador.registrar_comparacion(paso=0)
            if abs(float(A[i, j])) > m:
                m = abs(float(A[i, j]))
        s[i] = m
    return s


def _elegir_pivote(A: np.ndarray, k: int, pivoteo: str, escalas: Optional[np.ndarray],
                   contador: ContadorOperaciones, paso: int) -> int:
    n = A.shape[0]
    if pivoteo == "parcial":
        mejor, mejor_val = k, abs(float(A[k, k]))
        for i in range(k + 1, n):
            contador.registrar_comparacion(paso=paso)
            if abs(float(A[i, k])) > mejor_val:
                mejor, mejor_val = i, abs(float(A[i, k]))
        return mejor
    # escalado
    assert escalas is not None
    if escalas[k] == 0.0:
        raise ZeroDivisionError(f"Escala nula en fila {k + 1}.")
    mejor, mejor_val = k, abs(float(A[k, k])) / float(escalas[k])
    contador.registrar_pivoteo_div(paso=paso)
    for i in range(k + 1, n):
        if escalas[i] == 0.0:
            raise ZeroDivisionError(f"Escala nula en fila {i + 1}.")
        razon = abs(float(A[i, k])) / float(escalas[i])
        contador.registrar_pivoteo_div(paso=paso)
        contador.registrar_comparacion(paso=paso)
        if razon > mejor_val:
            mejor, mejor_val = i, razon
    return mejor


def _intercambiar(A, L, b, perm, escalas, a, c, contador, paso) -> None:
    A[[a, c], :] = A[[c, a], :]
    if a > 0:
        L[[a, c], :a] = L[[c, a], :a]
    b[[a, c]] = b[[c, a]]
    perm[a], perm[c] = perm[c], perm[a]
    if escalas is not None:
        escalas[[a, c]] = escalas[[c, a]]
    contador.registrar_intercambio(paso=paso)


def sustitucion_regresiva(U: np.ndarray, y: np.ndarray,
                          contador: Optional[ContadorOperaciones] = None) -> np.ndarray:
    n = U.shape[0]
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        suma = float(y[i])
        for j in range(i + 1, n):
            suma = suma - float(U[i, j]) * float(x[j])
            if contador is not None:
                contador.registrar_operacion("mult_div", paso=None)
                contador.registrar_operacion("sum_resta", paso=None)
        if U[i, i] == 0.0:
            raise ZeroDivisionError(f"Pivote nulo en la fila {i + 1}.")
        x[i] = suma / float(U[i, i])
        if contador is not None:
            contador.registrar_operacion("mult_div", paso=None)
    return x


def resolver_con_lu(res: ResultadoEliminacion, nuevo_b: np.ndarray,
                    contador: Optional[ContadorOperaciones] = None) -> np.ndarray:
    """Resuelve A z = nuevo_b reutilizando P, L, U (sin refactorizar): L y = P b, U z = y."""
    L, U, perm = res.matriz_l, res.matriz_u, res.permutacion
    n = U.shape[0]
    bp = np.array([float(nuevo_b[i]) for i in perm])
    y = np.zeros(n)
    for i in range(n):
        s = float(bp[i])
        for j in range(i):
            s -= float(L[i, j]) * float(y[j])
            if contador is not None:
                contador.registrar_operacion("mult_div", paso=None)
                contador.registrar_operacion("sum_resta", paso=None)
        y[i] = s
    return sustitucion_regresiva(U, y, contador)


def eliminacion_gaussiana(A_inicial: np.ndarray, b_inicial: np.ndarray,
                          pivoteo: str = "escalado",
                          contador: Optional[ContadorOperaciones] = None) -> ResultadoEliminacion:
    if pivoteo not in PIVOTEOS:
        raise ValueError(f"pivoteo debe ser uno de {PIVOTEOS}")
    A = np.array(A_inicial, dtype=np.float64, copy=True)
    b = np.array(b_inicial, dtype=np.float64, copy=True)
    if A.ndim != 2 or A.shape[0] != A.shape[1] or b.shape != (A.shape[0],):
        raise ValueError("Dimensiones incompatibles.")
    n = A.shape[0]
    contador = contador or ContadorOperaciones(num_pasos=n - 1)

    L = np.eye(n)
    perm = list(range(n))
    max_ini = float(np.max(np.abs(A)))
    max_dur = max_ini
    max_mult = 0.0
    bits_max = 0.0
    bits_pasos: List[float] = []
    escalas = _escalas_iniciales(A, contador) if pivoteo == "escalado" else None

    multiplicadores: List[Tuple[int, int, float]] = []
    intercambios: List[Tuple[int, int, int]] = []
    traza: List[str] = [
        "=" * 100, f"ESTRATEGIA: {NOMBRES[pivoteo].upper()}", "=" * 100,
        "Matriz inicial A^(0) (17 cifras significativas):", _fmt_matriz(A),
        "Vector b^(0):", _fmt_vector(b), "",
    ]

    for k in range(n - 1):
        paso = k + 1
        swap = k
        if pivoteo != "ninguno":
            p = _elegir_pivote(A, k, pivoteo, escalas, contador, paso)
            if p != k:
                _intercambiar(A, L, b, perm, escalas, k, p, contador, paso)
                swap = p
                intercambios.append((paso, k + 1, p + 1))
        pivote = float(A[k, k])
        if pivote == 0.0:
            raise ZeroDivisionError(f"Pivote nulo en el paso {paso}.")

        mults: List[str] = []
        bits_paso = 0.0
        for i in range(k + 1, n):
            m = float(A[i, k]) / pivote                      # multiplicador: UNA división por fila
            contador.registrar_operacion("mult_div", paso=paso)
            max_mult = max(max_mult, abs(m))
            L[i, k] = m
            multiplicadores.append((paso, i + 1, m))
            mults.append(f"m({i + 1},{k + 1})={m:.6e}")
            A[i, k] = 0.0                                    # anulado sin operar
            for j in range(k + 1, n):                        # barrido desde la columna k+1
                prod = m * float(A[k, j])
                contador.registrar_operacion("mult_div", paso=paso)
                # diagnóstico (no contado): bits del operando original absorbidos por |m*a_kj|
                if A[i, j] != 0.0 and abs(prod) > abs(float(A[i, j])):
                    bits_paso = max(bits_paso, math.log2(abs(prod) / abs(float(A[i, j]))))
                A[i, j] = float(A[i, j]) - prod
                contador.registrar_operacion("sum_resta", paso=paso)
                max_dur = max(max_dur, abs(float(A[i, j])))
            b[i] = float(b[i]) - m * float(b[k])
            contador.registrar_operacion("mult_div", paso=paso)
            contador.registrar_operacion("sum_resta", paso=paso)
        bits_pasos.append(bits_paso)
        bits_max = max(bits_max, bits_paso)

        traza += [
            f"--- Paso {paso}: columna pivote c{paso}, pivote = {pivote:.16e} ---",
            f"Intercambio de filas: " + (f"sí (fila {k + 1} <-> fila {swap + 1})" if swap != k else "no"),
            f"Multiplicadores: {', '.join(mults)}",
            f"Bits absorbidos máx. en este paso: {bits_paso:.2f}",
            f"Matriz tras el paso {paso}:", _fmt_matriz(A),
            f"Vector b tras el paso {paso}:", _fmt_vector(b), "",
        ]

    no_nulos = int(np.count_nonzero(np.tril(A, -1)))
    traza += [f"Elementos no nulos bajo la diagonal: {no_nulos} de {n * (n - 1) // 2}",
              f"Permutación de filas (índices originales, base 1): {[p + 1 for p in perm]}", ""]

    x = sustitucion_regresiva(A, b, contador)
    return ResultadoEliminacion(
        estrategia=NOMBRES[pivoteo], pivoteo=pivoteo, solucion_x=x, matriz_l=L, matriz_u=A,
        permutacion=perm, vector_b_mod=b, multiplicadores=multiplicadores,
        intercambios_realizados=intercambios, multiplicador_maximo=max_mult,
        factor_crecimiento=max_dur / max_ini, bits_absorbidos_max=bits_max,
        bits_absorbidos_por_paso=bits_pasos, no_nulos_bajo_diagonal=no_nulos,
        traza="\n".join(traza),
    )
