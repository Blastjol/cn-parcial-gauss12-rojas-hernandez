from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple
import numpy as np

from src.conteo import ContadorOperaciones


@dataclass
class ResultadoEliminacion:
    estrategia: str
    solucion_x: np.ndarray
    matriz_u: np.ndarray
    vector_b_mod: np.ndarray
    multiplicadores: List[Tuple[int, int, float]]
    intercambios_realizados: List[Tuple[int, int, int]]
    multiplicador_maximo: float
    factor_crecimiento: float
    traza: str


def _formatear_matriz(M: np.ndarray) -> str:
    filas_texto: List[str] = []
    for fila in M:
        elementos = "  ".join(f"{val:14.6e}" for val in fila)
        filas_texto.append(f"  [ {elementos} ]")
    return "\n".join(filas_texto)


def _calcular_factores_escala(
    A: np.ndarray,
    fila_inicio: int,
    contador: ContadorOperaciones,
    paso: int,
) -> np.ndarray:
    n = A.shape[0]
    escalas = np.zeros(n, dtype=np.float64)
    for i in range(fila_inicio, n):
        max_val = abs(float(A[i, fila_inicio]))
        for j in range(fila_inicio + 1, n):
            contador.registrar_comparacion(paso=paso, cantidad=1)
            candidato = abs(float(A[i, j]))
            if candidato > max_val:
                max_val = candidato
        escalas[i] = max_val
    return escalas


def _elegir_pivote_escalado(
    A: np.ndarray,
    escalas: np.ndarray,
    k: int,
    contador: ContadorOperaciones,
    paso: int,
) -> int:
    n = A.shape[0]
    mejor_fila = k
    escala_k = float(escalas[k])
    if escala_k == 0.0:
        raise ZeroDivisionError(f"Escala nula en fila {k + 1} durante el paso {paso}.")
    mejor_razon = abs(float(A[k, k])) / escala_k
    contador.registrar_pivoteo_div(paso=paso, cantidad=1)

    for i in range(k + 1, n):
        escala_i = float(escalas[i])
        if escala_i == 0.0:
            raise ZeroDivisionError(f"Escala nula en fila {i + 1} durante el paso {paso}.")
        razon_i = abs(float(A[i, k])) / escala_i
        contador.registrar_pivoteo_div(paso=paso, cantidad=1)
        contador.registrar_comparacion(paso=paso, cantidad=1)
        if razon_i > mejor_razon:
            mejor_razon = razon_i
            mejor_fila = i

    return mejor_fila


def _intercambiar_filas(
    A: np.ndarray,
    b: np.ndarray,
    escalas: Optional[np.ndarray],
    fila_a: int,
    fila_b: int,
    contador: ContadorOperaciones,
    paso: int,
) -> None:
    if fila_a == fila_b:
        return
    A[[fila_a, fila_b], :] = A[[fila_b, fila_a], :]
    b[[fila_a, fila_b]] = b[[fila_b, fila_a]]
    if escalas is not None:
        escalas[[fila_a, fila_b]] = escalas[[fila_b, fila_a]]
    contador.registrar_intercambio(paso=paso, cantidad=1)


def _sustitucion_regresiva(
    U: np.ndarray,
    b_mod: np.ndarray,
    contador: ContadorOperaciones,
) -> np.ndarray:
    n = U.shape[0]
    x = np.zeros(n, dtype=np.float64)

    for i in range(n - 1, -1, -1):
        suma = float(b_mod[i])
        for j in range(i + 1, n):
            producto = float(U[i, j]) * float(x[j])
            contador.registrar_operacion("mult_div", paso=None, cantidad=1)
            suma = suma - producto
            contador.registrar_operacion("sum_resta", paso=None, cantidad=1)

        pivote = float(U[i, i])
        if pivote == 0.0:
            raise ZeroDivisionError(f"Pivote nulo en sustitución regresiva en fila {i + 1}.")
        x[i] = suma / pivote
        contador.registrar_operacion("mult_div", paso=None, cantidad=1)

    return x


def eliminacion_gaussiana(
    A_inicial: np.ndarray,
    b_inicial: np.ndarray,
    usar_pivoteo_escalado: bool = False,
    contador: Optional[ContadorOperaciones] = None,
) -> ResultadoEliminacion:
    """Eliminación gaussiana densa sin pivoteo o con pivoteo parcial escalado."""
    A = np.array(A_inicial, dtype=np.float64, copy=True)
    b = np.array(b_inicial, dtype=np.float64, copy=True)

    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("La matriz A debe ser cuadrada.")
    n = A.shape[0]
    if b.shape != (n,):
        raise ValueError(f"El vector b debe tener forma ({n},).")

    if contador is None:
        contador = ContadorOperaciones(num_pasos=n - 1)

    estrategia = (
        "Con pivoteo parcial escalado" if usar_pivoteo_escalado else "Sin pivoteo"
    )

    max_abs_inicial = float(np.max(np.abs(A)))
    max_abs_durante = max_abs_inicial
    max_mult_abs = 0.0

    multiplicadores: List[Tuple[int, int, float]] = []
    intercambios_realizados: List[Tuple[int, int, int]] = []
    bloques_traza: List[str] = [
        "============================================================",
        f"ESTRATEGIA: {estrategia.upper()}",
        "============================================================",
        "Matriz inicial A^(0):",
        _formatear_matriz(A),
        "",
    ]

    for k in range(n - 1):
        paso = k + 1
        hubo_intercambio = False
        fila_intercambiada_con = k

        if usar_pivoteo_escalado:
            escalas = _calcular_factores_escala(A, k, contador, paso)
            fila_pivote = _elegir_pivote_escalado(A, escalas, k, contador, paso)
            if fila_pivote != k:
                _intercambiar_filas(A, b, escalas, k, fila_pivote, contador, paso)
                hubo_intercambio = True
                fila_intercambiada_con = fila_pivote
                intercambios_realizados.append((paso, k + 1, fila_pivote + 1))

        pivote = float(A[k, k])
        if pivote == 0.0:
            raise ZeroDivisionError(f"Pivote nulo encontrado en el paso {paso} (fila {k + 1}).")

        mults_paso: List[str] = []
        for i in range(k + 1, n):
            m_ik = float(A[i, k]) / pivote
            contador.registrar_operacion("mult_div", paso=paso, cantidad=1)

            abs_m = abs(m_ik)
            if abs_m > max_mult_abs:
                max_mult_abs = abs_m

            multiplicadores.append((paso, i + 1, m_ik))
            mults_paso.append(f"m({i + 1},{k + 1})={m_ik:.6e}")

            A[i, k] = 0.0

            for j in range(k + 1, n):
                producto = m_ik * float(A[k, j])
                contador.registrar_operacion("mult_div", paso=paso, cantidad=1)
                A[i, j] = float(A[i, j]) - producto
                contador.registrar_operacion("sum_resta", paso=paso, cantidad=1)

                abs_val = abs(float(A[i, j]))
                if abs_val > max_abs_durante:
                    max_abs_durante = abs_val

            prod_b = m_ik * float(b[k])
            contador.registrar_operacion("mult_div", paso=paso, cantidad=1)
            b[i] = float(b[i]) - prod_b
            contador.registrar_operacion("sum_resta", paso=paso, cantidad=1)

        bloque = [
            f"--- Paso {paso} (columna pivote c{k + 1}, pivote={pivote:.6e}) ---",
            f"Intercambio de filas: {'Sí (fila ' + str(k + 1) + ' <-> fila ' + str(fila_intercambiada_con + 1) + ')' if hubo_intercambio else 'No'}",
            f"Multiplicadores fila {k + 1}: {', '.join(mults_paso)}",
            f"Matriz tras paso {paso}:",
            _formatear_matriz(A),
            "",
        ]
        bloques_traza.extend(bloque)

    factor_crecimiento = (
        max_abs_durante / max_abs_inicial if max_abs_inicial != 0.0 else float("inf")
    )
    solucion_x = _sustitucion_regresiva(A, b, contador)

    return ResultadoEliminacion(
        estrategia=estrategia,
        solucion_x=solucion_x,
        matriz_u=A,
        vector_b_mod=b,
        multiplicadores=multiplicadores,
        intercambios_realizados=intercambios_realizados,
        multiplicador_maximo=max_mult_abs,
        factor_crecimiento=factor_crecimiento,
        traza="\n".join(bloques_traza),
    )