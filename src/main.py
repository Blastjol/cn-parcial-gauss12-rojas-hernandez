from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import List, Tuple

RAIZ_PROYECTO = Path(__file__).resolve().parent.parent
if str(RAIZ_PROYECTO) not in sys.path:
    sys.path.insert(0, str(RAIZ_PROYECTO))

import mpmath as mp
import numpy as np

from src.conteo import ContadorOperaciones, escribir_reporte_conteo
from src.gauss import ResultadoEliminacion, eliminacion_gaussiana
from src.verificacion import ResultadoVerificacion, verificar_solucion


def _cargar_datos_flotantes(
    ruta_matriz_a: Path,
    ruta_vector_b: Path,
) -> Tuple[np.ndarray, np.ndarray]:
    with open(ruta_matriz_a, mode="r", encoding="utf-8", newline="") as f_a:
        filas_a = [
            [float(celda.strip()) for celda in fila if celda.strip() != ""]
            for fila in csv.reader(f_a)
            if any(celda.strip() != "" for celda in fila)
        ]

    with open(ruta_vector_b, mode="r", encoding="utf-8", newline="") as f_b:
        valores_b = [
            float(fila[0].strip())
            for fila in csv.reader(f_b)
            if fila and fila[0].strip() != ""
        ]

    A = np.array(filas_a, dtype=np.float64)
    b = np.array(valores_b, dtype=np.float64)

    if A.shape != (12, 12):
        raise ValueError(f"Se esperaba matriz A de (12, 12), se obtuvo {A.shape}.")
    if b.shape != (12,):
        raise ValueError(f"Se esperaba vector b de (12,), se obtuvo {b.shape}.")
    if np.any(A == 0.0):
        raise ValueError("La matriz A no debe contener ceros según la especificación.")

    return A, b


def _estimar_condicion_kappa2(A: np.ndarray) -> float:
    valores_singulares = np.linalg.svd(A, compute_uv=False)
    sigma_max = float(valores_singulares[0])
    sigma_min = float(valores_singulares[-1])
    if sigma_min == 0.0:
        return float("inf")
    return sigma_max / sigma_min


def _formatear_digitos_defendibles(valor_mp: mp.mpf, digitos_correctos: int) -> str:
    cifras = max(1, min(15, digitos_correctos))
    return mp.nstr(valor_mp, n=cifras)


def _escribir_solucion_csv(
    ruta_salida: Path,
    res_sin_piv: ResultadoEliminacion,
    res_con_piv: ResultadoEliminacion,
    verif_con_piv: ResultadoVerificacion,
) -> None:
    ruta_salida.parent.mkdir(parents=True, exist_ok=True)
    lineas: List[str] = [
        "incognita,x_con_pivoteo_refinada,digitos_defendibles_refinada,x_con_pivoteo_inicial,digitos_defendibles_inicial,x_sin_pivoteo,x_numpy_control,referencia_alta_precision"
    ]
    for i in range(12):
        ref_i = verif_con_piv.referencia_alta_precision[i]
        x_ref_mp = mp.mpf(str(float(verif_con_piv.solucion_refinada[i])))
        x_ini_mp = mp.mpf(str(float(res_con_piv.solucion_x[i])))
        x_sin_mp = mp.mpf(str(float(res_sin_piv.solucion_x[i])))
        x_np_mp = mp.mpf(str(float(verif_con_piv.solucion_numpy[i])))

        d_ref = verif_con_piv.digitos_correctos_refinada[i]
        d_ini = verif_con_piv.digitos_correctos_inicial[i]

        rep_ref = _formatear_digitos_defendibles(x_ref_mp, d_ref)
        rep_ini = _formatear_digitos_defendibles(x_ini_mp, d_ini)

        lineas.append(
            f"x{i + 1},"
            f"{rep_ref},"
            f"{d_ref},"
            f"{rep_ini},"
            f"{d_ini},"
            f"{mp.nstr(x_sin_mp, n=12)},"
            f"{mp.nstr(x_np_mp, n=12)},"
            f"{mp.nstr(ref_i, n=20)}"
        )
    ruta_salida.write_text("\n".join(lineas) + "\n", encoding="utf-8")


def _escribir_traza_eliminacion(
    ruta_salida: Path,
    res_sin_piv: ResultadoEliminacion,
    res_con_piv: ResultadoEliminacion,
) -> None:
    ruta_salida.parent.mkdir(parents=True, exist_ok=True)
    contenido = (
        "TRAZA COMPLETA DE ELIMINACIÓN GAUSSIANA (11 PASOS PARA n = 12)\n\n"
        f"{res_sin_piv.traza}\n"
        f"{res_con_piv.traza}\n"
    )
    ruta_salida.write_text(contenido, encoding="utf-8")


def _escribir_residuo_txt(
    ruta_salida: Path,
    res_sin_piv: ResultadoEliminacion,
    res_con_piv: ResultadoEliminacion,
    verif_sin_piv: ResultadoVerificacion,
    verif_con_piv: ResultadoVerificacion,
    kappa_2: float,
) -> None:
    ruta_salida.parent.mkdir(parents=True, exist_ok=True)
    eps_mach = np.finfo(np.float64).eps
    cota_teorica = kappa_2 * eps_mach

    lineas: List[str] = [
        "VERIFICACIÓN INDEPENDIENTE Y REFINAMIENTO ITERATIVO (PRECISIÓN EXTENDIDA mpmath, dps=80)",
        "=" * 84,
        f"Epsilon de máquina IEEE 754 (float64): {eps_mach:.6e}",
        f"Número de condición estimado kappa_2(A): {kappa_2:.6e}",
        f"Cota teórica de error relativo kappa_2(A) * eps: {cota_teorica:.6e}",
        "",
        "--- 1. ESTRATEGIA CON PIVOTEO PARCIAL ESCALADO ---",
        f"Multiplicador máximo: {res_con_piv.multiplicador_maximo:.6e}",
        f"Factor de crecimiento: {res_con_piv.factor_crecimiento:.6e}",
        f"Norma infinito residuo inicial ||r_0||_inf: {mp.nstr(verif_con_piv.norma_inf_inicial, 12)}",
        "Vector residuo inicial r_0 = b - A x_0:",
    ]
    for i, r_val in enumerate(verif_con_piv.residuo_inicial):
        lineas.append(f"  r_0[{i + 1}] = {mp.nstr(r_val, 15)}")

    lineas.extend(
        [
            f"Pasadas de refinamiento iterativo ejecutadas: {verif_con_piv.pasadas_refinamiento}",
            "Corrección obtenida dx:",
        ]
    )
    for i, dx_val in enumerate(verif_con_piv.correccion_dx):
        lineas.append(f"  dx[{i + 1}] = {float(dx_val):.6e}")

    lineas.append("Solución corregida x_refinada vs inicial x_0:")
    for i in range(12):
        lineas.append(
            f"  x[{i + 1}]: inicial={float(res_con_piv.solucion_x[i]):.12e} | "
            f"refinada={float(verif_con_piv.solucion_refinada[i]):.12e} | "
            f"dígitos ini={verif_con_piv.digitos_correctos_inicial[i]} | "
            f"dígitos ref={verif_con_piv.digitos_correctos_refinada[i]}"
        )

    lineas.extend(
        [
            f"Norma infinito residuo refinado ||r_ref||_inf: {mp.nstr(verif_con_piv.norma_inf_refinada, 12)}",
            "Vector residuo refinado r_ref = b - A x_refinada:",
        ]
    )
    for i, r_val in enumerate(verif_con_piv.residuo_refinado):
        lineas.append(f"  r_ref[{i + 1}] = {mp.nstr(r_val, 15)}")

    lineas.append("Control con numpy.linalg.solve (diferencia absoluta |x_ref - x_numpy|):")
    for i in range(12):
        diff_np = abs(float(verif_con_piv.solucion_refinada[i]) - float(verif_con_piv.solucion_numpy[i]))
        lineas.append(f"  |x[{i + 1}]_ref - x[{i + 1}]_np| = {diff_np:.6e}")

    lineas.extend(
        [
            "",
            "--- 2. ESTRATEGIA SIN PIVOTEO ---",
            f"Multiplicador máximo: {res_sin_piv.multiplicador_maximo:.6e}",
            f"Factor de crecimiento: {res_sin_piv.factor_crecimiento:.6e}",
            f"Norma infinito residuo inicial ||r_0||_inf: {mp.nstr(verif_sin_piv.norma_inf_inicial, 12)}",
            f"Norma infinito residuo tras refinamiento ||r_ref||_inf: {mp.nstr(verif_sin_piv.norma_inf_refinada, 12)}",
            f"Dígitos defendibles iniciales por componente: {verif_sin_piv.digitos_correctos_inicial}",
            f"Dígitos defendibles refinados por componente: {verif_sin_piv.digitos_correctos_refinada}",
        ]
    )

    ruta_salida.write_text("\n".join(lineas) + "\n", encoding="utf-8")


def main() -> None:
    raiz = Path(__file__).resolve().parent.parent
    ruta_a = raiz / "datos" / "matriz_a.csv"
    ruta_b = raiz / "datos" / "vector_b.csv"
    carpeta_resultados = raiz / "resultados"

    A, b = _cargar_datos_flotantes(ruta_a, ruta_b)
    kappa_2 = _estimar_condicion_kappa2(A)

    contador_sin_piv = ContadorOperaciones(num_pasos=11)
    res_sin_piv = eliminacion_gaussiana(
        A, b, usar_pivoteo_escalado=False, contador=contador_sin_piv
    )

    contador_con_piv = ContadorOperaciones(num_pasos=11)
    res_con_piv = eliminacion_gaussiana(
        A, b, usar_pivoteo_escalado=True, contador=contador_con_piv
    )

    verif_con_piv = verificar_solucion(
        ruta_matriz_a=ruta_a,
        ruta_vector_b=ruta_b,
        x_inicial=res_con_piv.solucion_x,
        usar_pivoteo_escalado=True,
        pasadas_refinamiento=1,
        dps=80,
    )

    verif_sin_piv = verificar_solucion(
        ruta_matriz_a=ruta_a,
        ruta_vector_b=ruta_b,
        x_inicial=res_sin_piv.solucion_x,
        usar_pivoteo_escalado=False,
        pasadas_refinamiento=1,
        dps=80,
    )

    escribir_reporte_conteo(
        ruta_salida=carpeta_resultados / "conteo_flops.md",
        reportes_por_estrategia={
            "Sin pivoteo": contador_sin_piv,
            "Con pivoteo parcial escalado": contador_con_piv,
        },
        n=12,
    )
    _escribir_traza_eliminacion(
        ruta_salida=carpeta_resultados / "traza_eliminacion.txt",
        res_sin_piv=res_sin_piv,
        res_con_piv=res_con_piv,
    )
    _escribir_solucion_csv(
        ruta_salida=carpeta_resultados / "solucion.csv",
        res_sin_piv=res_sin_piv,
        res_con_piv=res_con_piv,
        verif_con_piv=verif_con_piv,
    )
    _escribir_residuo_txt(
        ruta_salida=carpeta_resultados / "residuo.txt",
        res_sin_piv=res_sin_piv,
        res_con_piv=res_con_piv,
        verif_sin_piv=verif_sin_piv,
        verif_con_piv=verif_con_piv,
        kappa_2=kappa_2,
    )

    tot_sin = contador_sin_piv.totales()
    tot_con = contador_con_piv.totales()
    print("=== EJECUCIÓN COMPLETA (python src/main.py) ===")
    print(f"kappa_2(A) estimado: {kappa_2:.6e}")
    print(
        f"Sin pivoteo -> FLOPs: {tot_sin['total_flops']} | mult_max: {res_sin_piv.multiplicador_maximo:.4e} | "
        f"crec: {res_sin_piv.factor_crecimiento:.4e} | ||r_0||_inf: {mp.nstr(verif_sin_piv.norma_inf_inicial, 6)}"
    )
    print(
        f"Pivoteo escalado -> FLOPs: {tot_con['total_flops']} | sobrecosto pivoteo: "
        f"{tot_con['pivoteo_div']} div_escala, {tot_con['comparaciones']} comp, {tot_con['intercambios']} inter | "
        f"mult_max: {res_con_piv.multiplicador_maximo:.4e} | crec: {res_con_piv.factor_crecimiento:.4e} | "
        f"||r_0||_inf: {mp.nstr(verif_con_piv.norma_inf_inicial, 6)} | ||r_ref||_inf: {mp.nstr(verif_con_piv.norma_inf_refinada, 6)}"
    )
    print("\n--- TABLA 10 DIGITOS PARA readme.md ---")
    for i in range(12):
        val_10 = mp.nstr(mp.mpf(str(float(verif_con_piv.solucion_refinada[i]))), n=10)
        print(f"| $x_{{{i + 1}}}$ | `{val_10}` |")

    print("\n--- TABLA DETALLADA PARA informe.md ---")
    for i in range(12):
        d_ini = verif_con_piv.digitos_correctos_inicial[i]
        d_ref = verif_con_piv.digitos_correctos_refinada[i]
        val_ini = _formatear_digitos_defendibles(mp.mpf(str(float(res_con_piv.solucion_x[i]))), d_ini)
        val_ref = _formatear_digitos_defendibles(mp.mpf(str(float(verif_con_piv.solucion_refinada[i]))), d_ref)
        print(f"| $x_{{{i + 1}}}$ | `{val_ini}` | {d_ini} | `{val_ref}` | {d_ref} |")


if __name__ == "__main__":
    main()