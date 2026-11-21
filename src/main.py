from __future__ import annotations

import csv
import math
import sys
from pathlib import Path
from typing import Dict, List, Tuple

RAIZ = Path(__file__).resolve().parent.parent
if str(RAIZ) not in sys.path:
    sys.path.insert(0, str(RAIZ))

import mpmath as mp  # noqa: E402
import numpy as np  # noqa: E402

from src.conteo import ContadorOperaciones, escribir_reporte_conteo  # noqa: E402
from src.gauss import NOMBRES, ResultadoEliminacion, eliminacion_gaussiana  # noqa: E402
from src.verificacion import DPS, ResultadoVerificacion, digitos_correctos, verificar_solucion  # noqa: E402

N = 12


def cargar_datos(ruta_a: Path, ruta_b: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Lee A y b SIN truncar (float(str) redondea correctamente al double más cercano)."""
    with open(ruta_a, encoding="utf-8", newline="") as f:
        A = [[float(c.strip()) for c in fila if c.strip()] for fila in csv.reader(f) if any(c.strip() for c in fila)]
    with open(ruta_b, encoding="utf-8", newline="") as f:
        b = [float(fila[0].strip()) for fila in csv.reader(f) if fila and fila[0].strip()]
    A, b = np.array(A), np.array(b)
    if A.shape != (N, N) or b.shape != (N,):
        raise ValueError(f"Dimensiones inesperadas: A{A.shape}, b{b.shape}")
    if np.any(A == 0.0):
        raise ValueError("A debe ser densa (sin ceros).")
    return A, b


def kappa2(A: np.ndarray) -> float:
    s = np.linalg.svd(A, compute_uv=False)
    return float(s[0] / s[-1])


def diagnostico_condicionamiento(A: np.ndarray) -> Dict[str, float]:
    equil = A / np.max(np.abs(A), axis=1)[:, None]          # D^-1 A, D = diag(s_i)
    B = np.vstack([A[2], A[4]]).T                            # columnas: filas 3 y 5
    w = A[6]
    peso = 1.0 / np.abs(w)                                   # error relativo por columna
    coef = np.linalg.lstsq(B * peso[:, None], w * peso, rcond=None)[0]
    desv = float(np.max(np.abs((B @ coef - w) / w)))
    return {"kappa2": kappa2(A), "kappa1": float(np.linalg.cond(A, 1)), "kappa_inf": float(np.linalg.cond(A, np.inf)),
            "kappa2_equilibrada": float(np.linalg.cond(equil, 2)),
            "alfa_f3": float(coef[0]), "beta_f5": float(coef[1]), "desv_rel_f7": desv,
            "det": float(np.linalg.det(A))}


def fmt(x, d: int) -> str:
    """x con exactamente d cifras significativas (conserva ceros finales)."""
    return mp.nstr(mp.mpf(float(x)), d, strip_zeros=False, min_fixed=-20, max_fixed=20)


def escribir_solucion(ruta: Path, res: Dict[str, ResultadoEliminacion],
                      ver: Dict[str, ResultadoVerificacion]) -> None:
    v = ver["escalado"]
    mp.mp.dps = DPS
    ref = v.referencia_alta_precision
    d_sin = digitos_correctos(res["ninguno"].solucion_x, ref)
    with open(ruta, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["incognita", "x_reportada", "digitos_correctos", "x_float64_17cifras",
                    "x_inicial_pivoteo_escalado", "digitos_inicial", "x_sin_pivoteo", "digitos_sin_pivoteo",
                    "x_numpy_control", "referencia_mpmath_20cifras"])
        for i in range(N):
            d = v.digitos_correctos_refinada[i]
            w.writerow([f"x{i + 1}", fmt(v.solucion_refinada[i], d), d, repr(float(v.solucion_refinada[i])),
                        fmt(res["escalado"].solucion_x[i], v.digitos_correctos_inicial[i]), v.digitos_correctos_inicial[i],
                        fmt(res["ninguno"].solucion_x[i], max(1, d_sin[i])), d_sin[i],
                        fmt(v.solucion_numpy[i], max(1, v.digitos_correctos_numpy[i])),
                        mp.nstr(ref[i], 20, strip_zeros=False, min_fixed=-20, max_fixed=20)])


def escribir_traza(ruta: Path, res: Dict[str, ResultadoEliminacion]) -> None:
    ruta.write_text("TRAZA DE ELIMINACIÓN GAUSSIANA — 11 pasos para n = 12 (A y b con 17 cifras significativas)\n\n"
                    + res["escalado"].traza + "\n\n" + res["ninguno"].traza + "\n", encoding="utf-8")


def _linea_hist(h) -> str:
    return f"  pasada {h.pasada}: ||r||_inf = {mp.nstr(h.norma_inf_residuo, 6)} | error relativo inf = {h.error_rel_inf:.3e}"


def escribir_residuo(ruta: Path, res: Dict[str, ResultadoEliminacion], ver: Dict[str, ResultadoVerificacion],
                     kap: float, diag: Dict[str, float]) -> None:
    eps = float(np.finfo(np.float64).eps)
    L: List[str] = [
        f"VERIFICACIÓN INDEPENDIENTE — residuo r = b - A x en precisión extendida (mpmath, dps = {DPS})",
        "=" * 100,
        f"epsilon de máquina (float64): {eps:.6e}",
        f"kappa_2(A) = sigma_max/sigma_min (SVD, solo diagnóstico): {kap:.6e}",
        f"cota clásica kappa_2 * eps: {kap * eps:.6e}  -> {-math.log10(kap * eps):.2f} dígitos garantizados",
        f"kappa_1(A) = {diag['kappa1']:.6e} | kappa_inf(A) = {diag['kappa_inf']:.6e} | det(A) = {diag['det']:.6e}",
        f"kappa_2 tras equilibrar filas (D^-1 A, s_i = max_j|a_ij|): {diag['kappa2_equilibrada']:.6e}",
        f"Casi dependencia: F7 ≈ {diag['alfa_f3']:.5f}*F3 + {diag['beta_f5']:.5f}*F5 (ajuste de mínimos cuadrados "
        f"relativo); máx. desviación relativa por columna = {diag['desv_rel_f7']:.3e}",
        "",
    ]
    for clave in ("escalado", "parcial", "ninguno"):
        r, v = res[clave], ver[clave]
        L += [f"--- {NOMBRES[clave].upper()} ---",
              f"Multiplicador máximo: {r.multiplicador_maximo:.6e} (2^{math.log2(r.multiplicador_maximo):.2f})",
              f"Factor de crecimiento rho: {r.factor_crecimiento:.6e} ({math.log2(r.factor_crecimiento):.2f} bits)",
              f"Bits absorbidos máximos en una resta (diagnóstico): {r.bits_absorbidos_max:.2f}",
              f"Elementos no nulos bajo la diagonal de U: {r.no_nulos_bajo_diagonal} de {N * (N - 1) // 2}",
              f"Error relativo inicial: inf = {v.error_relativo_inf_inicial:.4e} | 2 = {v.error_relativo_2_inicial:.4e}",
              f"||r_0||_inf = {mp.nstr(v.norma_inf_inicial, 8)}",
              f"Dígitos correctos iniciales por incógnita: {v.digitos_correctos_inicial}",
              f"Refinamiento iterado en float64 (reutilizando P, L, U):"]
        L += [_linea_hist(h) for h in v.historial_float]
        L += ["Refinamiento iterado en precisión extendida (diagnóstico de contracción):"]
        L += [_linea_hist(h) for h in v.historial_extendido]
        L.append("")

    v = ver["escalado"]
    L += ["=" * 100, "DETALLE DE LA SOLUCIÓN REPORTADA (pivoteo parcial escalado + 1 pasada de refinamiento)", "",
          "Vector residuo inicial r_0 = b - A x_0:"]
    L += [f"  r_0[{i + 1:2d}] = {mp.nstr(val, 15)}" for i, val in enumerate(v.residuo_inicial)]
    L += [f"||r_0||_inf = {mp.nstr(v.norma_inf_inicial, 12)}", "",
          f"Corrección dx (1 pasada, costo {v.flops_por_pasada} flops):"]
    L += [f"  dx[{i + 1:2d}] = {float(val):.6e}" for i, val in enumerate(v.correccion_dx)]
    L += ["", "Vector residuo refinado r = b - A x_refinada:"]
    L += [f"  r[{i + 1:2d}] = {mp.nstr(val, 15)}" for i, val in enumerate(v.residuo_refinado)]
    L += [f"||r||_inf = {mp.nstr(v.norma_inf_refinada, 12)}", "",
          f"Error relativo tras refinar: inf = {v.error_relativo_inf_refinado:.4e} | 2 = {v.error_relativo_2_refinado:.4e}",
          f"Dígitos correctos tras refinar: {v.digitos_correctos_refinada}", "",
          "Control con numpy.linalg.solve (solo control; no sustituye la implementación propia):",
          f"  error relativo inf de numpy = {v.error_relativo_inf_numpy:.4e}; dígitos correctos: {v.digitos_correctos_numpy}"]
    L += [f"  |x_refinada[{i + 1}] - x_numpy[{i + 1}]| = {abs(float(v.solucion_refinada[i]) - float(v.solucion_numpy[i])):.3e}"
          for i in range(N)]
    L += ["", "Nota: el residuo ||r||_inf no baja tras refinar porque x está almacenado en float64; su piso es",
          "u*|A|*|x| (≈ 1e-10 aquí). Lo que mejora es el error en x (forward error), no el residuo."]
    ruta.write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> None:
    carpeta = RAIZ / "resultados"
    carpeta.mkdir(exist_ok=True)
    ruta_a, ruta_b = RAIZ / "datos" / "matriz_a.csv", RAIZ / "datos" / "vector_b.csv"
    A, b = cargar_datos(ruta_a, ruta_b)
    diag = diagnostico_condicionamiento(A)
    kap = diag["kappa2"]

    contadores = {k: ContadorOperaciones(num_pasos=N - 1) for k in ("ninguno", "parcial", "escalado")}
    res = {k: eliminacion_gaussiana(A, b, pivoteo=k, contador=contadores[k]) for k in contadores}
    ver = {k: verificar_solucion(ruta_a, ruta_b, res[k], pasadas_refinamiento=1) for k in res}

    escribir_reporte_conteo(
        carpeta / "conteo_flops.md",
        {NOMBRES["ninguno"]: contadores["ninguno"], NOMBRES["parcial"]: contadores["parcial"],
         NOMBRES["escalado"]: contadores["escalado"]},
        n=N, flops_refinamiento_por_pasada=ver["escalado"].flops_por_pasada)
    escribir_traza(carpeta / "traza_eliminacion.txt", res)
    escribir_solucion(carpeta / "solucion.csv", res, ver)
    escribir_residuo(carpeta / "residuo.txt", res, ver, kap, diag)

    print("=== python src/main.py — ejecución completa ===")
    print(f"kappa_2(A) = {kap:.6e}")
    for k in ("ninguno", "parcial", "escalado"):
        t, r, v = contadores[k].totales(), res[k], ver[k]
        print(f"{NOMBRES[k]:32s} flops={t['total_flops']} sobrecosto=({t['pivoteo_div']} div, {t['comparaciones']} comp, "
              f"{t['intercambios']} int) m_max={r.multiplicador_maximo:.4e} rho={r.factor_crecimiento:.4e} "
              f"bits_abs={r.bits_absorbidos_max:.2f} err_ini={v.error_relativo_inf_inicial:.4e} "
              f"err_ref={v.error_relativo_inf_refinado:.4e} ||r0||={mp.nstr(v.norma_inf_inicial, 5)} "
              f"dig_ini={min(v.digitos_correctos_inicial)}-{max(v.digitos_correctos_inicial)}")
    print("Archivos regenerados en resultados/: conteo_flops.md, traza_eliminacion.txt, solucion.csv, residuo.txt")


if __name__ == "__main__":
    main()
