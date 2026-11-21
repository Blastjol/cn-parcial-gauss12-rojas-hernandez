from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List


@dataclass
class RegistroPaso:
    paso: int
    mult_div: int = 0
    sum_resta: int = 0
    pivoteo_div: int = 0
    comparaciones: int = 0
    intercambios: int = 0

    @property
    def flops_aritmeticos(self) -> int:
        return self.mult_div + self.sum_resta


@dataclass
class ContadorOperaciones:
    """Contador instrumentado de operaciones aritméticas y sobrecosto de pivoteo."""

    num_pasos: int = 11
    pasos: Dict[int, RegistroPaso] = field(init=False)
    extra_mult_div: int = 0
    extra_sum_resta: int = 0
    extra_pivoteo_div: int = 0
    extra_comparaciones: int = 0
    extra_intercambios: int = 0

    def __post_init__(self) -> None:
        self.pasos = {k: RegistroPaso(paso=k) for k in range(1, self.num_pasos + 1)}

    def _obtener_registro(self, paso: int | None) -> RegistroPaso | None:
        if paso is None:
            return None
        if paso not in self.pasos:
            raise ValueError(f"Paso fuera de rango (1..{self.num_pasos}): {paso}")
        return self.pasos[paso]

    def registrar_operacion(
        self, categoria: str, paso: int | None = None, cantidad: int = 1
    ) -> None:
        if cantidad < 0:
            raise ValueError("La cantidad de operaciones no puede ser negativa.")
        reg = self._obtener_registro(paso)
        if categoria == "mult_div":
            if reg is not None:
                reg.mult_div += cantidad
            else:
                self.extra_mult_div += cantidad
        elif categoria == "sum_resta":
            if reg is not None:
                reg.sum_resta += cantidad
            else:
                self.extra_sum_resta += cantidad
        else:
            raise ValueError(f"Categoría de operación no válida: {categoria}")

    def registrar_pivoteo_div(self, paso: int | None = None, cantidad: int = 1) -> None:
        if cantidad < 0:
            raise ValueError("La cantidad de divisiones de pivoteo no puede ser negativa.")
        reg = self._obtener_registro(paso)
        if reg is not None:
            reg.pivoteo_div += cantidad
        else:
            self.extra_pivoteo_div += cantidad

    def registrar_comparacion(self, paso: int | None = None, cantidad: int = 1) -> None:
        if cantidad < 0:
            raise ValueError("La cantidad de comparaciones no puede ser negativa.")
        reg = self._obtener_registro(paso)
        if reg is not None:
            reg.comparaciones += cantidad
        else:
            self.extra_comparaciones += cantidad

    def registrar_intercambio(self, paso: int | None = None, cantidad: int = 1) -> None:
        if cantidad < 0:
            raise ValueError("La cantidad de intercambios no puede ser negativa.")
        reg = self._obtener_registro(paso)
        if reg is not None:
            reg.intercambios += cantidad
        else:
            self.extra_intercambios += cantidad

    def totales(self) -> Dict[str, int]:
        total_mult_div = sum(r.mult_div for r in self.pasos.values()) + self.extra_mult_div
        total_sum_resta = sum(r.sum_resta for r in self.pasos.values()) + self.extra_sum_resta
        total_flops = total_mult_div + total_sum_resta
        total_pivoteo_div = sum(r.pivoteo_div for r in self.pasos.values()) + self.extra_pivoteo_div
        total_comparaciones = (
            sum(r.comparaciones for r in self.pasos.values()) + self.extra_comparaciones
        )
        total_intercambios = (
            sum(r.intercambios for r in self.pasos.values()) + self.extra_intercambios
        )

        return {
            "mult_div": total_mult_div,
            "sum_resta": total_sum_resta,
            "total_flops": total_flops,
            "pivoteo_div": total_pivoteo_div,
            "comparaciones": total_comparaciones,
            "intercambios": total_intercambios,
        }

    def filas_tabla(self) -> List[Dict[str, int | str]]:
        filas: List[Dict[str, int | str]] = []
        for k in range(1, self.num_pasos + 1):
            r = self.pasos[k]
            filas.append(
                {
                    "etapa": f"Paso {k}",
                    "mult_div": r.mult_div,
                    "sum_resta": r.sum_resta,
                    "flops_paso": r.flops_aritmeticos,
                    "pivoteo_div": r.pivoteo_div,
                    "comparaciones": r.comparaciones,
                    "intercambios": r.intercambios,
                }
            )
        if (
            self.extra_mult_div > 0
            or self.extra_sum_resta > 0
            or self.extra_pivoteo_div > 0
            or self.extra_comparaciones > 0
            or self.extra_intercambios > 0
        ):
            filas.append(
                {
                    "etapa": "Sustitución regresiva",
                    "mult_div": self.extra_mult_div,
                    "sum_resta": self.extra_sum_resta,
                    "flops_paso": self.extra_mult_div + self.extra_sum_resta,
                    "pivoteo_div": self.extra_pivoteo_div,
                    "comparaciones": self.extra_comparaciones,
                    "intercambios": self.extra_intercambios,
                }
            )
        return filas


def escribir_reporte_conteo(
    ruta_salida: str | Path,
    reportes_por_estrategia: Dict[str, ContadorOperaciones],
    n: int = 12,
) -> None:
    ruta = Path(ruta_salida)
    ruta.parent.mkdir(parents=True, exist_ok=True)

    cota_teorica_2n3_3 = (2 * (n**3)) / 3.0
    lineas: List[str] = [
        "# Conteo instrumentado de operaciones (FLOPs)",
        "",
        "## Convención de conteo",
        "- **1 flop = 1 operación aritmética elemental** (`+`, `-`, `*`, `/`) en eliminación, actualización de `b` y sustitución regresiva.",
        "- Cada multiplicador $m_{ik} = a_{ik}/a_{kk}$ se calcula **una sola vez por fila** y el barrido va desde $j = k+1$ hasta $n$ (sin operar sobre los 66 elementos subdiagonales anulados).",
        "- El sobrecosto del pivoteo parcial escalado (**divisiones de escalamiento** $|a_{ik}|/s_i$, **comparaciones** de búsqueda e **intercambios** de filas) se contabiliza por separado.",
        f"- Referencia asintótica dominante $2n^3/3$ para $n={n}$: `{cota_teorica_2n3_3:.2f}` FLOPs; cota de referencia del enunciado (incluyendo `b` y sustitución regresiva): `~1486` FLOPs (rango aceptado $\\pm 10\\%$: `[1337, 1635]`).",
        "",
    ]

    for nombre, contador in reportes_por_estrategia.items():
        tot = contador.totales()
        lineas.extend(
            [
                f"## Estrategia: {nombre}",
                "",
                "| Etapa / Paso | Mult / Div | Sum / Resta | FLOPs eliminación/sust. | Div. escalamiento | Comparaciones | Intercambios |",
                "|---|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for fila in contador.filas_tabla():
            lineas.append(
                f"| {fila['etapa']} | {fila['mult_div']} | {fila['sum_resta']} | "
                f"{fila['flops_paso']} | {fila['pivoteo_div']} | {fila['comparaciones']} | {fila['intercambios']} |"
            )
        lineas.extend(
            [
                f"| **Total** | **{tot['mult_div']}** | **{tot['sum_resta']}** | "
                f"**{tot['total_flops']}** | **{tot['pivoteo_div']}** | **{tot['comparaciones']}** | **{tot['intercambios']}** |",
                "",
            ]
        )

    ruta.write_text("\n".join(lineas), encoding="utf-8")
