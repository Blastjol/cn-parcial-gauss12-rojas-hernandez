from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List


@dataclass
class RegistroPaso:
    paso: int
    mult_div: int = 0
    sum_resta: int = 0
    comparaciones: int = 0
    intercambios: int = 0

    @property
    def flops_aritmeticos(self) -> int:
        return self.mult_div + self.sum_resta


@dataclass
class ContadorOperaciones:
    """Contador instrumentado de operaciones elementales, comparaciones e intercambios."""

    num_pasos: int = 11
    pasos: Dict[int, RegistroPaso] = field(init=False)
    extra_mult_div: int = 0
    extra_sum_resta: int = 0
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
        mult_div_pasos = sum(r.mult_div for r in self.pasos.values())
        sum_resta_pasos = sum(r.sum_resta for r in self.pasos.values())
        comp_pasos = sum(r.comparaciones for r in self.pasos.values())
        inter_pasos = sum(r.intercambios for r in self.pasos.values())

        total_mult_div = mult_div_pasos + self.extra_mult_div
        total_sum_resta = sum_resta_pasos + self.extra_sum_resta
        total_flops = total_mult_div + total_sum_resta
        total_comparaciones = comp_pasos + self.extra_comparaciones
        total_intercambios = inter_pasos + self.extra_intercambios

        return {
            "mult_div": total_mult_div,
            "sum_resta": total_sum_resta,
            "total_flops": total_flops,
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
                    "comparaciones": r.comparaciones,
                    "intercambios": r.intercambios,
                }
            )
        if (
            self.extra_mult_div > 0
            or self.extra_sum_resta > 0
            or self.extra_comparaciones > 0
            or self.extra_intercambios > 0
        ):
            filas.append(
                {
                    "etapa": "Sustitución regresiva",
                    "mult_div": self.extra_mult_div,
                    "sum_resta": self.extra_sum_resta,
                    "flops_paso": self.extra_mult_div + self.extra_sum_resta,
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
        "- **1 flop = 1 operación aritmética elemental** (`+`, `-`, `*`, `/`) ejecutada en el código.",
        "- Se contabilizan por separado las **multiplicaciones/divisiones** y las **sumas/restas**.",
        "- Las **comparaciones** e **intercambios** asociados al pivoteo se registran aparte como sobrecosto de estabilidad.",
        "- Se incluyen la actualización del vector `b` y la sustitución regresiva.",
        f"- Referencia teórica $2n^3/3$ para $n={n}$: `{cota_teorica_2n3_3:.2f}` (referencia práctica del enunciado: ~1486 flops + ~66 comparaciones con pivoteo).",
        "",
    ]

    for nombre, contador in reportes_por_estrategia.items():
        tot = contador.totales()
        lineas.extend(
            [
                f"## Estrategia: {nombre}",
                "",
                "| Etapa / Paso | Mult / Div | Sum / Resta | FLOPs paso | Comparaciones | Intercambios |",
                "|---|---:|---:|---:|---:|---:|",
            ]
        )
        for fila in contador.filas_tabla():
            lineas.append(
                f"| {fila['etapa']} | {fila['mult_div']} | {fila['sum_resta']} | "
                f"{fila['flops_paso']} | {fila['comparaciones']} | {fila['intercambios']} |"
            )
        lineas.extend(
            [
                f"| **Total** | **{tot['mult_div']}** | **{tot['sum_resta']}** | "
                f"**{tot['total_flops']}** | **{tot['comparaciones']}** | **{tot['intercambios']}** |",
                "",
            ]
        )

    ruta.write_text("\n".join(lineas), encoding="utf-8")