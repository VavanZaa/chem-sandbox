from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[3]
LOCAL_BIN = PROJECT_ROOT / ".tools" / "idaes-extensions"


def idaes_capabilities() -> dict[str, Any]:
    try:
        idaes_version = version("idaes-pse")
        framework_available = True
    except PackageNotFoundError:
        idaes_version = None
        framework_available = False
    ipopt = LOCAL_BIN / "ipopt.exe"
    return {
        "idaes_available": framework_available,
        "idaes_version": idaes_version,
        "ipopt_available": ipopt.is_file(),
        "ipopt_path": str(ipopt) if ipopt.is_file() else None,
        "native_balance_solver": "fixed-point recycle solver",
        "note": "The native solver remains available when IDAES/IPOPT is absent; cross-engine agreement is reported only when IPOPT runs.",
    }


def solve_idaes_algebraic(
    components: list[str],
    fresh: dict[str, float],
    stoichiometry: dict[str, float],
    limiting_species: str,
    conversion: float,
    recoveries: dict[str, float],
    recycle_fraction: float,
) -> dict[str, Any]:
    capability = idaes_capabilities()
    if not capability["idaes_available"] or not capability["ipopt_available"]:
        return {"status": "unavailable", "capabilities": capability}
    try:
        from idaes.core import FlowsheetBlock
        from pyomo.environ import ConcreteModel, Constraint, Set, SolverFactory, Var, value

        model = ConcreteModel()
        model.fs = FlowsheetBlock(dynamic=False)
        model.fs.components = Set(initialize=components)
        for name in ("mixed", "reactor_outlet", "product", "recycle", "purge"):
            setattr(model.fs, name, Var(model.fs.components, bounds=(0, None), initialize=1.0))
        model.fs.extent = Var(bounds=(0, None), initialize=1.0)

        @model.fs.Constraint(model.fs.components)
        def mixer_balance(block, component):
            return block.mixed[component] == fresh.get(component, 0.0) + block.recycle[component]

        limiting_coefficient = -stoichiometry[limiting_species]

        @model.fs.Constraint()
        def conversion_constraint(block):
            return block.extent == conversion * block.mixed[limiting_species] / limiting_coefficient

        @model.fs.Constraint(model.fs.components)
        def reactor_balance(block, component):
            return block.reactor_outlet[component] == block.mixed[component] + stoichiometry.get(component, 0.0) * block.extent

        @model.fs.Constraint(model.fs.components)
        def product_split(block, component):
            return block.product[component] == recoveries[component] * block.reactor_outlet[component]

        @model.fs.Constraint(model.fs.components)
        def recycle_split(block, component):
            return block.recycle[component] == recycle_fraction * (block.reactor_outlet[component] - block.product[component])

        @model.fs.Constraint(model.fs.components)
        def purge_split(block, component):
            return block.purge[component] == (1.0 - recycle_fraction) * (
                block.reactor_outlet[component] - block.product[component]
            )

        solver = SolverFactory("ipopt", executable=capability["ipopt_path"])
        result = solver.solve(model, tee=False)
        termination = str(result.solver.termination_condition)
        if termination.lower() not in {"optimal", "locallyoptimal", "globallyoptimal"}:
            return {"status": "failed", "termination_condition": termination, "capabilities": capability}
        return {
            "status": "completed",
            "termination_condition": termination,
            "extent_mol_s": float(value(model.fs.extent)),
            "mixed_mol_s": {key: float(value(model.fs.mixed[key])) for key in components},
            "product_mol_s": {key: float(value(model.fs.product[key])) for key in components},
            "purge_mol_s": {key: float(value(model.fs.purge[key])) for key in components},
            "capabilities": capability,
        }
    except Exception as error:
        return {"status": "failed", "error": f"{type(error).__name__}: {error}", "capabilities": capability}
