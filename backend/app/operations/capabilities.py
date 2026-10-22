from __future__ import annotations

from importlib import import_module
from importlib.metadata import PackageNotFoundError, version
from typing import Any


def _python_package(module: str, distribution: str, purpose: str) -> dict[str, Any]:
    try:
        import_module(module)
        package_version = version(distribution)
    except (ImportError, PackageNotFoundError) as error:
        return {
            "id": module,
            "available": False,
            "version": None,
            "purpose": purpose,
            "detail": f"{type(error).__name__}: {error}",
        }
    return {
        "id": module,
        "available": True,
        "version": package_version,
        "purpose": purpose,
        "detail": "Installed in the Chem Sandbox virtual environment",
    }


def operation_capabilities() -> dict[str, Any]:
    return {
        "installed": [
            _python_package(
                "cantera",
                "Cantera",
                "Homogeneous kinetics and batch/CSTR/PFR reactor networks after a mechanism is discovered",
            ),
            _python_package(
                "idaes",
                "idaes-pse",
                "Flowsheets, unit operations, optimization, and process-level constraints",
            ),
            _python_package(
                "libsbml",
                "python-libsbml",
                "SBML Level 3 model exchange, validation, and network inspection",
            ),
            _python_package(
                "roadrunner",
                "libroadrunner",
                "Independent CVODE time-course cross-checks for staged biochemical networks",
            ),
        ],
        "isolated_future_engines": [
            {
                "id": "rmg",
                "available": False,
                "purpose": "Automated detailed kinetic mechanism generation",
                "required_runtime": "Official Docker or Linux/WSL environment",
            },
            {
                "id": "reaktoro",
                "available": False,
                "purpose": "Aqueous/electrolyte equilibrium and reactive transport",
                "required_runtime": "Dedicated Conda/Linux environment",
            },
            {
                "id": "catmap",
                "available": False,
                "purpose": "Heterogeneous catalytic microkinetics",
                "required_runtime": "Dedicated compatible Python/Linux environment",
            },
        ],
        "boundary": (
            "Operation engines consume a frozen, provenance-labelled mechanism after blind discovery; "
            "their mechanism files and reference data do not feed candidate generation."
        ),
    }
