from __future__ import annotations

import json

import cantera as ct
import idaes

from app.operations.capabilities import operation_capabilities


def main() -> None:
    # Runtime smoke test only. This inert nitrogen state is not a reaction
    # recipe and is never visible to blind discovery.
    gas = ct.Solution("gri30.yaml")
    gas.TPX = 298.15, ct.one_atm, "N2:1"
    reactor = ct.IdealGasConstPressureReactor(gas, clone=True)
    network = ct.ReactorNet([reactor])
    network.advance(1e-6)
    payload = operation_capabilities()
    payload["smoke_test"] = {
        "cantera_temperature_K": reactor.T,
        "cantera_pressure_Pa": reactor.phase.P,
        "idaes_version": idaes.__version__,
        "status": "passed",
    }
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
