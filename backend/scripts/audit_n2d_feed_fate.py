"""Read-only N(2D) first-collision audit on the held-out-assay Barth feed.

The feed/geometry card contains no measured HCN amount. This audit neither
assigns an N2 dissociation branching fraction nor mutates apparatus chemistry.
"""

from __future__ import annotations

import json
from typing import Any

from app.prebiotic.n2d_fate import n2d_fate_diagnostic
from scripts.probe_barth_feed_pathway import make_probe_payload


def audit_barth_initial_n2d_fate() -> dict[str, Any]:
    payload = make_probe_payload(duration_s=60.0)
    gas = {
        item["species_id"]: float(item["amount"]["value"])
        for item in payload["initial_inventory"]
        if item["compartment"] == "gas_loop"
    }
    return {
        "scenario": "Barth 220725001 initial reported gas fractions; hypothetical 298.15 K, 0.95 L headspace and 1 bar",
        "feed_card_contains_hcn_assay": False,
        "electron_dissociation_state_branch_assigned": False,
        "n2d_fate": n2d_fate_diagnostic(
            gas,
            gas_volume_L=float(payload["apparatus"]["gas_loop_volume"]["value"]),
            gas_temperature_K=float(payload["apparatus"]["discharge_temperature"]["value"]),
        ),
    }


def main() -> None:
    print(json.dumps(audit_barth_initial_n2d_fate(), indent=2))


if __name__ == "__main__":
    main()
