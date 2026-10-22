from __future__ import annotations

import hashlib
import json
import os
from typing import Any

import httpx

from .schema import ElectrolyteEquilibriumRequest


class ReaktoroUnavailableError(RuntimeError):
    pass


def reaktoro_status() -> dict[str, Any]:
    url = os.getenv("CHEM_SANDBOX_REAKTORO_URL")
    if not url:
        return {
            "available": False,
            "mode": "isolated_service",
            "reason": "CHEM_SANDBOX_REAKTORO_URL is not configured",
            "fallback": "Henry-law gas-liquid transfer only; no electrolyte speciation is fabricated",
        }
    try:
        response = httpx.get(f"{url.rstrip('/')}/health", timeout=2.0)
        response.raise_for_status()
        remote = response.json()
    except Exception as error:
        return {
            "available": False,
            "mode": "isolated_service",
            "url": url,
            "reason": f"{type(error).__name__}: {error}",
        }
    if remote.get("status") != "healthy" or remote.get("engine") != "Reaktoro":
        return {
            "available": False,
            "mode": "isolated_service",
            "url": url,
            "reason": "The configured service did not report a healthy Reaktoro runtime",
            "remote": remote,
        }
    return {"available": True, "mode": "isolated_service", "url": url, "remote": remote}


def equilibrate_with_reaktoro(request: ElectrolyteEquilibriumRequest) -> dict[str, Any]:
    status = reaktoro_status()
    if not status.get("available"):
        raise ReaktoroUnavailableError(status.get("reason", "Reaktoro service is unavailable"))
    payload = request.model_dump(mode="json")
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    response = httpx.post(
        f"{status['url'].rstrip('/')}/equilibrate",
        json={"request_sha256": digest, "request": payload},
        timeout=120.0,
    )
    if response.status_code == 422:
        try:
            detail = response.json().get("detail", response.text)
        except ValueError:
            detail = response.text
        raise ValueError(f"Reaktoro rejected the equilibrium request: {detail}")
    response.raise_for_status()
    result = response.json()
    result["request_sha256"] = digest
    result["boundary"] = "Isolated post-discovery equilibrium service; results cannot feed blind candidate generation."
    return result
