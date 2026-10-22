from __future__ import annotations

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from app.api.routes import router


app = FastAPI(
    title="Chem Sandbox API",
    version="0.26.0",
    description=(
        "Blind molecular discovery, evidence-gated pathways, conditional kinetics, "
        "supported Cantera reactors, multiphase boundaries, explicit catalyst cycles, "
        "steady-state process flowsheets, recirculating prebiotic apparatus models, "
        "and leakage-resistant post-run validation."
    ),
)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
app.include_router(router)


@app.get("/")
def root() -> dict[str, str]:
    return {"name": "Chem Sandbox API", "docs": "/docs"}
