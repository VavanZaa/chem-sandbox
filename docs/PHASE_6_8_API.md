# Phase 6–8 API guide

These endpoints are post-discovery tools. They accept explicit mechanisms, operating evidence, or frozen predictions. None can add products, recipes, yields, or reference outcomes to blind molecular discovery.

## Phase 6: catalyst microkinetics

`POST /api/catalysis/simulate` integrates a reversible mean-field surface mechanism. A request supplies gas and surface species, exactly one vacant-site pseudo-species, initial coverages, fixed gas activities, active-site amount, and elementary adsorption/desorption/surface/diffusion steps. Every step requires sourced forward barrier, reverse barrier, and reaction Gibbs energy.

The schema enforces element, charge, and site conservation. Every elementary step must change at least one non-vacant bound intermediate; an overall reaction in which a catalyst is merely written on both sides is rejected. The solver checks microscopic reversibility, integrates coverages with BDF, and reports site closure, reaction rates, gas net changes, active fraction, and target turnover frequencies. Optional first-order deactivation and a sourced 0–1 transport effectiveness factor are explicit modifiers, not inferred corrections.

`POST /api/catalysis/catmap/export` returns a CatMAP setup fragment with reaction expressions, gas pressures, one site definition, and multisite occupancies. Formation energies and scaling relations remain a separate CatMAP input and are never invented.

Current limits: fixed gas reservoirs, one site family, mean-field activities, and no lateral interactions unless already embedded in the supplied free energies. Couple to a reactor model before interpreting gas production as bulk yield.

## Phase 7: process flowsheet

`GET /api/process/capabilities` reports IDAES and project-local IPOPT availability.

`POST /api/process/simulate` solves a steady fresh-feed → reactor → component separator → product/purge/recycle loop. Inputs include a balanced reaction, limiting species, sourced single-pass conversion, sourced recovery for every component, recycle fraction, optional utility per extent, and an optional internal-flow limit.

The response contains every stream, reaction extent, overall conversion, product rate, utility duty, recycle convergence, elemental balance, uncertainty endpoint scenarios, optimization status, and an independent IDAES/IPOPT algebraic cross-check when available. Recycle optimization searches only the submitted bounds and maximizes the submitted target flow under the internal-flow constraint; it cannot change chemistry.

Current limits: one reaction and one recycle topology, with no rigorous property package, phase equilibrium, pressure drop, equipment sizing, costing, or relief calculation. Those omissions are reported in every result.

## Phase 8: post-run calibration and validation

`POST /api/validation/freeze` hashes a prediction-only set containing means and standard deviations. The returned `prediction_sha256` binds all condition IDs and values.

`POST /api/validation/calibrate` accepts observations only after freeze and fits either an affine or offset-only calibration with an uncertainty scale. Its artifact records all training IDs and carries its own SHA-256.

`POST /api/validation/evaluate` computes raw and calibrated MAE, RMSE, bias, Gaussian mean negative log likelihood, standardized residuals, and 50/80/95% interval coverage. It rejects altered prediction/calibration artifacts and rejects any held-out condition ID used during calibration.

The reproducible Phase 3–5 chemistry gate is:

```powershell
.\.venv\Scripts\python.exe .\backend\scripts\run_phase_3_5_validation.py
```

It writes a reference-free prediction file first, then imports evaluator-only checks and writes a separately hashed unlocked evaluation under `backend/data/validation/phase_3_5_latest/`.
