"""
HTTP API for the website. Thin wrappers around the unchanged `flavorsol` package.

Run (from this folder):
    uvicorn api:app --host 127.0.0.1 --port 8000 --reload

Every endpoint calls the same functions as the CLI / notebook:
    /api/solubility  ->  python -m flavorsol solubility ...   (sle.solubility)
    /api/mixture     ->  python -m flavorsol mixture ...      (mixture.solve_mixture)
    /api/batch       ->  python -m flavorsol batch ...        (operating.batch)
    /api/predict     ->  python -m flavorsol predict x.json   (engine.predict)
Interactive docs: http://127.0.0.1:8000/docs
"""
import math
import warnings
from typing import Optional

import numpy as np
from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

from flavorsol import checks
from flavorsol.data.compounds import FLAVOURS
from flavorsol.data.oils import OILS, get_oil
from flavorsol.mixture import solve_mixture
from flavorsol.sle import solubility
from flavorsol.units import binary_mg_per_100g, to_kelvin

warnings.filterwarnings("ignore", message="UNIFAC: no published parameters")

MIX_ORDER = ("vanillin", "curcumin", "piperine")   # same order as the CLI `mixture` command

app = FastAPI(title="flavorsol API", version="1.0",
              description="Solubility of spice flavour compounds in cooking oils (CL 202, Group 4A)")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:8080", "http://127.0.0.1:8080"],
    allow_methods=["*"],
    allow_headers=["*"],
)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def clean(v):
    """numpy / inf / nan -> plain JSON values (inf and nan become null)."""
    if isinstance(v, dict):
        return {str(k): clean(x) for k, x in v.items()}
    if isinstance(v, (list, tuple, np.ndarray)):
        return [clean(x) for x in v]
    if isinstance(v, (np.bool_, bool)):
        return bool(v)
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        f = float(v)
        return f if math.isfinite(f) else None
    return v


def _solute(name):
    key = name.strip().lower()
    if key not in FLAVOURS:
        raise HTTPException(400, f"Unknown compound '{name}'. Known: {sorted(FLAVOURS)}")
    return FLAVOURS[key]


def _oil(name):
    try:
        return get_oil(name)
    except KeyError as e:
        raise HTTPException(400, str(e.args[0]))


def _check_T(T_C):
    if not (-50 <= T_C <= 300):
        raise HTTPException(400, f"T = {T_C} C is outside -50..300 C")


def _point(s, oil, T_C):
    r = solubility(s, oil, float(to_kelvin(T_C)))
    mg = binary_mg_per_100g(r.x, s, oil)
    return {"T_C": T_C, "status": r.status, "x": r.x, "gamma": r.gamma,
            "x_ideal": r.x_ideal, "mg_per_100g": mg,
            "miscible": r.status == "miscible" or not math.isfinite(mg)}


def _issues(lst):
    return [{"level": i.level, "code": i.code, "message": i.message} for i in lst]


def _run(fn):
    try:
        return fn()
    except HTTPException:
        raise
    except (ValueError, RuntimeError, KeyError, TypeError) as e:
        raise HTTPException(400, str(e))


# ---------------------------------------------------------------------------
# endpoints
# ---------------------------------------------------------------------------
@app.get("/api/health")
def health():
    return {"ok": True}


@app.get("/api/meta")
def meta():
    return clean({
        "compounds": [{"id": k, "M": c.M, "Tm_C": c.Tm - 273.15 if c.Tm else None,
                       "dHfus_kJ": c.dHfus / 1e3 if c.dHfus else None, "formula": c.formula}
                      for k, c in FLAVOURS.items()],
        "oils": [{"id": k, "M": o.M} for k, o in OILS.items()],
        "mixture_order": MIX_ORDER,
    })


@app.get("/api/solubility")
def api_solubility(solute: str, oil: str, T: list[float] = Query(default=[25, 60, 100])):
    """One compound in one oil at one or more temperatures (deg C)."""
    s, o = _solute(solute), _oil(oil)
    for t in T:
        _check_T(t)
    return clean(_run(lambda: {"solute": s.name, "oil": o.name,
                               "points": [_point(s, o, t) for t in T]}))


@app.get("/api/curve")
def api_curve(solute: str, oil: str, step: float = 5.0):
    """Solubility from 25 to 200 deg C, for the chart."""
    s, o = _solute(solute), _oil(oil)
    step = min(max(step, 1.0), 25.0)
    temps = [float(t) for t in np.arange(25.0, 200.0 + 1e-9, step)]
    return clean(_run(lambda: {"solute": s.name, "oil": o.name, "Tm_C": s.Tm - 273.15,
                               "points": [_point(s, o, t) for t in temps]}))


@app.get("/api/flags")
def api_flags(solute: str, oil: str, T: float, grade: str = "refined"):
    """Boundary / degradation / smoke-point checks from flavorsol.checks."""
    s = _solute(solute)
    return {"issues": _issues(checks.check_temperature(T, s.name, oil, grade, s.Tm))}


class MixtureIn(BaseModel):
    oil: str = "sunflower"
    T_C: float = 60.0
    # mg added per 100 g oil; None = every compound in excess (solid always present)
    loadings: Optional[dict[str, float]] = None
    compounds: list[str] = Field(default_factory=lambda: list(MIX_ORDER))


@app.post("/api/mixture")
def api_mixture(req: MixtureIn):
    """Several compounds in the same oil (notebook section 4)."""
    _check_T(req.T_C)
    o = _oil(req.oil)
    names = [n.strip().lower() for n in req.compounds]
    if not names:
        raise HTTPException(400, "Choose at least one compound.")
    sols = [_solute(n) for n in names]
    loads = None
    if req.loadings is not None:
        loads = [float(req.loadings.get(n, 0.0)) for n in names]
        if any(L < 0 for L in loads):
            raise HTTPException(400, "Loadings must be zero or positive.")
        if all(L == 0 for L in loads):
            raise HTTPException(400, "Add some of at least one compound.")
        # a zero loading breaks the log-space solver: drop those compounds
        keep = [i for i, L in enumerate(loads) if L > 0]
        names, sols, loads = [names[i] for i in keep], [sols[i] for i in keep], [loads[i] for i in keep]
    T = float(to_kelvin(req.T_C))

    def go():
        mix = solve_mixture(sols, o, T, loads)
        rows = []
        for i, n in enumerate(mix.solutes):
            alone = _point(sols[i], o, req.T_C)
            mg, mg_alone = float(mix.mg_per_100g[i]), alone["mg_per_100g"]
            change = (100 * (mg / mg_alone - 1)
                      if mix.saturated[i] and math.isfinite(mg_alone) and mg_alone > 0 else None)
            rows.append({
                "compound": n,
                "loading_mg_per_100g": loads[i] if loads else None,
                "mg_per_100g": mg,
                "alone_mg_per_100g": mg_alone,
                "change_vs_alone_percent": change,
                "x": mix.x[i],
                "gamma": mix.gamma[i],
                "gamma_alone": alone["gamma"],
                "saturated": mix.saturated[i],
                "state": "saturated (solid left)" if mix.saturated[i] else "fully dissolved",
                "extracted_fraction": mix.extracted_fraction[i],
            })
        return {"oil": o.name, "T_C": req.T_C, "rows": rows, "notes": mix.notes,
                "newton_iterations": mix.newton_iterations,
                "residual_history": mix.residual_history,
                "issues": _issues(checks.check_temperature(req.T_C, None, req.oil))}

    return clean(_run(go))


class BatchIn(BaseModel):
    oil: str = "coconut"
    T_C: float = 60.0
    oil_g: float = 15.0
    added_g: dict[str, float] = Field(default_factory=lambda: {"curcumin": 0.047, "piperine": 0.05})


@app.post("/api/batch")
def api_batch(req: BatchIn):
    """Grams of oil + grams of each compound -> dissolved / left solid (notebook section 5)."""
    from flavorsol.operating import batch
    _check_T(req.T_C)
    _oil(req.oil)
    if req.oil_g <= 0:
        raise HTTPException(400, "Oil mass must be positive.")
    added = {k.strip().lower(): float(v) for k, v in req.added_g.items() if v and v > 0}
    if not added:
        raise HTTPException(400, "Add some of at least one compound (grams).")
    for k in added:
        _solute(k)
    rows = _run(lambda: batch(req.oil, req.T_C, req.oil_g, added))
    return clean({"oil": req.oil, "T_C": req.T_C, "oil_g": req.oil_g, "rows": rows})


class MeasurementIn(BaseModel):
    value: float
    unit: str
    T_C: float
    solvent: Optional[str] = None
    kind: str = "solubility"
    sigma_rel: Optional[float] = None
    sigma_T: float = 0.5
    x_solution: Optional[float] = None
    p_sub_Pa: Optional[float] = None
    source: str = ""


class PredictIn(BaseModel):
    """Same fields as the JSON files in inputs/ (only solute, oil, T_C are mandatory)."""
    solute: str
    oil: str
    T_C: list[float]
    M: Optional[float] = None
    Tm: Optional[float] = None
    dHfus: Optional[float] = None
    dCp: Optional[float] = None
    groups: Optional[dict[str, int]] = None
    gamma: Optional[float] = None
    gamma_T_C: Optional[float] = None
    measurements: list[MeasurementIn] = Field(default_factory=list)
    oil_composition: Optional[dict[str, float]] = None
    oil_grade: str = "refined"
    uncertainties: dict[str, float] = Field(default_factory=dict)
    loading_mg_per_100g: Optional[float] = None
    monte_carlo: int = 0
    T_extrapolation: str = "unifac"


@app.post("/api/predict")
def api_predict(req: PredictIn):
    """Full workflow with routes R1-R5, checks and uncertainty (notebook section 6)."""
    from flavorsol.engine import Measurement, Request, format_report, predict
    d = req.model_dump()
    d["measurements"] = [Measurement(**m) for m in d["measurements"]]
    d["monte_carlo"] = min(max(int(d["monte_carlo"] or 0), 0), 5000)

    def go():
        preds, issues = predict(Request(**d))
        return {
            "report": format_report(preds, issues),
            "issues": _issues(issues),
            "predictions": [{
                "T_C": p.T_C, "route": p.route, "status": p.status, "x": p.x,
                "mg_per_100g": p.mg_per_100g, "gamma": p.gamma,
                "ci95_mg_per_100g": list(p.ci95_mg_per_100g),
                "issues": _issues(p.issues),
            } for p in preds],
        }

    return clean(_run(go))
