// Typed client for the optional Python model (Codes/python/api.py), used by the "Full report" tab.
// In development Vite proxies /api to 127.0.0.1:8000; the built single-file site calls it directly.
const BASE = import.meta.env.DEV ? "" : "http://127.0.0.1:8000";

export type Issue = { level: "error" | "warning" | "info"; code: string; message: string };

export type Point = {
  T_C: number;
  status: string;
  x: number;
  gamma: number;
  x_ideal: number;
  mg_per_100g: number | null;
  miscible: boolean;
};

export type CurveOut = { solute: string; oil: string; Tm_C: number; points: Point[] };

export type MixtureRow = {
  compound: string;
  loading_mg_per_100g: number | null;
  mg_per_100g: number | null;
  alone_mg_per_100g: number | null;
  change_vs_alone_percent: number | null;
  x: number;
  gamma: number;
  gamma_alone: number;
  saturated: boolean;
  state: string;
  extracted_fraction: number | null;
};

export type MixtureOut = {
  oil: string;
  T_C: number;
  rows: MixtureRow[];
  notes: string[];
  newton_iterations: number;
  residual_history: number[];
  issues: Issue[];
};

export type BatchRow = {
  compound: string;
  available_g: number;
  capacity_alone_g: number | null;
  dissolved_g: number;
  undissolved_g: number;
  dissolved_percent: number | null;
  limited_by: "solubility" | "amount added";
  mole_fraction: number;
  wt_percent: number | null;
  ppm_w: number | null;
  mg_per_100g_oil: number | null;
};

export type BatchOut = { oil: string; T_C: number; oil_g: number; rows: BatchRow[] };

export type PredictIn = {
  solute: string;
  oil: string;
  T_C: number[];
  gamma?: number;
  gamma_T_C?: number;
  oil_grade?: string;
  monte_carlo?: number;
  measurements?: { value: number; unit: string; T_C: number; solvent?: string }[];
};

export type PredictOut = {
  report: string;
  issues: Issue[];
  predictions: {
    T_C: number;
    route: string;
    status: string;
    x: number;
    mg_per_100g: number | null;
    gamma: number;
    ci95_mg_per_100g: (number | null)[];
    issues: Issue[];
  }[];
};

export class ApiError extends Error {}

async function call<T>(path: string, init?: RequestInit): Promise<T> {
  let res: Response;
  try {
    res = await fetch(BASE + path, {
      ...init,
      headers: { "content-type": "application/json", ...(init?.headers ?? {}) },
    });
  } catch {
    throw new ApiError("Could not reach the Python model at 127.0.0.1:8000.");
  }
  if (!res.ok) {
    let msg = `${res.status} ${res.statusText}`;
    try {
      const body = await res.json();
      if (body?.detail) msg = typeof body.detail === "string" ? body.detail : JSON.stringify(body.detail);
    } catch {
      // Not JSON: Vite's proxy error page (API down) or a Python crash.
      if (res.status >= 500) msg = "Python model is not reachable.";
    }
    throw new ApiError(msg);
  }
  return (await res.json()) as T;
}

const qs = (o: Record<string, string | number>) =>
  new URLSearchParams(Object.entries(o).map(([k, v]) => [k, String(v)])).toString();

export const api = {
  health: () => call<{ ok: boolean }>("/api/health"),
  curve: (solute: string, oil: string) => call<CurveOut>(`/api/curve?${qs({ solute, oil, step: 5 })}`),
  point: (solute: string, oil: string, T: number) =>
    call<{ points: Point[] }>(`/api/solubility?${qs({ solute, oil, T })}`).then((r) => r.points[0]),
  flags: (solute: string, oil: string, T: number, grade: string) =>
    call<{ issues: Issue[] }>(`/api/flags?${qs({ solute, oil, T, grade })}`).then((r) => r.issues),
  mixture: (body: {
    oil: string;
    T_C: number;
    loadings: Record<string, number> | null;
    compounds?: string[];
  }) => call<MixtureOut>("/api/mixture", { method: "POST", body: JSON.stringify(body) }),
  batch: (body: { oil: string; T_C: number; oil_g: number; added_g: Record<string, number> }) =>
    call<BatchOut>("/api/batch", { method: "POST", body: JSON.stringify(body) }),
  predict: (body: PredictIn) =>
    call<PredictOut>("/api/predict", { method: "POST", body: JSON.stringify(body) }),
};
