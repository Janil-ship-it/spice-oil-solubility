/*
 * Browser copy of the project's thermodynamic model (Codes/python/flavorsol).
 *
 *  - Single compound: solubility curves exported from the Python model
 *    (tools/export_web_data.py, 1 °C steps, 25-200 °C), interpolated in ln(mg).
 *  - Mixtures: UNIFAC + damped Newton-Raphson + active-set loop, a line-by-line
 *    port of flavorsol/unifac.py and flavorsol/mixture.py, using the same
 *    parameters (model.json). Results agree with Python to better than 0.1 %.
 *  - Saturation point: same rules as flavorsol/saturation.py.
 */
import model from "@/data/model.json";

export type Compound = "curcumin" | "vanillin" | "piperine";
export type Oil = "sunflower" | "coconut";

export const COMPOUNDS: { id: Compound; label: string; spice: string; note: string }[] = [
  { id: "curcumin", label: "Curcumin", spice: "Turmeric", note: "Yellow colour of turmeric. Melts at 183 °C." },
  { id: "vanillin", label: "Vanillin", spice: "Vanilla", note: "Sweet aroma of vanilla. Melts at 82 °C." },
  { id: "piperine", label: "Piperine", spice: "Black pepper", note: "Pungency of black pepper. Melts at 133 °C." },
];
export const OILS: { id: Oil; label: string }[] = [
  { id: "sunflower", label: "Sunflower" },
  { id: "coconut", label: "Coconut" },
];
export const MIX_ORDER: Compound[] = ["vanillin", "curcumin", "piperine"];
export const LABEL: Record<Compound, string> = { vanillin: "Vanillin", curcumin: "Curcumin", piperine: "Piperine" };

type CompoundData = { M: number; Tm: number; dHfus: number; dCp: number; formula: string; groups: Record<string, number> };
type OilData = { M: number; groups: Record<string, number> };
type Curve = { mg: (number | null)[]; x: number[]; gamma: number[]; status: string };

const M = model as unknown as {
  compounds: Record<Compound, CompoundData>;
  oils: Record<Oil, OilData>;
  subgroups: Record<string, { sid: number; main: number; R: number; Q: number }>;
  interactions: Record<string, number>;
  curves: Record<Oil, Record<Compound, Curve>>;
};

export const T_MIN = 25;
export const T_MAX = 200;
const R_GAS = 8.314462618;
const N = T_MAX - T_MIN + 1;

export const compoundData = (c: Compound) => M.compounds[c];
export const oilData = (o: Oil) => M.oils[o];
export const TmC = (c: Compound) => M.compounds[c].Tm - 273.15;

/* ------------------------------------------------------------------ single compound */
export type Status = "SLE" | "LLE-limited" | "miscible";
const ST: Record<string, Status> = { S: "SLE", L: "LLE-limited", M: "miscible" };

export type Point = { T: number; mgPer100g: number | null; x: number; gamma: number; status: Status; miscible: boolean };

/** Solubility of one compound in one oil at T (°C), from the Python curve. */
export function point(oil: Oil, c: Compound, T: number): Point {
  const cv = M.curves[oil][c];
  const t = Math.min(T_MAX, Math.max(T_MIN, T));
  const i0 = Math.min(N - 2, Math.floor(t - T_MIN));
  const f = t - T_MIN - i0;
  const s0 = ST[cv.status[i0]];
  const s1 = ST[cv.status[i0 + 1]];
  const molten = t + 273.15 >= M.compounds[c].Tm;
  const status: Status = f < 1e-9 ? s0 : molten ? s1 : s0 === "miscible" ? s1 : s0;
  const a = cv.mg[i0];
  const b = cv.mg[i0 + 1];
  const lerp = (p: number, q: number) => p + (q - p) * f;
  if (status === "miscible" || a == null) {
    return { T: t, mgPer100g: null, x: 1, gamma: 1, status: "miscible", miscible: true };
  }
  const mg = b == null ? a : Math.exp(lerp(Math.log(a), Math.log(b)));
  const x = b == null ? cv.x[i0] : Math.exp(lerp(Math.log(cv.x[i0]), Math.log(cv.x[i0 + 1])));
  return { T: t, mgPer100g: mg, x, gamma: lerp(cv.gamma[i0], cv.gamma[i0 + 1]), status, miscible: false };
}

export function curve(oil: Oil, c: Compound, step = 1) {
  const pts: { T: number; mg: number | null }[] = [];
  for (let t = T_MIN; t <= T_MAX; t += step) pts.push({ T: t, mg: M.curves[oil][c].mg[t - T_MIN] });
  return pts;
}

/** Temperature at which `loading` mg/100 g is exactly saturated (inverse of the curve). */
export function saturationTemperature(oil: Oil, c: Compound, loading: number): { T: number | null; kind: "root" | "below_range" | "above_range" | "none" } {
  if (!(loading > 0)) return { T: null, kind: "none" };
  const cv = M.curves[oil][c];
  const val = (i: number) => (cv.mg[i] == null ? Infinity : (cv.mg[i] as number));
  if (val(0) >= loading) return { T: null, kind: "below_range" };
  for (let i = 1; i < N; i++) {
    const b = val(i);
    if (b >= loading) {
      const a = val(i - 1);
      if (!isFinite(b)) return { T: Math.min(T_MIN + i, TmC(c)), kind: "root" };
      const f = (Math.log(loading) - Math.log(a)) / (Math.log(b) - Math.log(a));
      return { T: T_MIN + i - 1 + f, kind: "root" };
    }
  }
  return { T: null, kind: "above_range" };
}

/* ------------------------------------------------------------------ UNIFAC */
export class Unifac {
  nc: number;
  ng: number;
  nu: number[][];
  R: number[];
  Q: number[];
  A: number[][];
  r: number[];
  q: number[];
  private pureCache = new Map<number, number[][]>();

  constructor(groupsList: Record<string, number>[]) {
    const names = Array.from(new Set(groupsList.flatMap((g) => Object.keys(g)))).sort(
      (a, b) => M.subgroups[a].sid - M.subgroups[b].sid,
    );
    this.nc = groupsList.length;
    this.ng = names.length;
    this.nu = groupsList.map((gl) => names.map((g) => gl[g] ?? 0));
    this.R = names.map((g) => M.subgroups[g].R);
    this.Q = names.map((g) => M.subgroups[g].Q);
    const main = names.map((g) => M.subgroups[g].main);
    this.A = main.map((a) => main.map((b) => (a === b ? 0 : (M.interactions[`${a},${b}`] ?? 0))));
    this.r = this.nu.map((row) => row.reduce((s, v, k) => s + v * this.R[k], 0));
    this.q = this.nu.map((row) => row.reduce((s, v, k) => s + v * this.Q[k], 0));
  }

  private psi(T: number) {
    return this.A.map((row) => row.map((a) => Math.exp(-a / T)));
  }

  private lnGammaGroups(X: number[], Psi: number[][]) {
    const ng = this.ng;
    const th = X.map((v, k) => v * this.Q[k]);
    const tot = th.reduce((s, v) => s + v, 0);
    for (let k = 0; k < ng; k++) th[k] /= tot;
    const S = new Array(ng).fill(0);
    for (let k = 0; k < ng; k++) for (let m = 0; m < ng; m++) S[k] += th[m] * Psi[m][k];
    const out = new Array(ng);
    for (let k = 0; k < ng; k++) {
      let sum = 0;
      for (let m = 0; m < ng; m++) sum += (th[m] / S[m]) * Psi[k][m];
      out[k] = this.Q[k] * (1 - Math.log(S[k]) - sum);
    }
    return out as number[];
  }

  private pure(T: number) {
    const key = Math.round(T * 1e10) / 1e10;
    let v = this.pureCache.get(key);
    if (!v) {
      const Psi = this.psi(T);
      v = this.nu.map((row) => {
        const s = row.reduce((a, b) => a + b, 0);
        return this.lnGammaGroups(row.map((n) => n / s), Psi);
      });
      if (this.pureCache.size > 2000) this.pureCache.clear();
      this.pureCache.set(key, v);
    }
    return v;
  }

  lnGamma(xIn: number[], T: number) {
    const sx = xIn.reduce((a, b) => a + b, 0);
    const x = xIn.map((v) => v / sx);
    const xr = x.reduce((s, v, i) => s + v * this.r[i], 0);
    const xq = x.reduce((s, v, i) => s + v * this.q[i], 0);
    const Psi = this.psi(T);
    const Xg = new Array(this.ng).fill(0);
    for (let i = 0; i < this.nc; i++) for (let k = 0; k < this.ng; k++) Xg[k] += x[i] * this.nu[i][k];
    const sg = Xg.reduce((a, b) => a + b, 0);
    for (let k = 0; k < this.ng; k++) Xg[k] /= sg;
    const lnG = this.lnGammaGroups(Xg, Psi);
    const lnGp = this.pure(T);
    return x.map((_, i) => {
      const J = this.r[i] / xr;
      const L = this.q[i] / xq;
      const lnC = 1 - J + Math.log(J) - 5 * this.q[i] * (1 - J / L + Math.log(J / L));
      let lnR = 0;
      for (let k = 0; k < this.ng; k++) lnR += this.nu[i][k] * (lnG[k] - lnGp[i][k]);
      return lnC + lnR;
    });
  }
}

const unifacCache = new Map<string, Unifac>();
function unifacFor(names: Compound[], oil: Oil) {
  const key = names.join("+") + "|" + oil;
  let u = unifacCache.get(key);
  if (!u) {
    u = new Unifac([...names.map((n) => M.compounds[n].groups), M.oils[oil].groups]);
    unifacCache.set(key, u);
  }
  return u;
}

export function lnAsat(c: Compound, T_K: number) {
  const d = M.compounds[c];
  let v = -(d.dHfus / R_GAS) * (1 / T_K - 1 / d.Tm);
  if (d.dCp) v += (d.dCp / R_GAS) * (d.Tm / T_K - 1 - Math.log(d.Tm / T_K));
  return v;
}

/* ------------------------------------------------------------------ mixture solver */
function solveLinear(Ain: number[][], bIn: number[]) {
  const n = bIn.length;
  const A = Ain.map((r) => r.slice());
  const b = bIn.slice();
  for (let c = 0; c < n; c++) {
    let p = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(A[r][c]) > Math.abs(A[p][c])) p = r;
    [A[c], A[p]] = [A[p], A[c]];
    [b[c], b[p]] = [b[p], b[c]];
    for (let r = c + 1; r < n; r++) {
      const f = A[r][c] / A[c][c];
      for (let k = c; k < n; k++) A[r][k] -= f * A[c][k];
      b[r] -= f * b[c];
    }
  }
  const x = new Array(n).fill(0);
  for (let r = n - 1; r >= 0; r--) {
    let s = b[r];
    for (let k = r + 1; k < n; k++) s -= A[r][k] * x[k];
    x[r] = s / A[r][r];
  }
  return x as number[];
}

const xFromN = (n: number[]) => {
  const tot = 1 + n.reduce((a, b) => a + b, 0);
  return [...n, 1].map((v) => v / tot);
};
const norm = (v: number[]) => Math.sqrt(v.reduce((s, a) => s + a * a, 0));
const maxAbs = (v: number[]) => v.reduce((s, a) => Math.max(s, Math.abs(a)), 0);

export type MixtureResult = {
  T_C: number;
  names: Compound[];
  x: number[];
  n: number[];
  gamma: number[];
  saturated: boolean[];
  mgPer100g: number[];
  iterations: number;
  residuals: number[];
  notes: string[];
};

/** Port of flavorsol.mixture.solve_mixture. loads: mg per 100 g oil (Infinity = excess) or null (all in excess). */
export function solveMixture(names: Compound[], oil: Oil, T_C: number, loads: number[] | null): MixtureResult {
  const T = T_C + 273.15;
  const k = names.length;
  const u = unifacFor(names, oil);
  const Moil = M.oils[oil].M;
  const lnas = names.map((c) => lnAsat(c, T));
  const notes: string[] = [];
  const history: number[] = [];
  const nTot = names.map((c, i) => {
    const L = loads ? loads[i] : Infinity;
    return isFinite(L) ? L / 1e3 / M.compounds[c].M / (100 / Moil) : Infinity;
  });

  let n = new Array(k).fill(0);
  const active = new Array(k).fill(false);
  names.forEach((c, i) => {
    const p = point(oil, c, T_C);
    if (T >= M.compounds[c].Tm) {
      if (!isFinite(nTot[i])) throw new Error(`${c} is molten at ${T_C.toFixed(1)} °C; give an amount to treat it as a dissolved liquid`);
      const nCap = p.status === "miscible" ? Infinity : p.x / (1 - p.x);
      n[i] = Math.min(nTot[i], nCap);
      notes.push(
        nTot[i] > nCap
          ? `${c} is molten and exceeds its liquid-liquid limit; the excess forms a separate liquid phase (binary estimate)`
          : `${c} above its melting point: dissolved as a liquid`,
      );
      return;
    }
    const nSingle = p.x / (1 - p.x);
    if (nSingle < nTot[i]) {
      n[i] = nSingle;
      active[i] = true;
    } else n[i] = nTot[i];
  });
  n = n.map((v) => Math.max(v, 1e-30));

  const residualsFor = (nn: number[], idx: number[]) => {
    const x = xFromN(nn);
    const lg = u.lnGamma(x, T);
    return idx.map((i) => Math.log(x[i]) + lg[i] - lnas[i]);
  };

  const newton = (n0: number[]) => {
    const idx = active.map((a, i) => (a ? i : -1)).filter((i) => i >= 0);
    const nn = n0.slice();
    if (idx.length === 0) return { n: nn, it: 0 };
    let uu = idx.map((i) => Math.log(nn[i]));
    const F_of = (uv: number[]) => {
      const t = nn.slice();
      idx.forEach((i, j) => (t[i] = Math.exp(uv[j])));
      return residualsFor(t, idx);
    };
    let F = F_of(uu);
    for (let it = 1; it <= 50; it++) {
      history.push(maxAbs(F));
      if (maxAbs(F) < 1e-10) {
        idx.forEach((i, j) => (nn[i] = Math.exp(uu[j])));
        return { n: nn, it: it - 1 };
      }
      const h = 1e-6;
      const J = idx.map(() => new Array(idx.length).fill(0));
      for (let j = 0; j < idx.length; j++) {
        const up = uu.slice();
        const dn = uu.slice();
        up[j] += h;
        dn[j] -= h;
        const Fp = F_of(up);
        const Fm = F_of(dn);
        for (let r = 0; r < idx.length; r++) J[r][j] = (Fp[r] - Fm[r]) / (2 * h);
      }
      let du = solveLinear(J, F.map((v) => -v));
      const step = maxAbs(du);
      if (step > 2) du = du.map((v) => (v * 2) / step);
      let lam = 1;
      const normF = norm(F);
      let Fnew = F;
      while (lam > 1e-6) {
        Fnew = F_of(uu.map((v, j) => v + lam * du[j]));
        if (norm(Fnew) < normF) break;
        lam *= 0.5;
      }
      uu = uu.map((v, j) => v + lam * du[j]);
      F = Fnew;
    }
    throw new Error(`Newton-Raphson did not converge at ${T_C} °C`);
  };

  let total = 0;
  for (let pass = 0; pass < 4 * k + 2; pass++) {
    const r = newton(n);
    total += r.it;
    n = r.n;
    const x = xFromN(n);
    const lg = u.lnGamma(x, T);
    let changed = false;
    for (let i = 0; i < k; i++) {
      if (active[i] && n[i] > nTot[i]) {
        active[i] = false;
        n[i] = nTot[i];
        changed = true;
      } else if (!active[i] && T < M.compounds[names[i]].Tm && Math.log(x[i]) + lg[i] - lnas[i] > 1e-9) {
        active[i] = true;
        changed = true;
      }
    }
    if (!changed) break;
  }
  const x = xFromN(n);
  const gamma = u.lnGamma(x, T).map(Math.exp);
  const mg = names.map((c, i) => (1e5 * x[i] * M.compounds[c].M) / (x[k] * Moil));
  return { T_C, names, x, n, gamma, saturated: active.slice(), mgPer100g: mg, iterations: total, residuals: history, notes };
}

/* ------------------------------------------------------------------ saturation statements */
export const fmtMg = (mg: number | null | undefined) => {
  if (mg == null || !isFinite(mg)) return "no limit";
  if (mg >= 1000) return `${(mg / 1000).toFixed(2)} g`;
  if (mg >= 100) return `${mg.toFixed(0)} mg`;
  if (mg >= 10) return `${mg.toFixed(1)} mg`;
  return `${mg.toFixed(2)} mg`;
};
export const fmtT = (t: number) => `${t.toFixed(1)} °C`;

export type SingleSat = {
  S: number | null;
  status: Status;
  capacityMg: number;
  dissolvedMg: number;
  excessMg: number;
  saturated: boolean;
  loading: number;
  Tsat: number | null;
  TsatKind: string;
  statement: string;
  headline: string;
};

/** Port of flavorsol.saturation.single_saturation. */
export function singleSaturation(c: Compound, oil: Oil, T_C: number, oilG: number, addedMg: number): SingleSat {
  const p = point(oil, c, T_C);
  const S = p.mgPer100g;
  const loading = oilG > 0 ? (addedMg * 100) / oilG : Infinity;
  const cap = S == null ? Infinity : (S * oilG) / 100;
  const saturated = isFinite(cap) && addedMg > 0 && addedMg >= cap * (1 - 1e-9);
  const dissolved = Math.min(addedMg, cap);
  const excess = Math.max(addedMg - dissolved, 0);
  const ts = saturationTemperature(oil, c, loading);
  const name = c;
  let headline: string;
  let text: string;
  if (p.status === "miscible") {
    headline = "No saturation: fully miscible";
    text = `At ${fmtT(T_C)} ${name} is above its melting point (${TmC(c).toFixed(0)} °C) and mixes with ${oil} oil in all proportions, so there is no saturation point: all ${fmtMg(addedMg)} mixes in.`;
  } else if (saturated) {
    headline = `At ${fmtT(T_C)} saturation is achieved`;
    const where = p.status === "SLE" ? "stays as undissolved solid" : "separates as a second liquid layer (molten compound)";
    text = `At ${fmtT(T_C)} saturation is achieved: ${oilG} g of ${oil} oil can hold only ${fmtMg(cap)} of ${name} (${fmtMg(S)} per 100 g oil). More mixing cannot be done; the extra ${fmtMg(excess)} ${where}.`;
    if (ts.kind === "root" && ts.T != null) text += ` To dissolve all ${fmtMg(addedMg)}, heat to at least ${fmtT(ts.T)}.`;
    else if (ts.kind === "above_range") text += ` Even at ${T_MAX} °C this amount does not fully dissolve.`;
  } else if (addedMg <= 0) {
    headline = "Nothing added yet";
    text = `At ${fmtT(T_C)} the oil would be saturated at ${fmtMg(cap)} of ${name} (${fmtMg(S)} per 100 g oil). Nothing has been added yet.`;
  } else {
    headline = `At ${fmtT(T_C)} saturation is not reached`;
    text = `At ${fmtT(T_C)} saturation is not reached: all ${fmtMg(addedMg)} of ${name} dissolves. The oil can take ${fmtMg(cap - addedMg)} more before it saturates (saturation point ${fmtMg(cap)} in ${oilG} g oil).`;
    if (ts.kind === "root" && ts.T != null) text += ` On cooling, this mixture becomes saturated at ${fmtT(ts.T)}.`;
    else if (ts.kind === "below_range") text += ` It stays unsaturated down to ${T_MIN} °C.`;
  }
  return { S, status: p.status, capacityMg: cap, dissolvedMg: dissolved, excessMg: excess, saturated, loading, Tsat: ts.T, TsatKind: ts.kind, statement: text, headline };
}

export type MixRow = {
  compound: Compound;
  addedMg: number;
  loading: number;
  dissolvedMg: number;
  excessMg: number;
  capacityMg: number;
  gamma: number;
  saturated: boolean;
  aloneMg: number | null;
  statement: string;
};
export type MixSat = {
  rows: MixRow[];
  headline: string;
  allDissolvedT: number | null;
  allDissolvedKind: string;
  iterations: number;
  notes: string[];
};

function mixFlags(names: Compound[], oil: Oil, T_C: number, loads: number[]) {
  const mix = solveMixture(names, oil, T_C, loads);
  const flags = names.map((_, i) => mix.saturated[i] || mix.mgPer100g[i] < loads[i] * (1 - 1e-7));
  return { mix, flags };
}

export function allDissolvedTemperature(names: Compound[], oil: Oil, loads: number[]) {
  const anySat = (t: number) => mixFlags(names, oil, t, loads).flags.some(Boolean);
  if (!anySat(T_MIN)) return { T: null as number | null, kind: "below_range" };
  if (anySat(T_MAX)) return { T: null as number | null, kind: "above_range" };
  let prev = T_MIN;
  for (let t = T_MIN + 5; t <= T_MAX; t += 5) {
    if (!anySat(t)) {
      let lo = prev;
      let hi = t;
      for (let j = 0; j < 25; j++) {
        const mid = 0.5 * (lo + hi);
        if (anySat(mid)) lo = mid;
        else hi = mid;
      }
      return { T: hi as number | null, kind: "root" };
    }
    prev = t;
  }
  return { T: null as number | null, kind: "above_range" };
}

/** Port of flavorsol.saturation.mixture_saturation. addedMg: mg of each compound; oilG: grams of oil. */
export function mixtureSaturation(oil: Oil, T_C: number, oilG: number, addedMg: Partial<Record<Compound, number>>, withTemperature = true): MixSat {
  const names = MIX_ORDER.filter((c) => (addedMg[c] ?? 0) > 0);
  if (names.length === 0) throw new Error("Add some of at least one compound.");
  const loads = names.map((c) => ((addedMg[c] as number) * 100) / oilG);
  const { mix, flags } = mixFlags(names, oil, T_C, loads);
  const rows: MixRow[] = names.map((c, i) => {
    const added = addedMg[c] as number;
    const dissolved = Math.min((mix.mgPer100g[i] * oilG) / 100, added);
    const excess = Math.max(added - dissolved, 0);
    const alone = point(oil, c, T_C).mgPer100g;
    let cap: number;
    let line: string;
    const Cap = c.charAt(0).toUpperCase() + c.slice(1);
    if (flags[i]) {
      cap = dissolved;
      const molten = T_C + 273.15 >= M.compounds[c].Tm;
      const where = molten ? "separates as a second liquid layer" : "stays as undissolved solid";
      line = `${Cap}: at ${fmtT(T_C)} saturation is achieved in the mixture at ${fmtMg(cap)} in ${oilG} g oil. More ${c} cannot be mixed in; the extra ${fmtMg(excess)} ${where}.`;
    } else {
      let note = "";
      if (T_C + 273.15 >= M.compounds[c].Tm) {
        cap = alone == null ? Infinity : (alone * oilG) / 100;
        note = " (molten; binary liquid-liquid estimate)";
      } else {
        const l2 = loads.slice();
        l2[i] = Infinity;
        cap = (solveMixture(names, oil, T_C, l2).mgPer100g[i] * oilG) / 100;
      }
      line = isFinite(cap)
        ? `${Cap}: at ${fmtT(T_C)} not saturated, all ${fmtMg(added)} dissolves; about ${fmtMg(Math.max(cap - added, 0))} more could mix in before saturation${note}.`
        : `${Cap}: at ${fmtT(T_C)} it is molten and mixes in all proportions (no saturation limit).`;
    }
    return { compound: c, addedMg: added, loading: loads[i], dissolvedMg: dissolved, excessMg: excess, capacityMg: cap, gamma: mix.gamma[i], saturated: flags[i], aloneMg: alone, statement: line };
  });
  const td = withTemperature ? allDissolvedTemperature(names, oil, loads) : { T: null, kind: "skipped" };
  const sat = rows.filter((r) => r.saturated).map((r) => r.compound);
  let head: string;
  if (sat.length) {
    const listed = sat.length === 1 ? sat[0] : sat.slice(0, -1).join(", ") + " and " + sat[sat.length - 1];
    head = `At ${fmtT(T_C)} saturation is achieved for ${listed} in ${oil} oil: more mixing cannot be done for ${sat.length === 1 ? "this compound" : "these compounds"}.`;
    if (td.kind === "root" && td.T != null) head += ` The whole recipe dissolves only at ${fmtT(td.T)} or above.`;
    else if (td.kind === "above_range") head += ` The recipe does not fully dissolve even at ${T_MAX} °C.`;
  } else {
    head = `At ${fmtT(T_C)} no compound is saturated: the whole recipe dissolves in ${oilG} g of ${oil} oil.`;
    if (td.kind === "below_range") head += ` It stays fully dissolved down to ${T_MIN} °C.`;
    else if (td.kind === "root" && td.T != null && td.T <= T_C) head += ` On cooling, saturation starts at ${fmtT(td.T)}.`;
  }
  return { rows, headline: head, allDissolvedT: td.T, allDissolvedKind: td.kind, iterations: mix.iterations, notes: mix.notes };
}

/** Share of each added compound that dissolves, 25-200 °C, for the recipe chart. */
export function recipeSweep(oil: Oil, oilG: number, addedMg: Partial<Record<Compound, number>>, step = 5) {
  const names = MIX_ORDER.filter((c) => (addedMg[c] ?? 0) > 0);
  const loads = names.map((c) => ((addedMg[c] as number) * 100) / oilG);
  const out: Record<string, number>[] = [];
  for (let t = T_MIN; t <= T_MAX; t += step) {
    const mix = solveMixture(names, oil, t, loads);
    const row: Record<string, number> = { T: t };
    names.forEach((c, i) => (row[c] = Math.min(100, (100 * mix.mgPer100g[i]) / loads[i])));
    out.push(row);
  }
  return out;
}

/* ------------------------------------------------------------------ boundary checks (flavorsol/checks.py) */
export type Note = { level: "error" | "warning" | "info"; message: string };
export function checks(c: Compound | null, oil: Oil, T_C: number, grade: "refined" | "unrefined"): Note[] {
  const out: Note[] = [];
  if (T_C > 152) out.push({ level: "info", message: `T = ${T_C.toFixed(0)} °C is outside the 2-152 °C window where most UNIFAC parameters were fitted; treat predictions as extrapolations.` });
  if (c && T_C + 273.15 >= M.compounds[c].Tm)
    out.push({ level: "info", message: `T is above the melting point of ${c} (${TmC(c).toFixed(0)} °C): no solid can exist; the result is a miscibility (liquid-liquid) limit, not a solid solubility.` });
  if (c === "curcumin" && T_C >= 150) out.push({ level: "warning", message: "Curcumin degrades in hot oil (38.9 % loss measured in corn oil heated to 200 °C). The dissolved amount falls with time; the prediction is an upper bound." });
  if (c === "piperine" && T_C >= 100) out.push({ level: T_C >= 180 ? "warning" : "info", message: "Piperine was stable as a solid at 100 °C for 4 days; its stability in oil above 100 °C is not verified." });
  if (c === "vanillin" && T_C >= 147) out.push({ level: "warning", message: "Above the flash point of vanillin (147 °C) it evaporates from an open pan; equilibrium solubility ignores this loss." });
  const smoke = { sunflower: { refined: 232, unrefined: 107 }, coconut: { refined: 204, unrefined: 177 } }[oil][grade];
  if (T_C >= smoke) out.push({ level: "warning", message: `T is at or above the smoke point of ${grade} ${oil} oil (about ${smoke} °C): the oil oxidises and decomposes, so the solvent changes.` });
  if (oil === "coconut" && T_C < 24) out.push({ level: "warning", message: "Coconut oil is partly solid below about 24 °C; the prediction assumes liquid oil." });
  return out;
}
