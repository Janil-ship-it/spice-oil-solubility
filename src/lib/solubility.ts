import curves from "@/data/curves.json";

export type Compound = "curcumin" | "vanillin" | "piperine";
export type Oil = "sunflower" | "coconut" | "corn" | "olive" | "MCT";

export const COMPOUNDS: { id: Compound; label: string; spice: string; note: string }[] = [
  { id: "curcumin", label: "Curcumin", spice: "Turmeric", note: "Yellow colour. Melts at 183 °C." },
  { id: "vanillin", label: "Vanillin", spice: "Vanilla", note: "Sweet aroma. Melts at 82 °C." },
  { id: "piperine", label: "Piperine", spice: "Black pepper", note: "Pungency. Melts at 133 °C." },
];

export const OILS: { id: Oil; label: string }[] = [
  { id: "sunflower", label: "Sunflower" },
  { id: "coconut", label: "Coconut" },
  { id: "corn", label: "Corn" },
  { id: "olive", label: "Olive" },
  { id: "MCT", label: "MCT" },
];

type Series = {
  mg: (number | null)[];
  gamma: number[];
  x: number[];
  runs: { from: number; status: string }[];
};

type OilBlock = { M: number } & Record<Compound, Series>;

const oils = curves.oils as Record<Oil, OilBlock>;
const meta = curves.compounds as Record<
  Compound,
  { M: number; Tm_C: number; dHfus_kJ: number; formula: string }
>;

const N = 176; // 25 °C … 200 °C

function at(series: Series, T: number) {
  const t = Math.min(200, Math.max(25, T));
  const i0 = Math.min(N - 1, Math.max(0, Math.floor(t) - 25));
  const i1 = Math.min(N - 1, i0 + 1);
  const f = t - Math.floor(t);
  const mix = (a: number, b: number) => a + (b - a) * f;
  const aMg = series.mg[i0];
  const bMg = series.mg[i1];
  const miscible = aMg == null || bMg == null || series.x[i0] >= 0.999;
  return {
    T: t,
    mgPer100g: miscible ? null : mix(aMg as number, bMg as number),
    gamma: mix(series.gamma[i0], series.gamma[i1]),
    x: mix(series.x[i0], series.x[i1]),
    status: statusAt(series, t),
    miscible,
  };
}

function statusAt(series: Series, T: number) {
  let status = series.runs[0]?.status ?? "SLE";
  for (const run of series.runs) {
    if (T + 1e-9 >= run.from) status = run.status;
  }
  return status;
}

export function lookup(oil: Oil, solute: Compound, T: number) {
  return { ...at(oils[oil][solute], T), formula: meta[solute].formula, Tm: meta[solute].Tm_C };
}

export function curve(oil: Oil, solute: Compound) {
  const series = oils[oil][solute];
  const pts: { T: number; mg: number | null }[] = [];
  for (let i = 0; i < N; i += 2) {
    pts.push({ T: 25 + i, mg: series.mg[i] });
  }
  if ((N - 1) % 2 !== 0) pts.push({ T: 200, mg: series.mg[N - 1] });
  return pts;
}

export type Flags = { level: "stop" | "warn"; text: string }[];

export function flags(solute: Compound, oil: Oil, T: number, grade: "refined" | "unrefined"): Flags {
  const out: Flags = [];
  if (T > 200 || T < 25) out.push({ level: "stop", text: "Temperature must be between 25 and 200 °C." });
  if (solute === "curcumin" && T >= 150)
    out.push({ level: "warn", text: "Above 150 °C curcumin starts to break down in hot oil. Equilibrium is no longer the whole story." });
  if (solute === "vanillin" && T >= 147)
    out.push({ level: "warn", text: "Above its flash point (147 °C) vanillin evaporates from an open pan." });
  if (solute === "piperine" && T > 100)
    out.push({ level: "warn", text: "Piperine stability above 100 °C is not verified. Treat this as an upper bound." });
  if (oil === "coconut" && T < 24)
    out.push({ level: "warn", text: "Coconut oil is partly solid below about 24 °C." });
  if (oil === "sunflower" && grade === "unrefined" && T >= 107)
    out.push({ level: "warn", text: "Unrefined sunflower oil smokes near 107 °C." });
  if (T > 152)
    out.push({ level: "warn", text: "Above 152 °C is outside the range UNIFAC was fitted for." });
  return out;
}

export function batch(opts: {
  oil: Oil;
  solute: Compound;
  T: number;
  oilG: number;
  addedMg: number;
}) {
  const point = lookup(opts.oil, opts.solute, opts.T);
  const capacityMg =
    point.mgPer100g == null ? null : (point.mgPer100g * opts.oilG) / 100;
  const dissolved =
    capacityMg == null ? opts.addedMg : Math.min(opts.addedMg, Math.max(0, capacityMg));
  const left = Math.max(0, opts.addedMg - dissolved);
  const share = opts.addedMg > 0 ? dissolved / opts.addedMg : 0;
  return { point, capacityMg, dissolved, left, share };
}

export function formatMg(n: number | null) {
  if (n == null) return "fully miscible";
  if (n >= 1000) return `${(n / 1000).toFixed(2)} g`;
  if (n >= 100) return `${n.toFixed(0)} mg`;
  if (n >= 10) return `${n.toFixed(1)} mg`;
  return `${n.toFixed(2)} mg`;
}
