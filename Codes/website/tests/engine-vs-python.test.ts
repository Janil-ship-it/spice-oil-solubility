/*
 * Checks that the browser solver (src/lib/engine.ts) reproduces the Python
 * model (Codes/python/flavorsol). Reference numbers: tests/python_reference.json,
 * made by `python tools/make_web_reference.py`.   Run: npm test
 */
import { readFileSync } from "node:fs";
import { solveMixture, singleSaturation, mixtureSaturation, type Compound, type Oil } from "../src/lib/engine";

const ref = JSON.parse(readFileSync(new URL("./python_reference.json", import.meta.url), "utf8"));
let worst = 0;
let fails = 0;
const rel = (a: number, b: number) => Math.abs(a - b) / Math.max(Math.abs(b), 1e-12);
function check(label: string, ok: boolean, detail = "") {
  if (!ok) {
    fails++;
    console.log("FAIL", label, detail);
  }
}

for (const c of ref.mixture) {
  const r = solveMixture(c.names as Compound[], c.oil as Oil, c.T_C, c.loads);
  c.mg.forEach((v: number, i: number) => {
    const e = rel(r.mgPer100g[i], v);
    worst = Math.max(worst, e);
    check(`mixture ${c.oil} ${c.T_C} ${c.names[i]}`, e < 1e-3, `${r.mgPer100g[i]} vs ${v}`);
    check(`mixture flag ${c.oil} ${c.T_C} ${c.names[i]}`, r.saturated[i] === c.saturated[i]);
  });
}
console.log(`mixture cases: ${ref.mixture.length}, worst relative difference ${worst.toExponential(2)}`);

let worstT = 0;
for (const s of ref.single) {
  const r = singleSaturation(s.compound, s.oil, s.T_C, s.oil_g, s.added_mg);
  check(`single flag ${s.compound} ${s.oil} ${s.T_C}`, r.saturated === s.saturated);
  if (s.capacity_mg != null) check(`single capacity ${s.compound} ${s.oil} ${s.T_C}`, rel(r.capacityMg, s.capacity_mg) < 2e-3, `${r.capacityMg} vs ${s.capacity_mg}`);
  if (s.saturation_T_C != null && r.Tsat != null) {
    worstT = Math.max(worstT, Math.abs(r.Tsat - s.saturation_T_C));
    check(`single T_sat ${s.compound} ${s.oil}`, Math.abs(r.Tsat - s.saturation_T_C) < 0.1, `${r.Tsat} vs ${s.saturation_T_C}`);
  }
  check(`single T_sat kind ${s.compound} ${s.oil} ${s.T_C}`, r.TsatKind === s.saturation_T_kind, `${r.TsatKind} vs ${s.saturation_T_kind}`);
}
console.log(`single cases: ${ref.single.length}, worst saturation-temperature difference ${worstT.toFixed(3)} °C`);

for (const m of ref.mixture_saturation) {
  const r = mixtureSaturation(m.oil, m.T_C, 15, m.added_mg);
  m.rows.forEach((row: { compound: string; dissolved_mg: number; saturated: boolean }, i: number) => {
    check(`mixsat ${m.oil} ${m.T_C} ${row.compound}`, rel(r.rows[i].dissolvedMg, row.dissolved_mg) < 2e-3 && r.rows[i].saturated === row.saturated,
      `${r.rows[i].dissolvedMg} vs ${row.dissolved_mg}`);
  });
  check(`mixsat T ${m.oil} ${m.T_C}`, r.allDissolvedKind === m.kind && (m.all_dissolved_T_C == null || Math.abs((r.allDissolvedT ?? 0) - m.all_dissolved_T_C) < 0.3),
    `${r.allDissolvedT} vs ${m.all_dissolved_T_C}`);
}
console.log(`mixture saturation cases: ${ref.mixture_saturation.length}`);
if (fails) {
  console.log(`${fails} check(s) failed`);
  process.exit(1);
}
console.log("all checks passed: the website solver matches the Python model");
