import { useState } from "react";
import { api, type PredictOut } from "@/lib/api";
import { COMPOUNDS, OILS, type Compound, type Oil } from "@/lib/solubility";
import { Field, IssueList, Segmented } from "@/components/ui-bits";

const fmt = (n: number | null | undefined) =>
  n == null ? "—" : n >= 100 ? n.toFixed(0) : n >= 1 ? n.toFixed(2) : n.toPrecision(3);

/** Notebook §6: full workflow with routes R1–R5, checks and uncertainty (flavorsol.engine.predict). */
export function ReportPanel() {
  const [solute, setSolute] = useState<Compound>("piperine");
  const [oil, setOil] = useState<Oil>("coconut");
  const [temps, setTemps] = useState("40, 80, 120");
  const [grade, setGrade] = useState<"refined" | "unrefined">("refined");
  const [gamma, setGamma] = useState("");
  const [gammaT, setGammaT] = useState("60");
  const [meas, setMeas] = useState("");
  const [measT, setMeasT] = useState("25");
  const [mc, setMc] = useState(false);
  const [res, setRes] = useState<PredictOut | null>(null);
  const [err, setErr] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);

  const run = () => {
    const T_C = temps
      .split(/[,\s]+/)
      .map(Number)
      .filter((t) => Number.isFinite(t) && temps.trim() !== "");
    if (T_C.length === 0) {
      setErr("Enter at least one temperature.");
      return;
    }
    const g = parseFloat(gamma);
    const m = parseFloat(meas);
    setBusy(true);
    api
      .predict({
        solute,
        oil,
        T_C,
        oil_grade: grade,
        monte_carlo: mc ? 2000 : 0,
        ...(Number.isFinite(g) ? { gamma: g, gamma_T_C: parseFloat(gammaT) || 60 } : {}),
        ...(Number.isFinite(m)
          ? { measurements: [{ value: m, unit: "mg/100 g", T_C: parseFloat(measT) || 25 }] }
          : {}),
      })
      .then((r) => (setRes(r), setErr(null)))
      .catch((e: Error) => (setErr(e.message), setRes(null)))
      .finally(() => setBusy(false));
  };

  const input = "min-h-11 w-full rounded-xl border border-line bg-bg px-3 text-ink";

  return (
    <div className="grid gap-6 lg:grid-cols-[minmax(0,20rem)_1fr]">
      <form
        className="flex flex-col gap-5 rounded-2xl border border-line bg-surface p-5"
        onSubmit={(e) => {
          e.preventDefault();
          run();
        }}
      >
        <Field label="Compound">
          <Segmented value={solute} onChange={setSolute} options={COMPOUNDS.map((c) => ({ id: c.id, label: c.label }))} />
        </Field>
        <Field label="Oil">
          <select value={oil} onChange={(e) => setOil(e.target.value as Oil)} className={input}>
            {OILS.map((o) => (
              <option key={o.id} value={o.id}>
                {o.label}
              </option>
            ))}
          </select>
        </Field>
        <Field label="Temperatures (°C, comma separated)">
          <input value={temps} onChange={(e) => setTemps(e.target.value)} className={input} />
        </Field>
        <Field label="Oil grade">
          <Segmented
            value={grade}
            onChange={setGrade}
            options={[
              { id: "refined", label: "Refined" },
              { id: "unrefined", label: "Unrefined" },
            ]}
          />
        </Field>
        <details className="text-sm text-body">
          <summary className="cursor-pointer font-medium text-ink">Optional data (changes the route)</summary>
          <div className="mt-3 grid grid-cols-2 gap-3">
            <Field label="γ known (R1)">
              <input value={gamma} onChange={(e) => setGamma(e.target.value)} placeholder="e.g. 2.8" className={input} />
            </Field>
            <Field label="at T (°C)">
              <input value={gammaT} onChange={(e) => setGammaT(e.target.value)} className={input} />
            </Field>
            <Field label="Lab value mg/100 g (R2)">
              <input value={meas} onChange={(e) => setMeas(e.target.value)} placeholder="in this oil" className={input} />
            </Field>
            <Field label="at T (°C)">
              <input value={measT} onChange={(e) => setMeasT(e.target.value)} className={input} />
            </Field>
          </div>
          <label className="mt-3 flex items-center gap-2">
            <input type="checkbox" checked={mc} onChange={(e) => setMc(e.target.checked)} />
            Check uncertainty with Monte Carlo (2000 samples, slower)
          </label>
        </details>
        <button type="submit" disabled={busy} className="min-h-11 rounded-xl bg-accent text-sm font-medium text-bg disabled:opacity-60">
          {busy ? "Running…" : "Run prediction"}
        </button>
      </form>

      <section className="flex flex-col gap-4">
        {err && (
          <div className="rounded-2xl border border-accent bg-gold px-5 py-4 text-sm text-ink">
            <p className="font-medium">Couldn’t run the Python model</p>
            <p className="mt-1">{err}</p>
          </div>
        )}
        {!res && !err && (
          <p className="rounded-2xl border border-line bg-surface p-5 text-body">
            Runs <code>flavorsol.engine.predict</code>: input checks, route choice (R1 γ given, R2 anchored
            to a lab value, R5 pure UNIFAC), first-order uncertainty and breakdown warnings.
          </p>
        )}
        {res && res.predictions.length > 0 && (
          <div className="overflow-x-auto rounded-2xl border border-line bg-surface p-5">
            <table className="w-full min-w-[30rem] text-left text-sm">
              <thead className="text-xs text-muted">
                <tr>
                  <th className="py-2 pr-3 font-medium">T (°C)</th>
                  <th className="py-2 pr-3 font-medium">Route</th>
                  <th className="py-2 pr-3 font-medium">Status</th>
                  <th className="py-2 pr-3 font-medium">mg / 100 g</th>
                  <th className="py-2 pr-3 font-medium">95 % range</th>
                  <th className="py-2 font-medium">γ</th>
                </tr>
              </thead>
              <tbody className="text-ink">
                {res.predictions.map((p) => (
                  <tr key={p.T_C} className="border-t border-line">
                    <td className="py-2 pr-3">{p.T_C}</td>
                    <td className="py-2 pr-3">{p.route}</td>
                    <td className="py-2 pr-3">{p.status}</td>
                    <td className="py-2 pr-3">{fmt(p.mg_per_100g)}</td>
                    <td className="py-2 pr-3">
                      {p.ci95_mg_per_100g[0] == null ? "—" : `${fmt(p.ci95_mg_per_100g[0])}–${fmt(p.ci95_mg_per_100g[1])}`}
                    </td>
                    <td className="py-2">{p.gamma.toFixed(2)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        )}
        {res && <IssueList items={res.issues} />}
        {res && (
          <div className="rounded-2xl border border-line bg-surface p-4">
            <p className="mb-2 text-sm text-muted">Full report (same text as <code>python -m flavorsol predict</code>)</p>
            <pre className="overflow-x-auto text-xs leading-relaxed whitespace-pre text-ink">{res.report}</pre>
          </div>
        )}
      </section>
    </div>
  );
}
