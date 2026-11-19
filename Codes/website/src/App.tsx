import { useMemo, useState } from "react";
import {
  CartesianGrid,
  Line,
  LineChart,
  ReferenceDot,
  ReferenceLine,
  ResponsiveContainer,
  Tooltip,
  XAxis,
  YAxis,
} from "recharts";
import {
  COMPOUNDS,
  OILS,
  T_MAX,
  T_MIN,
  TmC,
  checks,
  curve,
  fmtMg,
  singleSaturation,
  type Compound,
  type Oil,
} from "@/lib/engine";
import { Field, IssueList, SaturationBanner, Segmented, Stat } from "@/components/ui-bits";
import { MixturePanel } from "@/components/mixture-panel";
import { MethodologyPanel } from "@/components/methodology-panel";
import { ReportPanel } from "@/components/report-panel";
import { REPO, TEAM } from "@/lib/meta";

type Tab = "single" | "mixture" | "method" | "report";


export default function App() {
  const [tab, setTab] = useState<Tab>("single");
  return (
    <main className="mx-auto min-h-screen max-w-5xl px-4 py-8 sm:px-6">
      <header className="mb-6 max-w-3xl">
        <p className="text-sm font-medium tracking-wide text-accent uppercase">
          CL 202 Chemical Engineering Thermodynamics · Group 4A · IIT Gandhinagar
        </p>
        <h1 className="mt-2 font-display text-4xl leading-tight text-ink sm:text-5xl">
          How much spice can hot oil hold?
        </h1>
        <p className="mt-3 text-lg text-body">
          Solid-liquid equilibrium with UNIFAC activity coefficients for curcumin (turmeric), vanillin (vanilla) and
          piperine (black pepper) in sunflower and coconut oil, 25 to 200 °C. For every case the page states
          whether saturation is reached at the chosen temperature, so you can see when no more compound can be mixed in.
        </p>
        <p className="mt-2 text-sm text-muted">{TEAM.join(", ")}</p>
      </header>
      <nav className="mb-6">
        <Segmented
          value={tab}
          onChange={setTab}
          options={[
            { id: "single", label: "Single compound" },
            { id: "mixture", label: "Mixtures" },
            { id: "method", label: "Methodology" },
            { id: "report", label: "Full report (Python)" },
          ]}
        />
      </nav>
      {tab === "single" && <SinglePanel />}
      {tab === "mixture" && <MixturePanel />}
      {tab === "method" && <MethodologyPanel />}
      {tab === "report" && <ReportPanel />}
      <footer className="mt-12 border-t border-line pt-4 text-sm text-muted">
        <p>
          Model: Codes/python/flavorsol (Python). This page runs a verified copy of the same equations in the browser,
          so it works offline. Source and data: <a className="text-accent underline" href={REPO}>{REPO.replace("https://", "")}</a>
        </p>
      </footer>
    </main>
  );
}

function SinglePanel() {
  const [solute, setSolute] = useState<Compound>("curcumin");
  const [oil, setOil] = useState<Oil>("coconut");
  const [temp, setTemp] = useState(60);
  const [oilG, setOilG] = useState(15);
  const [addedMg, setAddedMg] = useState(47);
  const [grade, setGrade] = useState<"refined" | "unrefined">("refined");

  const sat = useMemo(() => singleSaturation(solute, oil, temp, oilG, addedMg), [solute, oil, temp, oilG, addedMg]);
  const series = useMemo(() => curve(oil, solute, 1), [oil, solute]);
  const notes = checks(solute, oil, temp, grade);
  const spice = COMPOUNDS.find((c) => c.id === solute)!;
  const share = addedMg > 0 ? sat.dissolvedMg / addedMg : 0;
  const tm = TmC(solute);

  return (
    <div className="grid gap-6 lg:grid-cols-[minmax(0,20rem)_1fr]">
      <form className="flex flex-col gap-5 rounded-2xl border border-line bg-surface p-5" onSubmit={(e) => e.preventDefault()}>
        <Field label="Compound">
          <div className="grid grid-cols-3 gap-2">
            {COMPOUNDS.map((c) => (
              <button
                key={c.id}
                type="button"
                onClick={() => setSolute(c.id)}
                className={
                  "min-h-11 rounded-xl border px-2 py-2 text-sm font-medium " +
                  (solute === c.id ? "border-accent bg-accent text-bg" : "border-line bg-wash text-ink")
                }
              >
                {c.label}
              </button>
            ))}
          </div>
          <p className="text-sm font-normal text-muted">{spice.note}</p>
        </Field>

        <Field label="Oil">
          <div className="grid grid-cols-2 gap-2">
            {OILS.map((o) => (
              <button
                key={o.id}
                type="button"
                onClick={() => setOil(o.id)}
                className={
                  "min-h-11 rounded-full border px-4 text-sm " +
                  (oil === o.id ? "border-accent bg-gold text-ink" : "border-line bg-wash text-body")
                }
              >
                {o.label}
              </button>
            ))}
          </div>
        </Field>

        <Field label={`Temperature · ${temp} °C`}>
          <input type="range" min={T_MIN} max={T_MAX} value={temp} onChange={(e) => setTemp(Number(e.target.value))} className="w-full" />
          <div className="flex justify-between text-xs font-normal text-muted">
            <span>{T_MIN} °C</span>
            <span>{T_MAX} °C</span>
          </div>
        </Field>

        <div className="grid grid-cols-2 gap-3">
          <Field label="Oil in the pan (g)">
            <input
              type="number"
              min={1}
              step={1}
              value={oilG}
              onChange={(e) => setOilG(Math.max(0.1, Number(e.target.value) || 0.1))}
              className="min-h-11 w-full rounded-xl border border-line bg-bg px-3 text-ink"
            />
          </Field>
          <Field label={`${spice.label} added (mg)`}>
            <input
              type="number"
              min={0}
              step={1}
              value={addedMg}
              onChange={(e) => setAddedMg(Math.max(0, Number(e.target.value) || 0))}
              className="min-h-11 w-full rounded-xl border border-line bg-bg px-3 text-ink"
            />
          </Field>
        </div>

        <Field label="Oil grade (smoke-point check)">
          <Segmented
            value={grade}
            onChange={setGrade}
            options={[
              { id: "refined", label: "Refined" },
              { id: "unrefined", label: "Unrefined" },
            ]}
          />
        </Field>

        <button
          type="button"
          className="min-h-11 rounded-xl bg-wash text-sm text-body"
          onClick={() => {
            setSolute("curcumin");
            setOil("coconut");
            setTemp(60);
            setOilG(15);
            setAddedMg(47);
            setGrade("refined");
          }}
        >
          Reset to a typical tadka (1.5 g turmeric gives 47 mg curcumin in 15 g oil)
        </button>
      </form>

      <section className="flex flex-col gap-4">
        <SaturationBanner saturated={sat.saturated} title={sat.headline} text={sat.statement} />

        <div className="rounded-2xl border border-line bg-surface p-5">
          <p className="text-sm text-muted">Saturation solubility in {oil} oil at {temp} °C</p>
          <p className="mt-1 font-display text-4xl text-ink">
            {sat.S == null ? "Fully miscible" : fmtMg(sat.S)}
            {sat.S != null && <span className="ml-2 font-sans text-lg text-body">per 100 g oil</span>}
          </p>
          <p className="mt-2 text-sm text-body">
            {sat.status === "miscible"
              ? `${spice.label} has melted (${tm.toFixed(0)} °C) and mixes completely with the oil.`
              : sat.status === "SLE"
                ? "Limit set by the solid: excess compound stays as crystals (solid-liquid equilibrium)."
                : "Above the melting point: limit set by a second liquid layer of molten compound (liquid-liquid equilibrium)."}
          </p>
          <dl className="mt-5 grid grid-cols-2 gap-3 sm:grid-cols-4">
            <Stat k="Saturation point (this pan)" v={fmtMg(sat.capacityMg)} />
            <Stat k="Dissolves" v={fmtMg(sat.dissolvedMg)} />
            <Stat k="Cannot mix in" v={sat.status === "miscible" ? "none" : fmtMg(sat.excessMg)} />
            <Stat
              k="Saturation temperature"
              v={sat.Tsat != null ? `${sat.Tsat.toFixed(1)} °C` : sat.TsatKind === "below_range" ? `below ${T_MIN} °C` : sat.TsatKind === "above_range" ? `above ${T_MAX} °C` : "-"}
            />
          </dl>
          <div className="mt-4 h-2 overflow-hidden rounded-full bg-wash">
            <div className="h-full bg-accent" style={{ width: `${Math.round(share * 100)}%` }} />
          </div>
          <p className="mt-2 text-sm text-muted">
            {addedMg <= 0 ? "Add some compound to see how much the oil takes." : `${Math.round(share * 100)} % of what you added dissolves at ${temp} °C.`}
          </p>
        </div>

        <IssueList items={notes} />

        <div className="rounded-2xl border border-line bg-surface p-4">
          <p className="mb-2 px-1 text-sm text-muted">
            Solubility (mg per 100 g oil) vs temperature. Dashed line: your loading ({fmtMg(sat.loading)} per 100 g). Dot: where the curve meets it, the saturation temperature.
          </p>
          <div className="h-72 w-full">
            <ResponsiveContainer width="100%" height="100%">
              <LineChart data={series} margin={{ top: 8, right: 12, left: 0, bottom: 0 }}>
                <CartesianGrid stroke="#e6dccb" vertical={false} />
                <XAxis dataKey="T" type="number" domain={[T_MIN, T_MAX]} ticks={[25, 50, 75, 100, 125, 150, 175, 200]} tick={{ fill: "#7a6e60", fontSize: 12 }} stroke="#e6dccb" unit="°" />
                <YAxis tick={{ fill: "#7a6e60", fontSize: 12 }} stroke="#e6dccb" width={60} scale="log" domain={[1, 1e6]} allowDataOverflow tickFormatter={(v) => (v >= 1000 ? `${v / 1000}k` : `${v}`)} />
                <Tooltip
                  formatter={(v) => (v == null ? ["fully miscible", "solubility"] : [`${Number(v).toFixed(1)} mg/100 g`, "solubility"])}
                  labelFormatter={(l) => `${l} °C`}
                />
                <ReferenceLine x={temp} stroke="#a24e14" strokeDasharray="4 4" label={{ value: `${temp} °C`, fill: "#a24e14", fontSize: 11, position: "insideTopRight" }} />
                {sat.loading > 0 && isFinite(sat.loading) && (
                  <ReferenceLine y={sat.loading} stroke="#4a4038" strokeDasharray="6 3" />
                )}
                {tm <= T_MAX && <ReferenceLine x={tm} stroke="#7a6e60" strokeDasharray="2 4" label={{ value: "Tm", fill: "#7a6e60", fontSize: 11, position: "insideBottomLeft" }} />}
                {sat.Tsat != null && sat.loading > 0 && <ReferenceDot x={sat.Tsat} y={sat.loading} r={5} fill="#a24e14" stroke="#fffdf8" />}
                <Line type="monotone" dataKey="mg" stroke="#a24e14" strokeWidth={2.5} dot={false} connectNulls={false} isAnimationActive={false} />
              </LineChart>
            </ResponsiveContainer>
          </div>
          <p className="px-1 pt-2 text-xs text-muted">
            Values from the Python UNIFAC model in 1 °C steps. Below the dot the oil is saturated (solid left over); above it everything dissolves.
            A gap in the curve means the molten compound is fully miscible, so there is no saturation limit.
          </p>
        </div>
      </section>
    </div>
  );
}
