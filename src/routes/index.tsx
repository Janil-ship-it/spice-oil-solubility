import { createFileRoute } from "@tanstack/react-router";
import { useEffect, useMemo, useState } from "react";
import {
  CartesianGrid,
  Line,
  LineChart,
  ReferenceLine,
  ResponsiveContainer,
  Tooltip,
  XAxis,
  YAxis,
} from "recharts";
import {
  COMPOUNDS,
  OILS,
  batch,
  curve,
  flags,
  formatMg,
  type Compound,
  type Oil,
} from "@/lib/solubility";
import { api, type Issue, type Point } from "@/lib/api";
import { Field, IssueList, Segmented, Stat } from "@/components/ui-bits";
import { MixturePanel } from "@/components/mixture-panel";
import { ReportPanel } from "@/components/report-panel";

export const Route = createFileRoute("/")({
  component: Home,
});

type Tab = "single" | "mixture" | "report";

function Home() {
  const [tab, setTab] = useState<Tab>("single");
  return (
    <main className="mx-auto min-h-screen max-w-5xl px-4 py-8 sm:px-6">
      <header className="mb-6 max-w-2xl">
        <p className="text-sm font-medium tracking-wide text-accent uppercase">
          CL 202 · Group 4A · IIT Gandhinagar
        </p>
        <h1 className="mt-2 font-display text-4xl leading-tight text-ink sm:text-5xl">
          How much spice can hot oil hold?
        </h1>
        <p className="mt-3 text-lg text-body">
          UNIFAC solid–liquid equilibrium for curcumin, vanillin and piperine in cooking oils,
          computed live by the project’s Python model. One compound, all three together, or a full
          report with uncertainty.
        </p>
      </header>
      <div className="mb-6">
        <Segmented
          value={tab}
          onChange={setTab}
          options={[
            { id: "single", label: "Single compound" },
            { id: "mixture", label: "Mixture & recipe" },
            { id: "report", label: "Full report" },
          ]}
        />
      </div>
      {tab === "single" && <SinglePanel />}
      {tab === "mixture" && <MixturePanel />}
      {tab === "report" && <ReportPanel />}
    </main>
  );
}

function SinglePanel() {
  const [solute, setSolute] = useState<Compound>("curcumin");
  const [oil, setOil] = useState<Oil>("coconut");
  const [temp, setTemp] = useState(120);
  const [oilG, setOilG] = useState(15);
  const [addedMg, setAddedMg] = useState(47);
  const [grade, setGrade] = useState<"refined" | "unrefined">("refined");

  // Live values from the Python model; null = not loaded / API down (fallback to curves.json).
  const [online, setOnline] = useState<boolean | null>(null);
  const [liveCurve, setLiveCurve] = useState<{ key: string; pts: { T: number; mg: number | null }[] } | null>(null);
  const [livePoint, setLivePoint] = useState<{ key: string; p: Point } | null>(null);
  const [liveFlags, setLiveFlags] = useState<{ key: string; issues: Issue[] } | null>(null);

  const curveKey = `${solute}|${oil}`;
  const pointKey = `${solute}|${oil}|${temp}`;
  const flagKey = `${pointKey}|${grade}`;

  useEffect(() => {
    let live = true;
    api
      .curve(solute, oil)
      .then((r) => {
        if (!live) return;
        setOnline(true);
        setLiveCurve({ key: curveKey, pts: r.points.map((p) => ({ T: p.T_C, mg: p.mg_per_100g })) });
      })
      .catch(() => live && setOnline(false));
    return () => {
      live = false;
    };
  }, [curveKey, solute, oil]);

  useEffect(() => {
    if (online === false) return;
    let live = true;
    const t = setTimeout(() => {
      api
        .point(solute, oil, temp)
        .then((p) => live && setLivePoint({ key: pointKey, p }))
        .catch(() => live && setOnline(false));
      api
        .flags(solute, oil, temp, grade)
        .then((issues) => live && setLiveFlags({ key: flagKey, issues }))
        .catch(() => {});
    }, 150);
    return () => {
      live = false;
      clearTimeout(t);
    };
  }, [online, pointKey, flagKey, solute, oil, temp, grade]);

  const offline = useMemo(
    () => batch({ oil, solute, T: temp, oilG, addedMg }),
    [oil, solute, temp, oilG, addedMg],
  );
  const result = useMemo(() => {
    const lp = online && livePoint?.key === pointKey ? livePoint.p : null;
    if (!lp) return offline;
    const point = {
      ...offline.point,
      mgPer100g: lp.miscible ? null : lp.mg_per_100g,
      gamma: lp.gamma,
      x: lp.x,
      status: lp.status,
      miscible: lp.miscible,
    };
    const capacityMg = point.mgPer100g == null ? null : (point.mgPer100g * oilG) / 100;
    const dissolved = capacityMg == null ? addedMg : Math.min(addedMg, Math.max(0, capacityMg));
    return {
      point,
      capacityMg,
      dissolved,
      left: Math.max(0, addedMg - dissolved),
      share: addedMg > 0 ? dissolved / addedMg : 0,
    };
  }, [online, livePoint, pointKey, offline, oilG, addedMg]);

  const series = useMemo(
    () => (online && liveCurve?.key === curveKey ? liveCurve.pts : curve(oil, solute)),
    [online, liveCurve, curveKey, oil, solute],
  );
  const notes: { level: "error" | "warning" | "info"; message: string }[] =
    online && liveFlags?.key === flagKey
      ? liveFlags.issues
      : flags(solute, oil, temp, grade).map((f) => ({
          level: f.level === "stop" ? ("error" as const) : ("warning" as const),
          message: f.text,
        }));
  const spice = COMPOUNDS.find((c) => c.id === solute)!;
  const { point } = result;

  return (
    <>
      <p className="mb-4 text-sm">
        {online === null ? (
          <span className="rounded-full bg-wash px-3 py-1 text-muted">Connecting to the Python model…</span>
        ) : online ? (
          <span className="rounded-full bg-wash px-3 py-1 text-body">● Live: Python model (flavorsol)</span>
        ) : (
          <span className="rounded-full bg-gold px-3 py-1 text-ink">
            Offline: showing precomputed notebook values. Start the model with <code>npm run api</code>.
          </span>
        )}
      </p>
      <div className="grid gap-6 lg:grid-cols-[minmax(0,20rem)_1fr]">
        <form
          className="flex flex-col gap-5 rounded-2xl border border-line bg-surface p-5"
          onSubmit={(e) => e.preventDefault()}
        >
          <Field label="Compound">
            <div className="grid grid-cols-3 gap-2">
              {COMPOUNDS.map((c) => (
                <button
                  key={c.id}
                  type="button"
                  onClick={() => setSolute(c.id)}
                  className={
                    "min-h-11 rounded-xl border px-2 py-2 text-sm font-medium " +
                    (solute === c.id
                      ? "border-accent bg-accent text-bg"
                      : "border-line bg-wash text-ink")
                  }
                >
                  {c.label}
                </button>
              ))}
            </div>
            <p className="text-sm text-muted">{spice.note}</p>
          </Field>

          <Field label="Oil">
            <div className="flex flex-wrap gap-2">
              {OILS.map((o) => (
                <button
                  key={o.id}
                  type="button"
                  onClick={() => setOil(o.id)}
                  className={
                    "min-h-11 rounded-full border px-4 text-sm " +
                    (oil === o.id
                      ? "border-accent bg-gold text-ink"
                      : "border-line bg-wash text-body")
                  }
                >
                  {o.label}
                </button>
              ))}
            </div>
          </Field>

          <Field label={`Temperature · ${temp} °C`}>
            <input
              type="range"
              min={25}
              max={200}
              value={temp}
              onChange={(e) => setTemp(Number(e.target.value))}
              className="w-full"
            />
            <div className="flex justify-between text-xs text-muted">
              <span>25 °C</span>
              <span>200 °C</span>
            </div>
          </Field>

          <div className="grid grid-cols-2 gap-3">
            <Field label="Oil in the pan (g)">
              <input
                type="number"
                min={1}
                step={1}
                value={oilG}
                onChange={(e) => setOilG(Math.max(0, Number(e.target.value) || 0))}
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

          <Field label="Oil grade">
            <div className="flex gap-2">
              {(["refined", "unrefined"] as const).map((g) => (
                <button
                  key={g}
                  type="button"
                  onClick={() => setGrade(g)}
                  className={
                    "min-h-11 flex-1 rounded-xl border text-sm capitalize " +
                    (grade === g
                      ? "border-accent bg-wash text-ink"
                      : "border-line text-muted")
                  }
                >
                  {g}
                </button>
              ))}
            </div>
          </Field>

          <button
            type="button"
            className="min-h-11 rounded-xl bg-wash text-sm text-body"
            onClick={() => {
              setSolute("curcumin");
              setOil("coconut");
              setTemp(120);
              setOilG(15);
              setAddedMg(47);
              setGrade("refined");
            }}
          >
            Reset to a typical tadka (1.5 g turmeric → 47 mg curcumin)
          </button>
        </form>

        <section className="flex flex-col gap-4">
          <div className="rounded-2xl border border-line bg-surface p-5">
            <p className="text-sm text-muted">
              Saturation in {OILS.find((o) => o.id === oil)?.label.toLowerCase()} oil
            </p>
            <p className="mt-1 font-display text-4xl text-ink">
              {point.mgPer100g == null
                ? "Fully miscible"
                : `${point.mgPer100g >= 100 ? point.mgPer100g.toFixed(0) : point.mgPer100g.toFixed(2)} mg`}
              {point.mgPer100g != null && (
                <span className="ml-2 font-sans text-lg text-body">per 100 g oil</span>
              )}
            </p>
            <p className="mt-2 text-body">
              {point.miscible
                ? `${spice.label} has melted and mixes completely with the oil.`
                : `Mole fraction x = ${point.x.toExponential(2)} · activity coefficient γ = ${point.gamma.toFixed(1)} · ${point.status === "SLE" ? "solid still present at saturation" : point.status === "miscible" ? "one liquid phase" : "limited by a second liquid phase"}.`}
            </p>

            <dl className="mt-5 grid grid-cols-3 gap-3">
              <Stat k="Oil can hold" v={formatMg(result.capacityMg)} />
              <Stat k="Dissolves" v={formatMg(result.dissolved)} />
              <Stat k="Left solid" v={result.point.miscible ? "none" : formatMg(result.left)} />
            </dl>
            <div className="mt-4 h-2 overflow-hidden rounded-full bg-wash">
              <div
                className="h-full bg-accent"
                style={{ width: `${Math.round(result.share * 100)}%` }}
              />
            </div>
            <p className="mt-2 text-sm text-muted">
              {addedMg <= 0
                ? "Add some compound to see how much of it the oil takes."
                : `${Math.round(result.share * 100)}% of what you added dissolves at ${temp} °C.`}
            </p>
          </div>

          <IssueList items={notes} />

          <div className="rounded-2xl border border-line bg-surface p-4">
            <p className="mb-2 px-1 text-sm text-muted">
              Solubility vs temperature · marker at {temp} °C
            </p>
            <div className="h-64 w-full">
              <ResponsiveContainer width="100%" height="100%">
                <LineChart data={series} margin={{ top: 8, right: 8, left: 0, bottom: 0 }}>
                  <CartesianGrid stroke="#e6dccb" vertical={false} />
                  <XAxis
                    dataKey="T"
                    tick={{ fill: "#7a6e60", fontSize: 12 }}
                    stroke="#e6dccb"
                    unit="°"
                  />
                  <YAxis
                    tick={{ fill: "#7a6e60", fontSize: 12 }}
                    stroke="#e6dccb"
                    width={56}
                    scale="log"
                    domain={["auto", "auto"]}
                    allowDataOverflow
                  />
                  <Tooltip
                    formatter={(v) =>
                      v == null ? ["fully miscible", "solubility"] : [`${Number(v).toFixed(1)} mg/100 g`, "solubility"]
                    }
                    labelFormatter={(l) => `${l} °C`}
                  />
                  <ReferenceLine x={temp} stroke="#a24e14" strokeDasharray="4 4" />
                  <Line
                    type="monotone"
                    dataKey="mg"
                    stroke="#a24e14"
                    strokeWidth={2.5}
                    dot={false}
                    connectNulls={false}
                    name="mg / 100 g"
                  />
                </LineChart>
              </ResponsiveContainer>
            </div>
            <p className="px-1 pt-2 text-xs text-muted">
              {online
                ? "Computed live by the Python UNIFAC model (5 °C steps, 25–200 °C)."
                : "Precomputed notebook values (1 °C steps, 25–200 °C)."}{" "}
              A gap means the molten solute is fully miscible, so there is no solubility limit.
            </p>
          </div>
        </section>
      </div>
    </>
  );
}
