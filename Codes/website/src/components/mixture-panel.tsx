import { useEffect, useMemo, useState } from "react";
import {
  Bar,
  BarChart,
  CartesianGrid,
  Legend,
  Line,
  LineChart,
  ReferenceLine,
  ResponsiveContainer,
  Tooltip,
  XAxis,
  YAxis,
} from "recharts";
import {
  LABEL,
  MIX_ORDER,
  OILS,
  T_MAX,
  T_MIN,
  checks,
  fmtMg,
  mixtureSaturation,
  recipeSweep,
  solveMixture,
  point,
  TmC,
  type Compound,
  type Oil,
} from "@/lib/engine";
import { Field, IssueList, NumberInput, SaturationBanner, Segmented } from "@/components/ui-bits";

const COLOR: Record<Compound, string> = { vanillin: "#6b4fa0", curcumin: "#d08a00", piperine: "#2f6f5e" };

/** Debounce so typing in number fields does not recompute on every key press. */
function useDebounced<T>(value: T, ms = 200) {
  const [v, setV] = useState(value);
  useEffect(() => {
    const t = setTimeout(() => setV(value), ms);
    return () => clearTimeout(t);
  }, [value, ms]);
  return v;
}

export function MixturePanel() {
  const [mode, setMode] = useState<"recipe" | "excess">("recipe");
  return (
    <div className="flex flex-col gap-6">
      <p className="max-w-3xl text-body">
        When several compounds dissolve in the same oil they share one liquid phase, so the activity coefficient of each
        one depends on how much of the others has dissolved. The coupled equilibrium equations are solved with damped
        Newton-Raphson and an active-set loop that decides which compounds saturate (solid left over) and which dissolve
        completely. For each compound the page states whether saturation is achieved at the chosen temperature.
      </p>
      <Segmented
        value={mode}
        onChange={setMode}
        options={[
          { id: "recipe", label: "Recipe in grams" },
          { id: "excess", label: "All compounds in excess" },
        ]}
      />
      {mode === "recipe" ? <RecipeMode /> : <ExcessMode />}
    </div>
  );
}

function OilPicker({ oil, setOil }: { oil: Oil; setOil: (o: Oil) => void }) {
  return (
    <Field label="Oil">
      <div className="grid grid-cols-2 gap-2">
        {OILS.map((o) => (
          <button
            key={o.id}
            type="button"
            onClick={() => setOil(o.id)}
            className={"min-h-11 rounded-full border px-4 text-sm " + (oil === o.id ? "border-accent bg-gold text-ink" : "border-line bg-wash text-body")}
          >
            {o.label}
          </button>
        ))}
      </div>
    </Field>
  );
}

function TempSlider({ temp, setTemp }: { temp: number; setTemp: (t: number) => void }) {
  return (
    <Field label={`Temperature · ${temp} °C`}>
      <input type="range" min={T_MIN} max={T_MAX} value={temp} onChange={(e) => setTemp(Number(e.target.value))} className="w-full" />
      <div className="flex justify-between text-xs font-normal text-muted">
        <span>{T_MIN} °C</span>
        <span>{T_MAX} °C</span>
      </div>
    </Field>
  );
}

/* ------------------------------------------------------------------ recipe (grams) */
function RecipeMode() {
  const [oil, setOil] = useState<Oil>("coconut");
  const [temp, setTemp] = useState(60);
  const [oilG, setOilG] = useState(15);
  const [mg, setMg] = useState<Record<Compound, number>>({ vanillin: 15, curcumin: 47, piperine: 50 });
  const q = useDebounced({ oil, temp, oilG, mg });

  const result = useMemo(() => {
    try {
      return { r: mixtureSaturation(q.oil, q.temp, q.oilG, q.mg), err: null as string | null };
    } catch (e) {
      return { r: null, err: (e as Error).message };
    }
  }, [q]);
  const sweep = useMemo(() => {
    try {
      return recipeSweep(q.oil, q.oilG, q.mg, 5);
    } catch {
      return [];
    }
  }, [q.oil, q.oilG, q.mg]);
  const r = result.r;
  const names = MIX_ORDER.filter((c) => mg[c] > 0);
  const bars =
    r?.rows.map((row) => ({ name: LABEL[row.compound], Dissolved: +row.dissolvedMg.toFixed(2), "Cannot mix in": +row.excessMg.toFixed(2) })) ?? [];

  return (
    <div className="grid gap-6 lg:grid-cols-[minmax(0,20rem)_1fr]">
      <form className="flex flex-col gap-5 rounded-2xl border border-line bg-surface p-5" onSubmit={(e) => e.preventDefault()}>
        <OilPicker oil={oil} setOil={setOil} />
        <TempSlider temp={temp} setTemp={setTemp} />
        <Field label="Oil in the pan (g)">
          <NumberInput value={oilG} step={1} min={0.1} onChange={(v) => setOilG(Math.max(0.1, v))} />
        </Field>
        <div className="grid grid-cols-3 gap-2">
          {MIX_ORDER.map((m) => (
            <Field key={m} label={LABEL[m]}>
              <NumberInput value={mg[m]} step={5} onChange={(v) => setMg({ ...mg, [m]: v })} />
            </Field>
          ))}
          <p className="col-span-3 text-xs text-muted">mg of pure compound added (0 = leave it out)</p>
        </div>
        <div className="flex flex-col gap-2">
          <button
            type="button"
            className="min-h-11 rounded-xl bg-wash text-sm text-body"
            onClick={() => {
              setOil("coconut");
              setTemp(60);
              setOilG(15);
              setMg({ vanillin: 15, curcumin: 47, piperine: 50 });
            }}
          >
            Tadka: 15 g oil, 1.5 g turmeric, 1 g pepper, vanilla
          </button>
          <button
            type="button"
            className="min-h-11 rounded-xl bg-wash text-sm text-body"
            onClick={() => {
              setOil("coconut");
              setTemp(60);
              setOilG(15);
              setMg({ vanillin: 0, curcumin: 47, piperine: 50 });
            }}
          >
            Turmeric + pepper only (notebook example)
          </button>
        </div>
      </form>

      <section className="flex flex-col gap-4">
        {result.err && <ErrorBox msg={result.err} />}
        {r && (
          <>
            <SaturationBanner saturated={r.rows.some((x) => x.saturated)} title={r.rows.some((x) => x.saturated) ? `At ${q.temp.toFixed(1)} °C saturation is achieved` : `At ${q.temp.toFixed(1)} °C saturation is not reached`} text={r.headline} />
            <ul className="flex flex-col gap-2">
              {r.rows.map((row) => (
                <li key={row.compound} className={"rounded-xl border px-4 py-3 text-sm " + (row.saturated ? "border-accent bg-gold/30 text-ink" : "border-line bg-surface text-body")}>
                  {row.statement}
                </li>
              ))}
            </ul>
            <div className="rounded-2xl border border-line bg-surface p-5">
              <div className="overflow-x-auto">
                <table className="w-full min-w-[36rem] text-left text-sm">
                  <thead className="text-xs text-muted">
                    <tr>
                      <th className="py-2 pr-3 font-medium">Compound</th>
                      <th className="py-2 pr-3 font-medium">Added</th>
                      <th className="py-2 pr-3 font-medium">Dissolved</th>
                      <th className="py-2 pr-3 font-medium">Cannot mix in</th>
                      <th className="py-2 pr-3 font-medium">Saturation point in mixture</th>
                      <th className="py-2 pr-3 font-medium">γ</th>
                      <th className="py-2 font-medium">State</th>
                    </tr>
                  </thead>
                  <tbody className="text-ink">
                    {r.rows.map((row) => (
                      <tr key={row.compound} className="border-t border-line">
                        <td className="py-2 pr-3 font-medium">{LABEL[row.compound]}</td>
                        <td className="py-2 pr-3">{fmtMg(row.addedMg)}</td>
                        <td className="py-2 pr-3">{fmtMg(row.dissolvedMg)}</td>
                        <td className="py-2 pr-3">{fmtMg(row.excessMg)}</td>
                        <td className="py-2 pr-3">{fmtMg(row.capacityMg)}</td>
                        <td className="py-2 pr-3">{row.gamma.toFixed(2)}</td>
                        <td className="py-2">
                          <span className={"rounded-full px-2 py-1 text-xs " + (row.saturated ? "bg-gold text-ink" : "bg-wash text-body")}>
                            {row.saturated ? "saturated" : "fully dissolved"}
                          </span>
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
              <div className="mt-4 h-52 w-full">
                <ResponsiveContainer width="100%" height="100%">
                  <BarChart data={bars} margin={{ top: 8, right: 8, left: 0, bottom: 0 }}>
                    <CartesianGrid stroke="#e6dccb" vertical={false} />
                    <XAxis dataKey="name" tick={{ fill: "#7a6e60", fontSize: 12 }} stroke="#e6dccb" />
                    <YAxis tick={{ fill: "#7a6e60", fontSize: 12 }} stroke="#e6dccb" width={56} unit=" mg" />
                    <Tooltip formatter={(v) => `${Number(v).toFixed(1)} mg`} />
                    <Legend />
                    <Bar dataKey="Dissolved" stackId="a" fill="#a24e14" isAnimationActive={false} />
                    <Bar dataKey="Cannot mix in" stackId="a" fill="#f0a04b" isAnimationActive={false} />
                  </BarChart>
                </ResponsiveContainer>
              </div>
              <p className="mt-2 text-xs text-muted">
                Dissolved = min(added, saturation point), solved jointly for all compounds in the pan. Newton-Raphson iterations: {r.iterations}. This is the
                thermodynamic upper bound; release from real spice tissue is slower.
              </p>
            </div>
            {sweep.length > 0 && (
              <div className="rounded-2xl border border-line bg-surface p-4">
                <p className="mb-2 px-1 text-sm text-muted">Share of each compound that dissolves vs temperature (same recipe)</p>
                <div className="h-60 w-full">
                  <ResponsiveContainer width="100%" height="100%">
                    <LineChart data={sweep} margin={{ top: 8, right: 12, left: 0, bottom: 0 }}>
                      <CartesianGrid stroke="#e6dccb" vertical={false} />
                      <XAxis dataKey="T" type="number" domain={[T_MIN, T_MAX]} ticks={[25, 50, 75, 100, 125, 150, 175, 200]} tick={{ fill: "#7a6e60", fontSize: 12 }} stroke="#e6dccb" unit="°" />
                      <YAxis domain={[0, 100]} tick={{ fill: "#7a6e60", fontSize: 12 }} stroke="#e6dccb" width={44} unit="%" />
                      <Tooltip formatter={(v) => `${Number(v).toFixed(1)} %`} labelFormatter={(l) => `${l} °C`} />
                      <Legend />
                      <ReferenceLine x={temp} stroke="#a24e14" strokeDasharray="4 4" />
                      {r.allDissolvedT != null && <ReferenceLine x={r.allDissolvedT} stroke="#1f1a17" strokeDasharray="2 3" label={{ value: `all dissolved ${r.allDissolvedT.toFixed(0)} °C`, fontSize: 11, fill: "#1f1a17", position: "insideBottomRight" }} />}
                      {names.map((c) => (
                        <Line key={c} type="monotone" dataKey={c} name={LABEL[c]} stroke={COLOR[c]} strokeWidth={2.5} dot={false} isAnimationActive={false} />
                      ))}
                    </LineChart>
                  </ResponsiveContainer>
                </div>
              </div>
            )}
            <IssueList items={[...r.notes.map((n) => ({ level: "info" as const, message: n })), ...checks(null, oil, temp, "refined")]} />
          </>
        )}
      </section>
    </div>
  );
}

/* ------------------------------------------------------------------ all in excess (co-solvency) */
function ExcessMode() {
  const [oil, setOil] = useState<Oil>("sunflower");
  const [temp, setTemp] = useState(60);
  const [use, setUse] = useState<Record<Compound, boolean>>({ vanillin: true, curcumin: true, piperine: true });
  const names = MIX_ORDER.filter((c) => use[c]);
  const res = useMemo(() => {
    if (names.length === 0) return { r: null, err: "Choose at least one compound." };
    const molten = names.filter((c) => temp >= TmC(c));
    if (molten.length)
      return { r: null, err: `${molten.map((c) => LABEL[c]).join(" and ")} ${molten.length === 1 ? "is" : "are"} molten at ${temp} °C, so ${molten.length === 1 ? "it" : "they"} cannot be present as excess solid. Lower the temperature or switch to Recipe in grams.` };
    try {
      return { r: solveMixture(names, oil, temp, null), err: null as string | null };
    } catch (e) {
      return { r: null, err: (e as Error).message };
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [oil, temp, names.join()]);
  const r = res.r;

  return (
    <div className="grid gap-6 lg:grid-cols-[minmax(0,20rem)_1fr]">
      <form className="flex flex-col gap-5 rounded-2xl border border-line bg-surface p-5" onSubmit={(e) => e.preventDefault()}>
        <OilPicker oil={oil} setOil={setOil} />
        <TempSlider temp={temp} setTemp={setTemp} />
        <Field label="Compounds present as excess solid">
          <div className="flex flex-col gap-2 font-normal">
            {MIX_ORDER.map((c) => (
              <label key={c} className="flex items-center gap-2 text-body">
                <input type="checkbox" checked={use[c]} onChange={(e) => setUse({ ...use, [c]: e.target.checked })} />
                {LABEL[c]}
              </label>
            ))}
          </div>
        </Field>
        <p className="text-xs text-muted">Every chosen compound is added in excess, so the oil is saturated with all of them at once.</p>
      </form>
      <section className="flex flex-col gap-4">
        {res.err && <ErrorBox msg={res.err} />}
        {r && (
          <>
            <SaturationBanner
              saturated
              title={`At ${temp} °C the oil is saturated with ${names.length === 1 ? LABEL[names[0]].toLowerCase() : "all " + names.length + " compounds"}`}
              text={`Every compound is present as excess solid, so at ${temp} °C saturation is achieved for each of them and more mixing cannot be done. The values below are the saturation points in the mixture; the column "alone" is the single-compound value at the same temperature.`}
            />
            <div className="overflow-x-auto rounded-2xl border border-line bg-surface p-5">
              <table className="w-full min-w-[34rem] text-left text-sm">
                <thead className="text-xs text-muted">
                  <tr>
                    <th className="py-2 pr-3 font-medium">Compound</th>
                    <th className="py-2 pr-3 font-medium">Saturation in mixture (mg/100 g)</th>
                    <th className="py-2 pr-3 font-medium">Alone (mg/100 g)</th>
                    <th className="py-2 pr-3 font-medium">Mixture / alone</th>
                    <th className="py-2 font-medium">γ</th>
                  </tr>
                </thead>
                <tbody className="text-ink">
                  {r.names.map((c, i) => {
                    const alone = point(oil, c, temp).mgPer100g;
                    return (
                      <tr key={c} className="border-t border-line">
                        <td className="py-2 pr-3 font-medium">{LABEL[c]}</td>
                        <td className="py-2 pr-3">{r.mgPer100g[i].toPrecision(4)}</td>
                        <td className="py-2 pr-3">{alone == null ? "miscible" : alone.toPrecision(4)}</td>
                        <td className="py-2 pr-3">{alone ? (r.mgPer100g[i] / alone).toFixed(2) : "-"}</td>
                        <td className="py-2">{r.gamma[i].toFixed(2)}</td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
              <p className="mt-3 text-xs text-muted">
                A ratio above 1 means the other dissolved compounds make the oil a better solvent (co-solvency). Newton-Raphson iterations: {r.iterations}.
              </p>
            </div>
          </>
        )}
      </section>
    </div>
  );
}

function ErrorBox({ msg }: { msg: string }) {
  return (
    <div className="rounded-2xl border border-accent bg-gold px-5 py-4 text-sm text-ink">
      <p className="font-medium">Cannot compute this case</p>
      <p className="mt-1">{msg}</p>
    </div>
  );
}
