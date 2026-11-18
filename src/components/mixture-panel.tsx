import { useEffect, useState } from "react";
import {
  Bar,
  BarChart,
  CartesianGrid,
  Legend,
  ResponsiveContainer,
  Tooltip,
  XAxis,
  YAxis,
} from "recharts";
import { api, type BatchOut, type MixtureOut } from "@/lib/api";
import { OILS, type Oil } from "@/lib/solubility";
import { Field, IssueList, NumberInput, Segmented } from "@/components/ui-bits";

const MIX = ["vanillin", "curcumin", "piperine"] as const;
type Mix = (typeof MIX)[number];
const LABEL: Record<Mix, string> = { vanillin: "Vanillin", curcumin: "Curcumin", piperine: "Piperine" };

const fmt = (n: number | null | undefined, d = 2) =>
  n == null ? "—" : Math.abs(n) >= 1000 ? n.toFixed(0) : n.toPrecision(Math.max(3, d + 1));

/** Debounced value so typing in number fields doesn't fire a request per keystroke. */
function useDebounced(value: string, ms = 300) {
  const [v, setV] = useState(value);
  useEffect(() => {
    const t = setTimeout(() => setV(value), ms);
    return () => clearTimeout(t);
  }, [value, ms]);
  return v;
}

export function MixturePanel() {
  const [mode, setMode] = useState<"loading" | "recipe">("loading");
  return (
    <div className="flex flex-col gap-6">
      <div className="max-w-2xl">
        <p className="text-body">
          All compounds share one oil phase, so each one’s activity coefficient depends on how much of
          the others has dissolved. The Python model solves the coupled equations with damped
          Newton–Raphson and an active-set loop (notebook §4 and §5).
        </p>
      </div>
      <Segmented
        value={mode}
        onChange={setMode}
        options={[
          { id: "loading", label: "Per 100 g oil (§4)" },
          { id: "recipe", label: "Recipe in grams (§5)" },
        ]}
      />
      {mode === "loading" ? <LoadingMode /> : <RecipeMode />}
    </div>
  );
}

function OilPicker({ oil, setOil }: { oil: Oil; setOil: (o: Oil) => void }) {
  return (
    <Field label="Oil">
      <div className="flex flex-wrap gap-2">
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
  );
}

function TempSlider({ temp, setTemp }: { temp: number; setTemp: (t: number) => void }) {
  return (
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
  );
}

/* ---------------------------------------------------------------- §4 */
function LoadingMode() {
  const [oil, setOil] = useState<Oil>("sunflower");
  const [temp, setTemp] = useState(60);
  const [excess, setExcess] = useState(false);
  const [loads, setLoads] = useState<Record<Mix, number>>({ vanillin: 100, curcumin: 314, piperine: 333 });
  const [res, setRes] = useState<MixtureOut | null>(null);
  const [err, setErr] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);

  const key = useDebounced(JSON.stringify({ oil, temp, excess, loads }));
  useEffect(() => {
    const q = JSON.parse(key) as { oil: Oil; temp: number; excess: boolean; loads: Record<Mix, number> };
    let live = true;
    setBusy(true);
    api
      .mixture({ oil: q.oil, T_C: q.temp, loadings: q.excess ? null : q.loads, compounds: [...MIX] })
      .then((r) => live && (setRes(r), setErr(null)))
      .catch((e: Error) => live && (setErr(e.message), setRes(null)))
      .finally(() => live && setBusy(false));
    return () => {
      live = false;
    };
  }, [key]);

  return (
    <div className="grid gap-6 lg:grid-cols-[minmax(0,20rem)_1fr]">
      <form className="flex flex-col gap-5 rounded-2xl border border-line bg-surface p-5" onSubmit={(e) => e.preventDefault()}>
        <OilPicker oil={oil} setOil={setOil} />
        <TempSlider temp={temp} setTemp={setTemp} />
        <Field label="Amounts added">
          <Segmented
            value={excess ? "excess" : "given"}
            onChange={(v) => setExcess(v === "excess")}
            options={[
              { id: "given", label: "Given loadings" },
              { id: "excess", label: "All in excess" },
            ]}
          />
        </Field>
        {!excess && (
          <div className="grid grid-cols-3 gap-2">
            {MIX.map((m) => (
              <Field key={m} label={`${LABEL[m]}`}>
                <NumberInput value={loads[m]} step={10} onChange={(v) => setLoads({ ...loads, [m]: v })} />
              </Field>
            ))}
            <p className="col-span-3 text-xs text-muted">mg added per 100 g oil (0 = leave it out)</p>
          </div>
        )}
        <button
          type="button"
          className="min-h-11 rounded-xl bg-wash text-sm text-body"
          onClick={() => {
            setOil("sunflower");
            setTemp(60);
            setExcess(false);
            setLoads({ vanillin: 100, curcumin: 314, piperine: 333 });
          }}
        >
          Reset to notebook example (sunflower, 60 °C, 100 / 314 / 333)
        </button>
      </form>

      <section className="flex flex-col gap-4">
        {err && <ErrorBox msg={err} />}
        {res && (
          <div className={"rounded-2xl border border-line bg-surface p-5 " + (busy ? "opacity-60" : "")}>
            <p className="text-sm text-muted">
              {OILS.find((o) => o.id === oil)?.label} oil at {res.T_C} °C · Newton–Raphson iterations:{" "}
              {res.newton_iterations}
            </p>
            <div className="mt-3 overflow-x-auto">
              <table className="w-full min-w-[34rem] text-left text-sm">
                <thead className="text-xs text-muted">
                  <tr>
                    <th className="py-2 pr-3 font-medium">Compound</th>
                    {!excess && <th className="py-2 pr-3 font-medium">Added</th>}
                    <th className="py-2 pr-3 font-medium">Dissolved</th>
                    <th className="py-2 pr-3 font-medium">Alone</th>
                    <th className="py-2 pr-3 font-medium">Mixture effect</th>
                    <th className="py-2 pr-3 font-medium">γ</th>
                    <th className="py-2 font-medium">State</th>
                  </tr>
                </thead>
                <tbody className="text-ink">
                  {res.rows.map((r) => (
                    <tr key={r.compound} className="border-t border-line">
                      <td className="py-2 pr-3 font-medium">{LABEL[r.compound as Mix] ?? r.compound}</td>
                      {!excess && <td className="py-2 pr-3">{fmt(r.loading_mg_per_100g)}</td>}
                      <td className="py-2 pr-3">{fmt(r.mg_per_100g)}</td>
                      <td className="py-2 pr-3 text-body">
                        {r.alone_mg_per_100g == null ? "miscible" : fmt(r.alone_mg_per_100g)}
                      </td>
                      <td className="py-2 pr-3">
                        {r.change_vs_alone_percent == null
                          ? "—"
                          : `${r.change_vs_alone_percent >= 0 ? "+" : ""}${r.change_vs_alone_percent.toFixed(1)} %`}
                      </td>
                      <td className="py-2 pr-3">{r.gamma.toFixed(2)}</td>
                      <td className="py-2">
                        <span
                          className={
                            "rounded-full px-2 py-1 text-xs " +
                            (r.saturated ? "bg-gold text-ink" : "bg-wash text-body")
                          }
                        >
                          {r.saturated ? "saturated, solid left" : "fully dissolved"}
                        </span>
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
            <p className="mt-3 text-xs text-muted">
              Values in mg per 100 g oil. “Alone” is the single-compound solubility at the same
              temperature; “mixture effect” shows how the other dissolved compounds change it.
            </p>
          </div>
        )}
        {res && res.notes.length > 0 && (
          <IssueList items={res.notes.map((n) => ({ level: "info" as const, code: "NOTE", message: n }))} />
        )}
        {res && <IssueList items={res.issues} />}
      </section>
    </div>
  );
}

/* ---------------------------------------------------------------- §5 */
function RecipeMode() {
  const [oil, setOil] = useState<Oil>("coconut");
  const [temp, setTemp] = useState(60);
  const [oilG, setOilG] = useState(15);
  const [grams, setGrams] = useState<Record<Mix, number>>({ vanillin: 0, curcumin: 0.047, piperine: 0.05 });
  const [res, setRes] = useState<BatchOut | null>(null);
  const [err, setErr] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);

  const key = useDebounced(JSON.stringify({ oil, temp, oilG, grams }));
  useEffect(() => {
    const q = JSON.parse(key) as { oil: Oil; temp: number; oilG: number; grams: Record<Mix, number> };
    let live = true;
    setBusy(true);
    api
      .batch({ oil: q.oil, T_C: q.temp, oil_g: q.oilG, added_g: q.grams })
      .then((r) => live && (setRes(r), setErr(null)))
      .catch((e: Error) => live && (setErr(e.message), setRes(null)))
      .finally(() => live && setBusy(false));
    return () => {
      live = false;
    };
  }, [key]);

  const chart =
    res?.rows.map((r) => ({
      name: LABEL[r.compound as Mix] ?? r.compound,
      Dissolved: +(r.dissolved_g * 1000).toFixed(3),
      "Left solid": +(r.undissolved_g * 1000).toFixed(3),
    })) ?? [];

  return (
    <div className="grid gap-6 lg:grid-cols-[minmax(0,20rem)_1fr]">
      <form className="flex flex-col gap-5 rounded-2xl border border-line bg-surface p-5" onSubmit={(e) => e.preventDefault()}>
        <OilPicker oil={oil} setOil={setOil} />
        <TempSlider temp={temp} setTemp={setTemp} />
        <Field label="Oil in the pan (g)">
          <NumberInput value={oilG} step={1} min={0.1} onChange={setOilG} />
        </Field>
        <div className="grid grid-cols-3 gap-2">
          {MIX.map((m) => (
            <Field key={m} label={LABEL[m]}>
              <NumberInput value={grams[m]} step={0.01} onChange={(v) => setGrams({ ...grams, [m]: v })} />
            </Field>
          ))}
          <p className="col-span-3 text-xs text-muted">grams of pure compound added (0 = leave it out)</p>
        </div>
        <button
          type="button"
          className="min-h-11 rounded-xl bg-wash text-sm text-body"
          onClick={() => {
            setOil("coconut");
            setTemp(60);
            setOilG(15);
            setGrams({ vanillin: 0, curcumin: 0.047, piperine: 0.05 });
          }}
        >
          Reset to notebook example (coconut, 60 °C, 15 g oil)
        </button>
      </form>

      <section className="flex flex-col gap-4">
        {err && <ErrorBox msg={err} />}
        {res && (
          <div className={"rounded-2xl border border-line bg-surface p-5 " + (busy ? "opacity-60" : "")}>
            <div className="overflow-x-auto">
              <table className="w-full min-w-[34rem] text-left text-sm">
                <thead className="text-xs text-muted">
                  <tr>
                    <th className="py-2 pr-3 font-medium">Compound</th>
                    <th className="py-2 pr-3 font-medium">Added</th>
                    <th className="py-2 pr-3 font-medium">Dissolved</th>
                    <th className="py-2 pr-3 font-medium">Left solid</th>
                    <th className="py-2 pr-3 font-medium">% dissolved</th>
                    <th className="py-2 pr-3 font-medium">ppm</th>
                    <th className="py-2 font-medium">Limited by</th>
                  </tr>
                </thead>
                <tbody className="text-ink">
                  {res.rows.map((r) => (
                    <tr key={r.compound} className="border-t border-line">
                      <td className="py-2 pr-3 font-medium">{LABEL[r.compound as Mix] ?? r.compound}</td>
                      <td className="py-2 pr-3">{(r.available_g * 1000).toFixed(1)} mg</td>
                      <td className="py-2 pr-3">{(r.dissolved_g * 1000).toFixed(1)} mg</td>
                      <td className="py-2 pr-3">{(r.undissolved_g * 1000).toFixed(1)} mg</td>
                      <td className="py-2 pr-3">{r.dissolved_percent == null ? "—" : `${r.dissolved_percent.toFixed(1)} %`}</td>
                      <td className="py-2 pr-3">{r.ppm_w == null ? "—" : r.ppm_w.toFixed(0)}</td>
                      <td className="py-2">{r.limited_by}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
            <div className="mt-4 h-56 w-full">
              <ResponsiveContainer width="100%" height="100%">
                <BarChart data={chart} margin={{ top: 8, right: 8, left: 0, bottom: 0 }}>
                  <CartesianGrid stroke="#e6dccb" vertical={false} />
                  <XAxis dataKey="name" tick={{ fill: "#7a6e60", fontSize: 12 }} stroke="#e6dccb" />
                  <YAxis tick={{ fill: "#7a6e60", fontSize: 12 }} stroke="#e6dccb" width={56} unit=" mg" />
                  <Tooltip formatter={(v) => `${Number(v).toFixed(1)} mg`} />
                  <Legend />
                  <Bar dataKey="Dissolved" stackId="a" fill="#a24e14" />
                  <Bar dataKey="Left solid" stackId="a" fill="#f0a04b" />
                </BarChart>
              </ResponsiveContainer>
            </div>
            <p className="mt-2 text-xs text-muted">
              Dissolved = min(added, capacity), solved jointly for all compounds in the pan. This is
              the thermodynamic upper bound; release from real spice tissue is slower.
            </p>
          </div>
        )}
      </section>
    </div>
  );
}

function ErrorBox({ msg }: { msg: string }) {
  return (
    <div className="rounded-2xl border border-accent bg-gold px-5 py-4 text-sm text-ink">
      <p className="font-medium">Couldn’t run the Python model</p>
      <p className="mt-1">{msg}</p>
    </div>
  );
}
