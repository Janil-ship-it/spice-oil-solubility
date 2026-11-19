import type { ReactNode } from "react";


export function Field({ label, children }: { label: string; children: ReactNode }) {
  return (
    <label className="flex flex-col gap-2 text-sm font-medium text-ink">
      {label}
      {children}
    </label>
  );
}

export function Stat({ k, v }: { k: string; v: string }) {
  return (
    <div className="rounded-xl bg-wash px-3 py-3">
      <dt className="text-xs text-muted">{k}</dt>
      <dd className="mt-1 text-base font-medium text-ink">{v}</dd>
    </div>
  );
}

export function NumberInput({
  value,
  onChange,
  step = 1,
  min = 0,
}: {
  value: number;
  onChange: (v: number) => void;
  step?: number;
  min?: number;
}) {
  return (
    <input
      type="number"
      min={min}
      step={step}
      value={value}
      onChange={(e) => onChange(Math.max(0, Number(e.target.value) || 0))}
      className="min-h-11 w-full rounded-xl border border-line bg-bg px-3 text-ink"
    />
  );
}

export function Segmented<T extends string>({
  value,
  onChange,
  options,
}: {
  value: T;
  onChange: (v: T) => void;
  options: { id: T; label: string }[];
}) {
  return (
    <div className="inline-flex flex-wrap gap-1 rounded-xl border border-line bg-wash p-1">
      {options.map((o) => (
        <button
          key={o.id}
          type="button"
          onClick={() => onChange(o.id)}
          className={
            "min-h-10 rounded-lg px-4 text-sm font-medium " +
            (value === o.id ? "bg-accent text-bg" : "text-body")
          }
        >
          {o.label}
        </button>
      ))}
    </div>
  );
}

export function IssueList({ items }: { items: { level: "error" | "warning" | "info"; message: string }[] }) {
  if (items.length === 0) return null;
  return (
    <ul className="flex flex-col gap-2">
      {items.map((n, i) => (
        <li
          key={i + n.message}
          className={
            "rounded-xl border px-4 py-3 text-sm " +
            (n.level === "error" ? "border-accent bg-gold text-ink" : "border-line bg-wash text-body")
          }
        >
          {n.level !== "info" && <span className="mr-1 font-medium capitalize">{n.level}:</span>}
          {n.message}
        </li>
      ))}
    </ul>
  );
}

/** The saturation statement, shown prominently above every result. */
export function SaturationBanner({ saturated, title, text }: { saturated: boolean; title: string; text: string }) {
  return (
    <div
      role="status"
      className={
        "rounded-2xl border-2 px-5 py-4 " + (saturated ? "border-accent bg-gold/40 text-ink" : "border-line bg-surface text-ink")
      }
    >
      <p className="text-xs font-semibold tracking-wide uppercase text-accent">Saturation check</p>
      <p className="mt-1 font-display text-2xl">{title}</p>
      <p className="mt-2 text-body">{text}</p>
    </div>
  );
}
