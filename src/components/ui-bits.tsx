import type { ReactNode } from "react";
import type { Issue } from "@/lib/api";

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

export function IssueList({ items }: { items: Pick<Issue, "level" | "message">[] }) {
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
