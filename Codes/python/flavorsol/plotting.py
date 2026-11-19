"""Figures for the report and slides (Member 5 - Analysis)."""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

# Validated categorical slots (first three pass all-pairs CVD checks)
COLORS = {"vanillin": "#2a78d6", "curcumin": "#eb6834", "piperine": "#1baf7a"}
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8985", "#e4e3de"
OIL_STYLE = {"sunflower": "-", "coconut": (0, (5, 2.2))}

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 8.5,
    "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "legend.fontsize": 7, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "lines.linewidth": 1.5, "legend.frameon": False,
    "savefig.dpi": 300, "figure.dpi": 150, "axes.titleweight": "bold",
    "axes.titlelocation": "left",
})


def save(fig, outdir, name):
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(os.path.join(outdir, name + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(outdir, name + ".png"), bbox_inches="tight")
    plt.close(fig)


def fig_validation(points, oil_rows, outdir):
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 3.1))
    # (a) analogues in pure solvents
    alcohols = {"ethanol", "1-propanol", "2-propanol"}
    xe = np.array([p["x_exp"] for p in points])
    a.plot([3e-3, 1], [3e-3, 1], color=INK2, lw=1)
    for k in (2, 0.5):
        a.plot([3e-3, 1], [3e-3 * k, k], color=MUTED, lw=0.8, ls=(0, (3, 2)))
    a.scatter(xe, [p["x_ideal"] for p in points], s=14, facecolors="none",
              edgecolors=MUTED, lw=0.8, label="ideal solution (γ = 1)", zorder=2)
    for cls, col, mk in (("aprotic (esters, ketone)", COLORS["vanillin"], "o"),
                         ("alcohols", COLORS["curcumin"], "s")):
        sel = [p for p in points if (p["solvent"] in alcohols) == (cls == "alcohols")]
        a.scatter([p["x_exp"] for p in sel], [p["x_calc"] for p in sel], s=16, marker=mk,
                  color=col, edgecolors="white", lw=0.4, label=f"UNIFAC, {cls}", zorder=3)
    a.set(xscale="log", yscale="log", xlim=(7e-3, 0.6), ylim=(7e-3, 0.9),
          xlabel="measured solubility, $x_s$", ylabel="predicted solubility, $x_s$",
          title="(a) Vanillin analogues in pure solvents")
    a.text(0.30, 0.62, "×2", color=MUTED, fontsize=7)
    a.text(0.30, 0.13, "×½", color=MUTED, fontsize=7)
    a.legend(loc="upper left", handletextpad=0.3)

    # (b) published solubility in edible oils
    lo, hi = 0.01, 30
    b.plot([lo, hi], [lo, hi], color=INK2, lw=1)
    for k in (3, 1 / 3):
        b.plot([lo, hi], [lo * k, hi * k], color=MUTED, lw=0.8, ls=(0, (3, 2)))
    b.text(0.75, 3.9, "×3", color=MUTED, fontsize=7)
    b.text(2.2, 0.48, "×⅓", color=MUTED, fontsize=7)
    offsets = {("curcumin", "MCT", 25): (5, -2, "left"), ("curcumin", "coconut", 25): (5, -1, "left"),
               ("curcumin", "olive", 25): (5, -7, "left"), ("curcumin", "corn", 25): (-5, 1, "right"),
               ("curcumin", "corn", 100): (5, -3, "left"), ("piperine", "MCT", 37): (-6, 4, "right"),
               ("curcumin", "MCT", 100): (5, -10, "left")}
    for r in oil_rows:
        col = COLORS[r["solute"]]
        tc = round(r["T"] - 273.15)
        lab = f"{r['oil']} {tc} °C"
        dx, dy, ha = offsets.get((r["solute"], r["oil"], tc), (5, -3, "left"))
        if r["kind"] == "oil":
            b.scatter(r["exp_mg_mL"], r["pred_mg_mL"], s=26, color=col, edgecolors="white",
                      lw=0.5, zorder=3)
            b.annotate(lab, (r["exp_mg_mL"], r["pred_mg_mL"]), xytext=(dx, dy), ha=ha,
                       textcoords="offset points", fontsize=6.3, color=INK2)
        elif r["kind"] == "oil-lower-bound":
            b.scatter(r["exp_mg_mL"], r["pred_mg_mL"], s=30, color=col, marker=">", zorder=3)
            b.annotate(lab + " (≥)", (r["exp_mg_mL"], r["pred_mg_mL"]), xytext=(dx, dy), ha=ha,
                       textcoords="offset points", fontsize=6.3, color=INK2)
        else:
            b.scatter(r["exp_mg_mL"], r["pred_mg_mL"], s=24, facecolors="none",
                      edgecolors=col, lw=1, marker="D", zorder=3)
        if r["solute"] == "curcumin" and r["kind"] == "oil":
            b.scatter(r["exp_mg_mL"], r["ideal_mg_mL"], s=12, facecolors="none",
                      edgecolors=MUTED, lw=0.8, zorder=2)
    b.annotate("reports above the\nideal-solution limit", (8.0, 0.045), xytext=(0, 14),
               textcoords="offset points", fontsize=6.3, color=INK2, ha="center")
    b.scatter([], [], s=26, color=COLORS["curcumin"], label="curcumin")
    b.scatter([], [], s=26, color=COLORS["piperine"], label="piperine")
    b.scatter([], [], s=24, facecolors="none", edgecolors=COLORS["curcumin"], marker="D",
              label="curcumin, high-value reports")
    b.scatter([], [], s=12, facecolors="none", edgecolors=MUTED, label="ideal solution (curcumin)")
    b.set(xscale="log", yscale="log", xlim=(lo, hi), ylim=(lo, 200),
          xlabel="published solubility (mg per mL oil)", ylabel="predicted (mg per mL oil)",
          title="(b) Flavour compounds in edible oils")
    b.legend(loc="upper left", handletextpad=0.3)
    fig.tight_layout(w_pad=2.0)
    save(fig, outdir, "fig2_validation")


def fig_solubility_T(rows, tm_c, lit_points, outdir):
    solutes = ["vanillin", "curcumin", "piperine"]
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.75), sharey=True)
    for ax, s in zip(axes, solutes):
        col = COLORS[s]
        for oil, ls in OIL_STYLE.items():
            sel = [r for r in rows if r["solute"] == s and r["oil"] == oil]
            T = np.array([r["T_C"] for r in sel])
            y = np.array([r["mg_per_100g"] for r in sel])
            st = np.array([r["status"] for r in sel])
            sle = st == "SLE"
            ax.plot(T[sle], y[sle], color=col, ls=ls, label=f"{oil} oil")
            lle = st == "LLE-limited"
            if lle.any():
                ax.plot(T[lle], y[lle], color=col, ls=ls, lw=1.0, alpha=0.6)
            if oil == "sunflower":
                yi = np.array([r["mg_per_100g_ideal"] for r in sel])
                ok = np.isfinite(yi) & (T < tm_c[s])
                ax.plot(T[ok], yi[ok], color=MUTED, ls=(0, (1, 1.5)), lw=1.2,
                        label="ideal (sunflower)")
        ax.axvline(tm_c[s], color=INK2, lw=0.8, ls=(0, (4, 2)))
        ax.text(tm_c[s] - 3, 1.3, f"$T_m$ = {tm_c[s]:.0f} °C", rotation=90, fontsize=6.5,
                color=INK2, ha="right", va="bottom")
        for lp in lit_points.get(s, []):
            ax.scatter(lp["T_C"], lp["mg_per_100g"], s=22, color=col, edgecolors=INK,
                       lw=0.6, marker="o", zorder=4)
            ax.annotate(lp["label"], (lp["T_C"], lp["mg_per_100g"]), xytext=(4, -9),
                        textcoords="offset points", fontsize=6.2, color=INK2)
        if s == "vanillin":
            ax.text(132, 1.1e4, "molten: liquid–liquid\nlimit", fontsize=6.3, color=INK2,
                    ha="center")
        if s == "piperine":
            ax.annotate("fully miscible\nabove $T_m$", (tm_c[s] + 3, 2.2e4), fontsize=6.3,
                        color=INK2)
        ax.set(yscale="log", xlim=(20, 185), ylim=(1, 1e5), xlabel="temperature (°C)",
               title=s.capitalize())
        ax.set_xticks([25, 60, 100, 140, 180])
    axes[0].set_ylabel("solubility (mg per 100 g oil)")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=INK2, ls="-"), Line2D([], [], color=INK2, ls=OIL_STYLE["coconut"]),
         Line2D([], [], color=MUTED, ls=(0, (1, 1.5)), lw=1.2)]
    l = ["sunflower oil", "coconut oil", "ideal solution (sunflower)"]
    fig.legend(h, l, loc="lower center", ncol=3, handlelength=2.6, bbox_to_anchor=(0.5, -0.04))
    fig.tight_layout(w_pad=0.8, rect=(0, 0.07, 1, 1))
    save(fig, outdir, "fig3_solubility_vs_T")


def fig_gamma(decomp, rows, outdir):
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.7), gridspec_kw={"width_ratios": [1, 1.25]})
    labels, C, R = [], [], []
    for d in decomp:
        labels.append(f"{d['solute']}\n{d['oil']}")
        C.append(d["lnC"])
        R.append(d["lnR"])
    xpos = np.arange(len(labels))
    a.bar(xpos, R, width=0.62, color=INK2, label="residual (energy)", zorder=2)
    a.bar(xpos, C, width=0.62, color="#b4b2a9", label="combinatorial (size/shape)", zorder=2)
    a.scatter(xpos, np.array(C) + np.array(R), color=INK, s=14, zorder=3, label="ln γ (total)")
    a.axhline(0, color=INK2, lw=0.8)
    a.set_xticks(xpos, ["sunfl.", "coco."] * 3, fontsize=6.5)
    for k, sname in enumerate(("vanillin", "curcumin", "piperine")):
        a.text(2 * k + 0.5, -0.2, sname, transform=a.get_xaxis_transform(), ha="center",
               fontsize=7, color=INK)
    a.set_ylim(-1.5, 6.3)
    a.set(ylabel="ln γ$_s^\\infty$ at 25 °C", title="(a) Where the non-ideality comes from")
    a.legend(loc="upper right", fontsize=6.5)
    a.grid(axis="x", visible=False)
    for s in ("vanillin", "curcumin", "piperine"):
        for oil, ls in OIL_STYLE.items():
            sel = [r for r in rows if r["solute"] == s and r["oil"] == oil and r["status"] == "SLE"]
            b.plot([r["T_C"] for r in sel], [r["gamma"] for r in sel], color=COLORS[s], ls=ls)
        last = [r for r in rows if r["solute"] == s and r["oil"] == "sunflower" and r["status"] == "SLE"][-1]
        b.text(last["T_C"] + 3, last["gamma"], s, color=INK2, fontsize=6.8, va="center")
    b.plot([], [], color=INK2, ls="-", label="sunflower oil")
    b.plot([], [], color=INK2, ls=OIL_STYLE["coconut"], label="coconut oil")
    b.set(yscale="log", xlabel="temperature (°C)", ylabel="γ$_s$ at saturation",
          title="(b) Activity coefficient at the solubility limit", xlim=(20, 205))
    b.legend(loc="upper right")
    fig.tight_layout(w_pad=1.5)
    save(fig, outdir, "fig4_activity_coefficients")


def fig_recipe(rows, cross, outdir):
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.7))
    for s in ("curcumin", "piperine", "vanillin"):
        for oil, ls in OIL_STYLE.items():
            sel = [r for r in rows if r["solute"] == s and r["oil"] == oil]
            a.plot([r["T_C"] for r in sel], [r["extracted_pct"] for r in sel],
                   color=COLORS[s], ls=ls)
    a.text(162, 90, "all fully dissolved", fontsize=6.5, color=INK2, ha="center")
    a.text(30, 82, "piperine", color=INK2, fontsize=6.8)
    a.text(33, 10, "curcumin", color=INK2, fontsize=6.8)
    a.text(100, 103, "vanillin (always dissolved)", color=INK2, fontsize=6.8, ha="center")
    a.plot([], [], color=INK2, ls="-", label="sunflower oil")
    a.plot([], [], color=INK2, ls=OIL_STYLE["coconut"], label="coconut oil")
    a.legend(loc="center right")
    a.set(xlabel="temperature (°C)", ylabel="share of added compound in the oil (%)",
          ylim=(0, 110), title="(a) Recipe loading: share that dissolves")
    for s in ("vanillin", "curcumin", "piperine"):
        for oil, ls in OIL_STYLE.items():
            sel = [r for r in cross if r["solute"] == s and r["oil"] == oil]
            b.plot([r["T_C"] for r in sel], [r["ratio"] for r in sel], color=COLORS[s],
                   ls=ls, marker="o", ms=3.5)
        last = [r for r in cross if r["solute"] == s and r["oil"] == "coconut"][-1]
        b.text(last["T_C"] + 1.5, last["ratio"], s, color=INK2, fontsize=6.8, va="center")
    b.axhline(1, color=INK2, lw=0.8)
    b.set(xlabel="temperature (°C)", ylabel="solubility in mixture / alone",
          title="(b) All three saturated together", xlim=(22, 90))
    fig.tight_layout(w_pad=1.5)
    save(fig, outdir, "fig5_recipe_and_mixture")


PRETTY = {"dHfus -10%": "$\\Delta H_{fus}$ −10%", "dHfus +10%": "$\\Delta H_{fus}$ +10%",
          "Tm -5 K": "$T_m$ −5 K", "Tm +5 K": "$T_m$ +5 K",
          "dCp = dHfus/Tm": "$\\Delta C_p = \\Delta H_{fus}/T_m$",
          "a_mn -10%": "all $a_{mn}$ −10%", "a_mn +10%": "all $a_{mn}$ +10%",
          "high-oleic sunflower": "high-oleic sunflower", "keto-enol groups": "keto–enol groups"}


def fig_sensitivity(sens, t_c, outdir):
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.6))
    for ax, (s, (base, cases)) in zip(axes, sens.items()):
        items = sorted(cases.items(), key=lambda kv: abs(kv[1]))
        y = np.arange(len(items))
        vals = [v for _, v in items]
        ax.barh(y, vals, color=COLORS[s], height=0.62, zorder=2)
        for yi, v in zip(y, vals):
            ax.text(v + (2 if v >= 0 else -2), yi, f"{v:+.0f}%", va="center", fontsize=6.3,
                    ha="left" if v >= 0 else "right", color=INK2)
        ax.set_yticks(y, [PRETTY.get(k, k) for k, _ in items], fontsize=6.5)
        ax.axvline(0, color=INK2, lw=0.8)
        lim = max(abs(v) for v in vals) * 1.75
        ax.set_xlim(-lim, lim)
        ax.grid(axis="y", visible=False)
        ax.set(title=f"{s.capitalize()} ({base:.0f} mg/100 g)", xlabel="change in solubility (%)")
    fig.tight_layout(w_pad=0.6)
    save(fig, outdir, "fig6_sensitivity")


def fig_convergence(brent_hist, newton_hist, outdir):
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.4))
    for lab, hist, col in brent_hist:
        a.semilogy(range(1, len(hist) + 1), np.maximum(np.abs(hist), 1e-16), color=col,
                   marker="o", ms=3.5, label=lab)
    a.set(xlabel="function evaluation", ylabel="|ln x + ln γ − ln a$_{sat}$|",
          title="(a) Brent's method, one solute")
    a.legend(loc="upper right")
    for lab, hist, ls in newton_hist:
        b.semilogy(range(len(hist)), np.maximum(hist, 1e-16), color=INK2, ls=ls, marker="o",
                   ms=3.5, label=lab)
    b.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    b.set(xlabel="Newton iteration", ylabel="max |F$_i$|",
          title="(b) Newton–Raphson, three solutes")
    b.legend(loc="upper right")
    fig.tight_layout(w_pad=1.5)
    save(fig, outdir, "fig7_convergence")


def fig_extrapolation(r5, r2, lit, ext_T, ext_S, outdir):
    """(a) data-anchored vs pure prediction with 95 % bands; (b) validation of the routes."""
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.8), gridspec_kw={"width_ratios": [1.1, 1]})
    col = COLORS["curcumin"]
    for preds, c, lab, alpha in ((r5, MUTED, "pure UNIFAC (R5)", 0.18),
                                 (r2, col, "anchored at 25 °C (R2)", 0.22)):
        T = [p.T_C for p in preds]
        a.fill_between(T, [p.ci95_mg_per_100g[0] for p in preds],
                       [p.ci95_mg_per_100g[1] for p in preds], color=c, alpha=alpha, lw=0)
        a.plot(T, [p.mg_per_100g for p in preds], color=c, label=lab)
    a.scatter([lit["anchor"][0]], [lit["anchor"][1]], s=30, color=col, edgecolors=INK, zorder=4,
              label="anchor (0.07 mg/mL)")
    a.errorbar([lit["hot"][0]], [np.mean(lit["hot"][1:])], yerr=[[np.mean(lit["hot"][1:]) - lit["hot"][1]],
               [lit["hot"][2] - np.mean(lit["hot"][1:])]], fmt="s", color=INK, ms=4.5, capsize=3,
               label="published at 100 °C (3–4 mg/mL)", zorder=4)
    a.set(yscale="log", xlabel="temperature (°C)", ylabel="curcumin in corn oil (mg/100 g)",
          title="(a) One measurement narrows the prediction")
    a.legend(loc="lower right", fontsize=6.4)
    a.set_xlim(20, 145)
    labels = list(ext_T) + list(ext_S)
    vals = [ext_T[k]["mean_abs_ln"] for k in ext_T] + [ext_S[k]["mean_abs_ln"] for k in ext_S]
    n_T, n_S = len(ext_T), len(ext_S)
    y = np.array([n_S + 1.0 + (n_T - 1 - i) for i in range(n_T)] + [float(n_S - 1 - i) for i in range(n_S)])
    colors = [INK2 if "pure" in k else COLORS["vanillin"] for k in labels]
    b.barh(y, vals, color=colors, height=0.6, zorder=2)
    for yi, v in zip(y, vals):
        b.text(v + 0.01, yi, f"×{np.exp(v):.2f}", va="center", fontsize=6.5, color=INK2)
    nice = {"anchored UNIFAC (R2)": "anchored UNIFAC (R2)", "regular solution, no groups": "regular solution, no groups",
            "constant gamma": "constant γ", "pure UNIFAC (R5)": "pure UNIFAC (R5)",
            "anchored transfer (R3)": "anchored transfer (R3)"}
    b.set_yticks(y, [nice.get(k, k) for k in labels], fontsize=6.6)
    b.text(0.005, y[0] + 0.55, f"new temperature, same solvent ({list(ext_T.values())[0]['n']} points)",
           fontsize=6.5, color=INK, va="bottom")
    b.text(0.005, y[n_T] + 0.55, f"new solvent ({list(ext_S.values())[0]['n']} points)",
           fontsize=6.5, color=INK, va="bottom")
    b.set_ylim(-0.6, y[0] + 1.2)
    b.set_xlim(0, max(vals) * 1.3)
    b.grid(axis="y", visible=False)
    b.set(xlabel="mean |ln(x$_{calc}$/x$_{exp}$)|  (label: typical factor)",
          title="(b) Error of each route, analogue data")
    fig.tight_layout(w_pad=1.2)
    save(fig, outdir, "fig8_extrapolation_uncertainty")


def fig_saturation(curve_rows, tsat_rows, example, outdir):
    """(a) How the saturation point is read off: solubility curve meets the loading line.
    (b) Saturation temperature versus amount added for every compound and oil."""
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.8))
    for oil, ls in OIL_STYLE.items():
        sel = [r for r in curve_rows if r["oil"] == oil]
        a.plot([r["T_C"] for r in sel], [r["mg"] for r in sel], color=COLORS["curcumin"], ls=ls,
               label=f"curcumin, {oil} oil")
    L = example["loading"]
    a.axhline(L, color=INK2, lw=0.9, ls=(0, (1.5, 1.5)))
    a.text(140, L * 0.55, f"loading {L:.0f} mg/100 g\n(47 mg in 15 g oil)", fontsize=6.5, color=INK2, va="top")
    for oil, ts in example["T_sat"].items():
        a.plot([ts], [L], "o", color=COLORS["curcumin"], ms=4.5, mec="white")
        off = (-30, 6) if oil == "coconut" else (5, -11)
        a.annotate(f"{ts:.0f} °C", (ts, L), xytext=off, textcoords="offset points", fontsize=6.5, color=INK2)
    a.axvspan(25, example["T_sat"]["coconut"], color=COLORS["curcumin"], alpha=0.07, lw=0)
    a.text(45, 4, "saturated:\nmore cannot mix in", fontsize=6.5, color=INK2)
    a.set_yscale("log")
    a.set(xlabel="temperature (°C)", ylabel="saturation solubility (mg per 100 g oil)",
          title="(a) Saturation point of a recipe", xlim=(25, 180))
    a.legend(loc="upper left")
    for s in ("vanillin", "curcumin", "piperine"):
        for oil, ls in OIL_STYLE.items():
            sel = [r for r in tsat_rows if r["solute"] == s and r["oil"] == oil and r["T_sat"] is not None]
            b.plot([r["loading"] for r in sel], [r["T_sat"] for r in sel], color=COLORS[s], ls=ls, marker="o", ms=2.5)
        b.plot([], [], color=COLORS[s], label=s)
    b.plot([], [], color=INK2, ls="-", label="sunflower")
    b.plot([], [], color=INK2, ls=OIL_STYLE["coconut"], label="coconut")
    b.set_xscale("log")
    b.set(xlabel="amount added (mg per 100 g oil)", ylabel="saturation temperature (°C)",
          title="(b) Temperature needed to dissolve all", ylim=(20, 205))
    b.legend(loc="upper left", ncol=2)
    fig.tight_layout(w_pad=1.5)
    save(fig, outdir, "fig9_saturation_point")
