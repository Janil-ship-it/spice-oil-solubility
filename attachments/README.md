# Flavour compounds in hot cooking oil — solubility from thermodynamics

CL 202 Chemical Engineering Thermodynamics, Project 4 (Group 4A), IIT Gandhinagar.

This Python tool predicts how much **curcumin** (turmeric), **vanillin** (vanilla) and
**piperine** (black pepper) can dissolve in **sunflower** and **coconut oil** between
25 and 200 °C, alone or together, and how much of a recipe's spice compound actually ends
up in the oil.

## Quick start

```bash
pip install -r requirements.txt        # numpy, scipy, matplotlib (pytest optional)
python demo_m3.py                      # step-by-step demo for the meeting (~10 s)
python run_all.py                      # every number, table and figure (~2 min)
python -m pytest -q                    # 30 tests
python -m flavorsol database           # FlavorDB2/SpiceRx records loaded
python -m flavorsol predict inputs/03_anchor_same_oil.json
python -m flavorsol solubility curcumin coconut --T 25 100 150
python -m flavorsol mixture --oil sunflower --T 60 --loadings 100 314 333
python -m flavorsol batch --oil coconut --T 60 --oil-g 15 --curcumin 0.047 --piperine 0.05
```

Google Colab: upload the folder, then `!pip install -r requirements.txt` and `!python demo_m3.py`.

## The physics in four lines

1. Solid–liquid equilibrium: `ln(x·γ) = −(ΔH_fus/R)(1/T − 1/T_m) + (ΔC_p/R)[T_m/T − 1 − ln(T_m/T)]`
2. γ from **UNIFAC** (original, Hansen 1991 parameters from DDBST). Each oil is one average
   triglyceride built from its fatty-acid composition.
3. γ depends on x, so the equation is solved by **Brent's method** (one compound) or
   **damped Newton–Raphson** (several compounds in one oil, with an active-set loop that decides
   which compounds saturate and which dissolve completely).
4. Results are converted to **mg per 100 g oil** and compared with sensory thresholds.

## Inputs: mandatory vs optional (from the meeting 2 feedback)

| Mandatory | Optional — supplying it changes the route |
|---|---|
| `solute` (library name, or `M` + `Tm` + `groups` for a new compound) | `dHfus`, `dCp` — else Walden's rule (ΔS_fus = 56.5 J/mol K, ±25 %) and ΔC_p = 0 |
| `oil` (name, or `oil_composition` in % fatty acids) | `gamma` (+ `gamma_T_C`) → **R1**: UNIFAC is skipped |
| `T_C` (one value or a list) | `measurements` in the same oil → **R2** (anchored UNIFAC) |
| | `measurements` in another solvent → **R3** (transferred with UNIFAC) |
| | headspace `measurements` (partial pressure) → **R4** (Henry's law) |
| | `uncertainties`, `oil_grade`, `loading_mg_per_100g`, `monte_carlo`, `T_extrapolation` |

Nothing optional given → **R5**, pure UNIFAC prediction. Examples of every route are in `inputs/`
(including deliberately bad inputs in `07_bad_inputs.json`).

Every run also:
- **checks inputs** (units, ranges, Walden's rule, mole fractions, composition sums, data above the
  ideal-solution limit) and stops on errors;
- **flags breakdown instead of dissolution**: curcumin degradation in hot oil, vanillin above its
  flash point, unverified piperine stability, oil smoke points, solute melting, extrapolation outside
  UNIFAC's fitted range or far from a measurement;
- **propagates uncertainties** (temperature, T_m, ΔH_fus, ΔC_p, measurement and model error) with
  numerical derivatives (first-order), optionally checked by Monte Carlo, and reports which input
  dominates;
- reports **Gibbs–Helmholtz derivatives**: partial molar excess enthalpy and entropy, d ln x/dT
  (implicit-function theorem) and the apparent enthalpy of dissolution.

## Folder map and who owns what

| Path | What it does | Owner |
|---|---|---|
| `flavorsol/data/compounds.py` | T_m, ΔH_fus, M, UNIFAC groups, sources | Member 1 – data |
| `flavorsol/data/oils.py` | fatty-acid compositions → triglyceride pseudo-component | Member 1 |
| `flavorsol/data/unifac_params.py` | UNIFAC R, Q and a_mn (DDBST) | Member 1 |
| `flavorsol/data/spices.py` | spice → compound map; reads the FlavorDB2 records | Member 1 |
| `data/flavordb_*.csv`, `data/spicerx_links.csv`, `tools/fetch_flavordb.py` | FlavorDB2 and SpiceRx records (IIIT-Delhi) and the refresh script | Member 1 |
| `flavorsol/unifac.py` | UNIFAC activity coefficients | Member 2 – single-solute solver |
| `flavorsol/sle.py` | ideal solubility, Brent solver, liquid–liquid split above/near T_m | Member 2 |
| `flavorsol/derivatives.py` | Gibbs–Helmholtz, d ln x/dT, excess enthalpy/entropy | Member 2 |
| `flavorsol/mixture.py` | multi-solute Newton–Raphson + mass-limited active set | Member 3 – mixture solver |
| `flavorsol/engine.py` | mandatory/optional inputs, routes R1–R5, fallbacks, uncertainty | Members 3 & 4 |
| `flavorsol/operating.py` | grams of oil and compound → grams dissolved / left solid; spice extraction upper bound | Member 3 |
| `flavorsol/hsp.py` | Hansen distance, RED, alternative HSP activity coefficient — screening only, never combined with UNIFAC | Member 1 |
| `flavorsol/checks.py` | boundary and degradation checks | Member 4 – validation |
| `flavorsol/validation.py` | benchmark, literature and lab comparisons, route validation | Member 4 |
| `data/literature_solubility.csv`, `data/lab_measurements_TEMPLATE.csv` | validation data | Member 4 |
| `flavorsol/analysis.py`, `flavorsol/plotting.py`, `data/sensory.py` | sweeps, recipe, thresholds, sensitivity, figures | Member 5 – analysis |
| `run_all.py`, `demo_m3.py`, `tests/` | reproduce everything / meeting demo / tests | all |

## FlavorDB and SpiceRx (IIIT-Delhi), as advised in meeting 1

FlavorDB2 and SpiceRx (Complex Systems Lab, G. Bagler) give each compound's identity, flavour
profile and the spices it occurs in. They hold no melting or solubility data, so they are the front
end of the tool. All 935 FlavorDB2 entries were scanned and SpiceRx's chemical index searched on
2026-11-09; the results are in three data files that the code reads:

| File | Contents |
|---|---|
| `data/flavordb_records.csv` | FlavorDB2 molecule records: CAS, FEMA number and profile, FooDB profile, odour, taste, functional groups, molar mass, XLogP |
| `data/flavordb_spice_entities.csv` | FlavorDB2 entries for our spices: Turmeric 341, Pepper 339, Vanilla 245, Carom Seed (ajwain) 381, White Pepper 399, Fenugreek 258, Dill 256 |
| `data/spicerx_links.csv` | SpiceRx MeSH IDs and the spices each compound occurs in |

| Compound | FlavorDB2 | SpiceRx |
|---|---|---|
| vanillin | record from Vanilla (245): FEMA 3107, vanilla/chocolate/sweet/creamy, XLogP 1.2 | not listed |
| piperine | record from Pepper (339): FEMA 2909, "pepper", XLogP 3.5; not in White Pepper | Black Pepper, Long Pepper |
| thymol | record from Carom Seed (381): FEMA 3066, herbal/thyme/phenolic, XLogP 3.3 | Ajowan + 12 other spices |
| curcumin | **absent from all 935 entries** (Turmeric lists curcumene, curcumenol) | Turmeric, Zedoary, Curry Leaf |

How the code uses them:
- every prediction prints the PubChem CID and flavour profile with its source, and warns if the
  molar mass disagrees with FlavorDB2 (all three agree to 0.01 g/mol);
- `validation.database_check()` compares FlavorDB2's XLogP with our UNIFAC log Kow
  (`results/validation_flavordb.csv`): vanillin 1.35 vs 1.2, thymol 3.75 vs 3.3, but piperine
  5.79 vs 3.5, which flags piperine's amide group assignment as the weakest;
- `python -m flavorsol database` shows everything loaded.

Two cautions: FlavorDB2 covers flavour volatiles, so non-volatile colour/pungency compounds
(curcumin) must come from SpiceRx; and SpiceRx's formula and molar mass for piperine and thymol
contradict its own InChI, so only its identity fields and spice links are used.
`tools/fetch_flavordb.py` repeats the FlavorDB2 scan from any computer with normal internet access.

## Three quantities, kept separate

1. **Saturation solubility** x_sat(T, oil) — needs only temperature and oil (`solubility`, `predict`).
2. **Dissolved amount at a loading** — needs grams of oil and of each compound; dissolved = min(added, capacity),
   solved jointly for mixtures (`batch`).
3. **Extraction yield from a spice** — needs spice mass, compound content and a release factor η (0–1) that only
   experiments can give; η = 1 is the thermodynamic upper bound (`operating.batch(..., spice=...)`).

## Adding your lab results

Copy `data/lab_measurements_TEMPLATE.csv` to `data/lab_measurements.csv`, one row per measurement
(`solute, oil, T_C, solubility_mg_per_100g, uncertainty_mg_per_100g, ...`) and rerun `run_all.py`.
It writes `results/validation_lab.csv` and `results/tables/lab.tex`; switch on `\labdatatrue` in
`manuscript/main.tex` to print the comparison. A single lab point can also be used as an R2 anchor
(see `inputs/03_anchor_same_oil.json`).

## Adding a compound or an oil

- **Compound**: find it in SpiceRx/FlavorDB (identity, PubChem CID, flavour), then add `M`, `Tm`,
  UNIFAC `groups` (and `dHfus` if known) to `data/compounds.py` or pass them in a JSON input.
- **Oil**: add its fatty-acid composition to `OIL_COMPOSITIONS` in `data/oils.py`, or pass
  `oil_composition` in the JSON.

## Main limitations

UNIFAC is a group-contribution estimate (typical error a factor 1.5–2.5 here); one triglyceride
stands in for the whole oil; the ACOH–DMF interaction parameter is not published (set to 0);
curcumin's tautomerism and polymorphism are not resolved; heating also degrades curcumin and
evaporates vanillin, which equilibrium alone cannot capture.
