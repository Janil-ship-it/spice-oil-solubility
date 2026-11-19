# Spice–Oil Solubility Model (Group 4A)

Predicts the solubility of curcumin, vanillin and piperine in cooking oils
(sunflower, coconut, …) at cooking temperatures using solid–liquid equilibrium
+ UNIFAC activity coefficients.

## IIIT-Delhi databases (FlavorDB2 / SpiceRx)

These databases (Complex Systems Lab, Prof. Ganesh Bagler, IIIT-Delhi) supply
**identity and flavour profiles**, not melting points or oil solubilities.

| Compound  | FlavorDB2                         | SpiceRx                          |
|-----------|-----------------------------------|----------------------------------|
| vanillin  | Full record (FEMA 3107, XLogP 1.2)| Not listed                       |
| piperine  | Full record (FEMA 2909, XLogP 3.5)| Black Pepper, Long Pepper        |
| thymol    | Full record (FEMA 3066, XLogP 3.3)| Ajowan + 12 other spices         |
| curcumin  | **Absent** (not a volatile)       | Turmeric, Zedoary, Curry Leaf    |

Thermodynamic inputs (Tm, ΔHfus, UNIFAC groups) come from NIST / DSC literature.
See `python -m flavorsol database`.

## Quick start

```bash
cd python
pip install -r requirements.txt
python demo_m3.py                          # live M3 walkthrough
python -m flavorsol solubility curcumin coconut --T 25 100 150
python -m flavorsol mixture --oil sunflower --T 60
python -m flavorsol batch --oil coconut --T 60 --oil-g 15 --curcumin 0.047 --piperine 0.05
python -m flavorsol database               # IIIT-Delhi records
python -m flavorsol validate
python run_all.py                          # full study → results/
```

## Gaps / honesty notes

1. **Lab data**: `data/lab_measurements.csv` is a template — replace with your measurements.
2. **Curcumin ΔHfus**: literature range 37–49 kJ/mol; we use 37.19 kJ/mol (Ranjbar 2024). Sensitive.
3. **Two clusters of published curcumin–oil solubilities**: model matches the low cluster (~0.1 mg/mL at 25 °C). High values (~5–12 mg/mL) may reflect supersaturation or degraded material.
4. **Headspace route**: only synthetic examples; no real headspace data yet.
5. **Hansen parameters**: screening only; verified values not available for all compounds.

## Structure

```
python/
  flavorsol/          # package
    data/             # compounds, oils, spices, UNIFAC params
    engine.py         # decision pathways (mandatory vs optional inputs)
    sle.py            # solid–liquid equilibrium + Brent solver
    unifac.py         # UNIFAC activity coefficients
    mixture.py        # Newton–Raphson multi-solute
    checks.py         # thermal degradation / smoke-point flags
    ...
  data/               # FlavorDB / SpiceRx CSVs, literature solubility
  inputs/             # example JSON requests (routes R1–R5)
  results/            # generated tables
  demo_m3.py
  run_all.py
```

## Licence

MIT (see LICENSE in repository root).
