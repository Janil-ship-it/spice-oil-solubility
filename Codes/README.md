# Code package: Spice compounds in hot oil (Group 4A)

CL 202 Chemical Engineering Thermodynamics, IIT Gandhinagar, Project 4 (Flavour extraction).
Authors: Rohan Jaybhave, Janil Jain, Preet Patel, Saksham Malav, Manish Kumawat.

The code predicts how much **curcumin** (turmeric), **vanillin** (vanilla) and **piperine**
(black pepper) dissolve in **sunflower** and **coconut oil** from 25 to 200 °C, alone or
together, and finds the **saturation point**: the temperature and amount at which the oil
can take no more of a compound ("at this temperature saturation is achieved and more mixing
cannot be done").

## Folder structure

```
Codes/
├── README.md                 this file
├── python/                   the thermodynamic model (Python package "flavorsol")
│   ├── flavorsol/
│   │   ├── data/             compound data (Tm, dHfus, UNIFAC groups), oils, UNIFAC parameters, spices
│   │   ├── unifac.py         UNIFAC activity coefficients
│   │   ├── sle.py            solid-liquid equilibrium, Brent solver, liquid-liquid split above Tm
│   │   ├── mixture.py        several compounds in one oil: Newton-Raphson + active set
│   │   ├── saturation.py     SATURATION POINT: capacity, saturated or not, saturation temperature
│   │   ├── operating.py      grams of oil and compound -> grams dissolved / undissolved
│   │   ├── engine.py         mandatory vs optional inputs, routes R1-R5, uncertainty
│   │   ├── checks.py         degradation, flash point, smoke point, range checks
│   │   ├── derivatives.py    Gibbs-Helmholtz, d ln x/dT, excess enthalpy
│   │   ├── validation.py     comparison with literature data
│   │   ├── analysis.py       sweeps, recipes, sensitivity
│   │   ├── plotting.py       figures of the manuscript
│   │   └── __main__.py       command-line interface
│   ├── data/                 FlavorDB2 / SpiceRx records, literature solubilities, lab template
│   ├── inputs/               example JSON inputs for routes R1-R5 (one deliberately wrong)
│   ├── tests/                37 automated tests (pytest)
│   ├── tools/                FlavorDB refresh, export of model data for the website
│   ├── run_all.py            reproduces every table, figure and the Data folder
│   ├── export_data.py        writes ../../Data (inputs and outputs of all scenarios)
│   ├── demo_m3.py            short step-by-step demo
│   ├── api.py                optional HTTP API (FastAPI) used by the website's "Full report" tab
│   └── requirements.txt
├── website/                  the web application (React + TypeScript, Vite)
│   ├── dist/index.html       READY-TO-USE WEBSITE: one file, open it in any browser
│   ├── src/App.tsx           page layout, Single compound tab
│   ├── src/components/       Mixtures, Methodology and Full report tabs
│   ├── src/lib/engine.ts     browser copy of UNIFAC + mixture solver + saturation rules
│   ├── src/data/model.json   model data exported from Python (tools/export_web_data.py)
│   └── tests/                check that the browser solver equals the Python model
└── notebook/Group4A_solubility_model.ipynb   Colab notebook walkthrough
```

## How the saturation point is found

1. **Saturation solubility S(T).** Solve the SLE equation
   `ln(x·γ(x,T)) = -(ΔHfus/R)(1/T - 1/Tm)` for x (Brent's method; γ from UNIFAC) and convert
   to mg per 100 g oil.
2. **Saturation point of the pan:** capacity = S(T) × grams of oil / 100.
3. **Compare with the amount added.** If added ≥ capacity, saturation is achieved at this
   temperature and more mixing cannot be done; the extra stays undissolved (solid, or a second
   liquid layer above the melting point). Otherwise everything dissolves and the remaining
   capacity is reported.
4. **Saturation temperature:** S(T) rises monotonically with T (dissolution is endothermic), so
   the root of S(T_sat) = loading is unique; it is found on a 1 °C grid and refined with Brent's
   method. Below T_sat the oil is saturated, above it everything dissolves.
5. **Mixtures:** the coupled equations for all compounds are solved together (Newton-Raphson +
   active set). A compound is saturated when its equation is active and less than the added
   amount dissolves. The temperature at which the whole recipe dissolves is found by bisection.

## Running the Python model

Python 3.10 or newer.

```bash
cd Codes/python
pip install -r requirements.txt
```

| Command | What it does |
|---|---|
| `python -m flavorsol saturation curcumin coconut --T 60 --oil-g 15 --added-g 0.047` | saturation statement, one compound |
| `python -m flavorsol mixsat --oil coconut --T 60 --oil-g 15 --curcumin 0.047 --piperine 0.05 --vanillin 0.015` | saturation statements, mixture |
| `python -m flavorsol solubility curcumin sunflower --T 25 60 100 140` | solubility table |
| `python -m flavorsol mixture --oil sunflower --T 60 --loadings 100 314 333` | three compounds together (mg per 100 g oil) |
| `python -m flavorsol batch --oil coconut --T 60 --oil-g 15 --curcumin 0.047 --piperine 0.05` | grams dissolved / left |
| `python -m flavorsol predict inputs/03_anchor_same_oil.json` | full workflow with routes R1-R5 and uncertainty |
| `python -m flavorsol validate` | comparison with literature |
| `python -m flavorsol database` | FlavorDB2 / SpiceRx records loaded |
| `python demo_m3.py` | step-by-step demo |
| `python -m pytest -q` | 37 tests |
| `python run_all.py` | everything: results/, figures/ and the ../../Data folder (about 10 min) |

Example output:

```
$ python -m flavorsol saturation curcumin coconut --T 60 --oil-g 15 --added-g 0.047
At 60.0 °C saturation is achieved: 15 g of coconut oil can hold only 6.14 mg of curcumin
(40.9 mg per 100 g oil). More mixing cannot be done; the extra 40.9 mg stays as undissolved
solid. To dissolve all 47.0 mg, heat to at least 115.2 °C.
```

## Website

**Open `Codes/website/dist/index.html` in any browser** (double-click; no installation, works
offline). Online copy (after this branch is merged into main):
https://htmlpreview.github.io/?https://github.com/Janil-ship-it/spice-oil-solubility/blob/main/Codes/website/dist/index.html

Tabs:

| Tab | Content |
|---|---|
| Single compound | compound, oil (sunflower or coconut), temperature, grams of oil and compound; saturation statement, capacity, amount that cannot mix in, saturation temperature marked on the solubility curve |
| Mixtures | recipe in grams with all three compounds: saturation statement for each compound, temperature at which the whole recipe dissolves, share dissolved vs temperature; and an "all in excess" mode showing co-solvency |
| Methodology | data used, equations, how the saturation point is found, thermodynamic consistency checks, validation and limitations |
| Full report (Python) | calls the Python engine (routes R1-R5, uncertainty); needs `uvicorn api:app --port 8000` running in `Codes/python` |

The Single compound and Mixtures tabs run entirely in the browser: single-compound values are
the Python solubility curves (1 °C steps) and the mixture equations are solved by a line-by-line
TypeScript port of the Python UNIFAC and Newton-Raphson code (`src/lib/engine.ts`). `npm test`
checks it against the Python model (48 mixture cases agree to 1e-13, saturation temperatures to
0.001 °C).

Rebuilding the website (Node.js 20+):

```bash
cd Codes/python && python tools/export_web_data.py && python tools/make_web_reference.py
cd ../website
npm install
npm test            # browser solver = Python model
npm run dev         # development server on http://localhost:8080
npm run build       # writes dist/index.html (single file)
```

## Data

All inputs and outputs (CSV/TXT) for every scenario are in the top-level `Data/` folder; see
`Data/README.txt`. They are regenerated by `python run_all.py`.

## Main limitations

UNIFAC is a group-contribution estimate (typical error a factor 1.5 to 2.5 here); one average
triglyceride stands in for each oil; the ACOH-DMF interaction parameter is unpublished (set to
0); frying temperatures are above the range where most UNIFAC parameters were fitted; heating
also degrades curcumin and evaporates vanillin, which equilibrium cannot capture (flagged).
The saturation point is therefore a thermodynamic upper bound on what the oil takes up.

Licence: MIT (see LICENSE in the repository root).
