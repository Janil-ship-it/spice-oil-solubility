# Spice Compounds in Hot Oil: Solubility and Saturation Point (Group 4A)

CL 202 Chemical Engineering Thermodynamics, IIT Gandhinagar, Project 4 (Flavour extraction).

**Authors:** Rohan Jaybhave, Janil Jain, Preet Patel, Saksham Malav, Manish Kumawat

We predict how much curcumin (turmeric), vanillin (vanilla) and piperine (black pepper) can
dissolve in sunflower and coconut oil from 25 to 200 °C, alone and in mixtures, using
solid-liquid equilibrium with UNIFAC activity coefficients. For every case the code and the
website state whether **saturation is achieved at the chosen temperature, so that more mixing
cannot be done**, how much stays undissolved, and the temperature above which everything dissolves.

## Submission folder layout

```
.
├── Group4A_Manuscript.pdf     final manuscript (JCED / ACS format)
├── Group4A_Manuscript.tex     editable LaTeX source (pdflatex, achemso class)
├── Group4A_Manuscript.docx    editable Word version
├── manuscript_figures/        figures used by the manuscript (PDF for LaTeX, PNG for Word)
├── Presentations/             project slides
├── Codes/                     complete code package + README (Python model, website, notebook)
├── Data/                      input and output data of every scenario (CSV / TXT) + README
├── Literature/                reference PDFs + README (list with DOIs)
├── LICENSE
└── README.md                  this file
```

Upload `Codes/`, `Data/` and `Literature/` to the shared subfolders of the same names; the
manuscript files (PDF, TEX, DOCX, manuscript_figures) and this README go to the main folder.

## Website

Open **`Codes/website/dist/index.html`** in any browser (single file, works offline).
Tabs: Single compound, Mixtures, Methodology, Full report.
Online copy (after merging into main):
https://htmlpreview.github.io/?https://github.com/Janil-ship-it/spice-oil-solubility/blob/main/Codes/website/dist/index.html

## Quick start (Python)

```bash
cd Codes/python
pip install -r requirements.txt
python -m flavorsol saturation curcumin coconut --T 60 --oil-g 15 --added-g 0.047
python -m flavorsol mixsat --oil coconut --T 60 --curcumin 0.047 --piperine 0.05 --vanillin 0.015
python -m pytest -q          # 37 tests
python run_all.py            # all results, figures and the Data folder
```

## Key results

| | sunflower oil | coconut oil |
|---|---|---|
| Curcumin solubility at 25 / 100 / 180 °C (mg per 100 g) | 3.3 / 97 / 1430 | 8.8 / 186 / 2450 |
| Tadka (47 mg curcumin in 15 g oil): saturated below | 132.8 °C | 115.2 °C |
| Tadka (50 mg piperine in 15 g oil): saturated below | 30.7 °C | 32.8 °C |
| Tadka vanillin (15 mg in 15 g oil) | never saturated | never saturated |
| Whole tadka recipe (V + C + P together) dissolves above | 131.6 °C | 114.1 °C |

Validation: the UNIFAC code reproduces the DDBST reference to 6e-11; published oil
solubilities are matched within a factor of 2.2 on average (ideal solution: factor of about 20).
Every prediction satisfies γ ≥ 1, x ≤ x_ideal and a solubility that rises with temperature.

See `Codes/README.md` for the code, `Data/README.txt` for the data files and the Methodology
tab of the website for the method.
