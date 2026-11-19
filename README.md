# Spice–Oil Solubility (Group 4A)

**Authors:** Rohan Jaybhave, Janil Jain, Preet Patel, Manish Kumawat, Saksham Malav  
**Topic:** Flavor Extraction – Solubility of curcumin, vanillin and piperine in cooking oils  
**Deadline:** 19 November 2026, 11:59 PM

## Submission layout (matches course guidelines)

```
├── Codes/                 ← complete code package + README
│   ├── flavorsol/         # Python package
│   ├── data/              # FlavorDB2 / SpiceRx CSVs, lab template
│   ├── inputs/            # example JSON routes
│   ├── results/           # generated tables
│   ├── demo_m3.py
│   ├── run_all.py
│   ├── requirements.txt
│   ├── README.md
│   └── GAPS_AND_AUDIT.md
├── Data/                  ← input & output data (CSV / TXT)
├── Literature/            ← place PDF copies of all references here
├── manuscript/            ← LaTeX source of the final manuscript
│   └── main.tex
├── docs/website/          ← static checker page
└── README.md              ← this file
```

Upload:
- `Codes/` → shared folder **Codes**
- `Data/` → shared folder **Data**
- `Literature/` → shared folder **Literature**
- `manuscript/` (PDF + editable source) and this README → **main project folder**

## Quick start

```bash
cd Codes
pip install -r requirements.txt
python demo_m3.py
python -m flavorsol database          # IIIT-Delhi FlavorDB2 / SpiceRx status
python -m flavorsol solubility curcumin coconut --T 25 100 150
python -m flavorsol batch --oil coconut --T 60 --oil-g 15 --curcumin 0.047 --piperine 0.05
python -m flavorsol validate
python run_all.py                     # regenerates all results/
```

## IIIT-Delhi data

FlavorDB2 and SpiceRx (Complex Systems Lab, Prof. Ganesh Bagler, **IIIT-Delhi**) supply molecule identity and flavour profiles. They do **not** supply melting points or oil solubilities. See `Codes/GAPS_AND_AUDIT.md` and `python -m flavorsol database`.

## Manuscript

Compile with any recent LaTeX distribution:

```bash
cd manuscript
pdflatex main.tex
```

Authors are already set to: Rohan Jaybhave, Janil Jain, Preet Patel, Manish Kumawat, Saksham Malav.
