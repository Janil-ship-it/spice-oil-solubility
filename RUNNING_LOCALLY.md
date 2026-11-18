# Running the website locally (Windows)

The site has two parts that run side by side:

| Part | Folder | Port | What it is |
|---|---|---|---|
| Python model (API) | `python/` | 8000 | The project's `flavorsol` package (UNIFAC, SLE, mixture solver, engine) behind a small FastAPI server (`python/api.py`) |
| Website | this folder | 8080 | React / TanStack Start app. Its dev server forwards every `/api/...` request to port 8000 |

If the Python model is not running, the **Single compound** tab still works from precomputed values
(`src/data/curves.json`) and says "Offline". **Mixture & recipe** and **Full report** need the model.

## Prerequisites

- **Node.js 20 or newer**: check with `node -v`
- **Python 3.10 or newer**: check with `python --version`

## 1. One-time setup

Open PowerShell:

```powershell
cd C:\Users\Admin\Downloads\spice-oil-solubility-website

# Python model
cd python
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
cd ..

# Website
npm install
```

If `Activate.ps1` is blocked ("running scripts is disabled"), run this once and try again:
`Set-ExecutionPolicy -Scope CurrentUser RemoteSigned`

## 2. Start it (two terminals)

**Terminal 1: Python model**

```powershell
cd C:\Users\Admin\Downloads\spice-oil-solubility-website\python
.\.venv\Scripts\Activate.ps1
uvicorn api:app --host 127.0.0.1 --port 8000 --reload
```

Wait for `Uvicorn running on http://127.0.0.1:8000`. To check it, open
http://127.0.0.1:8000/api/health (should show `{"ok":true}`). The interactive API docs are at
http://127.0.0.1:8000/docs.

> Shortcut: from the website folder with the venv active, `npm run api` does the same.

**Terminal 2: website**

```powershell
cd C:\Users\Admin\Downloads\spice-oil-solubility-website
npm run dev
```

Open **http://localhost:8080**. Stop either one with `Ctrl+C`.

## What's on the page

| Tab | Python function | Notebook section |
|---|---|---|
| Single compound | `sle.solubility`, `checks.check_temperature` | §3 `python -m flavorsol solubility curcumin coconut --T 25 100 150` |
| Mixture & recipe: Per 100 g oil | `mixture.solve_mixture` | §4 `python -m flavorsol mixture --oil sunflower --T 60 --loadings 100 314 333` |
| Mixture & recipe: Recipe in grams | `operating.batch` | §5 `python -m flavorsol batch --oil coconut --T 60 --oil-g 15 --curcumin 0.047 --piperine 0.05` |
| Full report | `engine.predict` + `format_report` | §6 `python -m flavorsol predict my_case.json` |

## The original command-line tools still work

With the venv active, from `python/`:

```powershell
python -m pytest -q                      # the 30 tests
python demo_m3.py                        # step-by-step demo
python -m flavorsol mixture --oil sunflower --T 60 --loadings 100 314 333
python -m flavorsol predict inputs\03_anchor_same_oil.json
python run_all.py                        # every table and figure (~2 min, writes results\ and figures\)
```

## Troubleshooting

- **"Port 8080 is already in use"**: another dev server is running. Close it, or find it with
  `netstat -ano | findstr :8080` and stop it with `taskkill /PID <pid> /F`.
- **Page says "Offline" or "Python model is not reachable"**: Terminal 1 isn't running, or it is
  on a different port. The website expects `127.0.0.1:8000`.
- **`npm run dev` fails with an error about a native binding (`rolldown`, `lightningcss`,
  `@tailwindcss/oxide`)**: `node_modules` came from another OS. Delete it and reinstall:
  `Remove-Item -Recurse -Force node_modules; npm install`.
- **`uvicorn` is not recognised**: the venv isn't active. Run `.\.venv\Scripts\Activate.ps1`
  first, or use `.\.venv\Scripts\python.exe -m uvicorn api:app --port 8000`.
- **`pip install` fails building numpy/scipy**: your Python is too new for the available wheels.
  Install Python 3.12 from python.org and create the venv with `py -3.12 -m venv .venv`.
