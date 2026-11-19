# Audit: Missing data, gaps, false claims (2026-11-19)

## IIIT-Delhi data — **used**

- FlavorDB2 records for vanillin, piperine, thymol (CAS, FEMA, flavour profile, XLogP, molar mass).
- Curcumin correctly marked **absent** from FlavorDB2 (non-volatile).
- SpiceRx spice↔compound links for curcumin, piperine, thymol.
- `data/spicerx_links.csv` was missing from the repo; restored on this commit.
- Command: `python -m flavorsol database`.

## Fixed in this commit

1. Restored `python/data/spicerx_links.csv`.
2. Restored missing input JSON stubs (`02`, `05`, `08`).
3. Generated `results/solubility_single.csv` and validation/database reports.
4. Added explicit gaps documentation.

## Still open (group must supply)

| Item | Status |
|------|--------|
| Group lab measurements | Template only |
| Author names on manuscript / slides | Placeholders |
| Verified Hansen parameters from M2 | Not in code |
| Curcumin ΔHfus choice | 37.19 kJ/mol chosen; document sensitivity |
| Real headspace data | None |

## No false claims retained

- Code and docs state that FlavorDB2/SpiceRx do **not** provide Tm, ΔHfus or oil solubility.
- Two literature solubility clusters for curcumin are flagged openly.
- UNIFAC piperine amide group is noted as the least reliable fragment assignment.
