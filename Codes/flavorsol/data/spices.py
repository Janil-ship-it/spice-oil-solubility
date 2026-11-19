"""
Spice -> flavour compound map, linked to FlavorDB and SpiceRx (Member 1 - Data).

FlavorDB (25,595 flavour molecules; Garg et al., Nucleic Acids Res. 2018) and
SpiceRx (188 spices/herbs and 866 phytochemicals; Rakhi et al., PLoS ONE 2018)
are maintained by Ganesh Bagler's Complex Systems Lab at IIIT-Delhi
(https://cosylab.iiitd.edu.in). They give each molecule's identity (PubChem
CID, SMILES), its flavour profile and which spices contain it, but not
oil-solubility or melting data. We therefore use them as the FRONT END:

    spice (SpiceRx)  ->  compound + PubChem CID + flavour profile (FlavorDB)
                     ->  Tm, dHfus, UNIFAC groups (this project, compounds.py)
                     ->  solubility in oil (engine.py)

Anything a database record lacks is reported by `missing_for_prediction`,
and the engine fills what it safely can (e.g. dHfus by Walden's rule).

What we took from the two databases (scanned 2026-11-09, all 935 FlavorDB2
entities and the SpiceRx chemical index):

  data/flavordb_records.csv       FlavorDB2 molecule records: CAS, FEMA number and
                                  flavour profile, FooDB profile, odour, taste,
                                  functional groups, molar mass, XLogP
                                  (vanillin, thymol, piperine). Curcumin is in
                                  none of the 935 entities.
  data/flavordb_spice_entities.csv FlavorDB2 entities for our spices
                                  (Turmeric 341, Pepper 339, Vanilla 245, Carom
                                  Seed 381, White Pepper 399, ...)
  data/spicerx_links.csv          SpiceRx: which spices contain each compound,
                                  MeSH IDs (curcumin, piperine, thymol; vanillin
                                  is not in SpiceRx)

Findings that matter for the model:
  * FlavorDB2 catalogues flavour volatiles. Its Turmeric entry lists curcumene
    and curcumenol but not curcumin; White Pepper lacks piperine (Pepper has it).
    Curcumin's spice link therefore comes from SpiceRx.
  * SpiceRx's formula and molar mass for piperine (C20H16O6, 352.3) and thymol
    (C6H14S2) contradict its own InChI; only its identity and spice links are
    used. FlavorDB2's molar masses agree with our library and are cross-checked
    on every prediction.
  * FlavorDB2 XLogP (vanillin 1.2, thymol 3.3, piperine 3.5) is a second
    independent check of the UNIFAC-based log Kow.

tools/fetch_flavordb.py refreshes the FlavorDB2 file. FALLBACK_FLAVOR holds
general descriptors used only where FlavorDB2 has no record (curcumin);
`flavor_profile()` always names the source.
"""
import csv
import os

from .compounds import ALL_SOLUTES, Compound

_DATA = os.path.join(os.path.dirname(__file__), "..", "..", "data")
RECORDS_CSV = os.path.join(_DATA, "flavordb_records.csv")
SPICERX_CSV = os.path.join(_DATA, "spicerx_links.csv")
ENTITIES_CSV = os.path.join(_DATA, "flavordb_spice_entities.csv")

PUBCHEM_CID = {"vanillin": 1183, "curcumin": 969516, "piperine": 638024, "thymol": 6989}

# General descriptors, used only when no FlavorDB record has been retrieved.
FALLBACK_FLAVOR = {
    "vanillin": ["vanilla", "sweet", "creamy"],
    "curcumin": ["earthy", "slightly bitter", "colour (yellow)"],
    "piperine": ["pungent", "peppery", "bitter"],
    "thymol":   ["herbal", "thyme-like", "phenolic"],
}


def load_records(path=RECORDS_CSV):
    """Read data/flavordb_records.csv -> {compound: row dict}. Empty dict if absent."""
    if not os.path.exists(path):
        return {}
    with open(path, newline="", encoding="utf-8") as f:
        return {r["compound"].strip().lower(): r for r in csv.DictReader(f)}


def _read(path, key):
    if not os.path.exists(path):
        return {}
    with open(path, newline="", encoding="utf-8") as f:
        return {r[key].strip().lower(): r for r in csv.DictReader(f)}


FLAVORDB = load_records()
SPICERX = _read(SPICERX_CSV, "compound")
FLAVORDB_ENTITIES = _read(ENTITIES_CSV, "name")


def spices_containing(name):
    """Spices that contain a compound, from SpiceRx and FlavorDB2 (each source listed)."""
    key = name.lower()
    out = {"SpiceRx": [], "FlavorDB2": []}
    if SPICERX.get(key, {}).get("spices_in_spicerx"):
        out["SpiceRx"] = [x.strip() for x in SPICERX[key]["spices_in_spicerx"].split(";")]
    out["FlavorDB2"] = [e["name"] for e in FLAVORDB_ENTITIES.values()
                        if key in [t.strip() for t in e["target_compounds"].split(";")]]
    return out


def flavor_profile(name):
    """Flavour descriptors for a compound and where they came from.

    Returns (descriptors, source). FooDB/FEMA profiles from a retrieved FlavorDB
    record take priority; otherwise the general fallback list is returned and
    the source says so, so a report never passes off our words as FlavorDB's.
    """
    key = name.lower()
    r = FLAVORDB.get(key, {})
    if r.get("status") == "retrieved":
        words = []
        for field in ("fooddb_flavor_profile", "fema_flavor_profile"):
            for w in (r.get(field) or "").replace("@", ",").split(","):
                w = w.strip().lower()
                if w and w not in words:
                    words.append(w)
        if words:
            return words, r.get("flavordb_source") or "FlavorDB2"
    why = ("not in FlavorDB2" if r.get("status") == "absent_from_flavordb"
           else "no FlavorDB2 record retrieved")
    return FALLBACK_FLAVOR.get(key, []), f"general descriptor ({why})"


COMPOUND_IDS = {
    k: {"pubchem_cid": cid,
        "flavor": flavor_profile(k)[0],
        "flavor_source": flavor_profile(k)[1],
        "flavordb_status": FLAVORDB.get(k, {}).get("status", "not_retrieved"),
        "spicerx_mesh": SPICERX.get(k, {}).get("spicerx_mesh_id", "")}
    for k, cid in PUBCHEM_CID.items()
}

SPICES = {
    "turmeric": {"compounds": {"curcumin": 0.0314},
                 "content_source": "Tayyem et al., Nutr. Cancer 2006, 55, 126 (mean of turmeric powders)"},
    "black pepper": {"compounds": {"piperine": 0.05},
                     "content_source": "Gorgani et al., Compr. Rev. Food Sci. Food Saf. 2017, 16, 124 (about 2-7 %)"},
    "vanilla": {"compounds": {"vanillin": None},
                "content_source": "content varies with curing; enter your own value"},
    "ajwain": {"compounds": {"thymol": None},
               "content_source": "thymol is ajwain's main volatile; enter your own content"},
}


def database_status():
    """One line per compound: what the IIIT-Delhi databases gave us so far."""
    lines = []
    for k, info in COMPOUND_IDS.items():
        r = FLAVORDB.get(k, {})
        lines.append(f"{k:<9} CID {info['pubchem_cid']:<7} {info['flavordb_status']:<13} "
                     f"flavour: {', '.join(info['flavor'])}  [{info['flavor_source']}]")
        if r.get("fema_number"):
            lines.append(f"{'':<9} FEMA {r['fema_number']}; CAS {r.get('cas', '')}; M {r.get('molecular_weight', '')}; "
                         f"XLogP {r.get('xlogp', '')}; odour: {r.get('odor') or '-'}")
        where = spices_containing(k)
        lines.append(f"{'':<9} in spices - SpiceRx: {', '.join(where['SpiceRx']) or 'not listed'}; "
                     f"FlavorDB2: {', '.join(where['FlavorDB2']) or 'not listed'}")
    return "\n".join(lines)

# Example of a compound added from FlavorDB/SpiceRx with only partial data:
# Tm and groups known, dHfus NOT known -> engine applies Walden's rule.
# 2-isopropyl-5-methylphenol; Tm ~ 49.6 C (PubChem; check before use).
THYMOL = Compound(
    name="thymol", formula="C10H14O", M=150.22, Tm=322.75, dHfus=None,
    groups={"ACH": 3, "ACCH3": 1, "ACCH": 1, "ACOH": 1, "CH3": 2},
    source="Tm from PubChem CID 6989 (verify); dHfus deliberately left empty to demonstrate the fallback.",
)

EXTRA_COMPOUNDS = {"thymol": THYMOL}

REQUIRED_FOR_UNIFAC = ("M", "Tm", "groups")


def compounds_in_spice(spice: str):
    return SPICES[spice.lower()]["compounds"]


def missing_for_prediction(c: Compound):
    """List which thermodynamic inputs a compound record is missing."""
    missing = [k for k in REQUIRED_FOR_UNIFAC if not getattr(c, k, None)]
    if c.dHfus is None:
        missing.append("dHfus (optional: Walden's rule fallback)")
    return missing


def load_flavordb_csv(path):
    """Read a FlavorDB/SpiceRx export saved as CSV.

    Expected (case-insensitive) columns, any extra columns are ignored:
        common name | pubchem id | molecular weight | flavor profile
    Returns a list of dicts with the fields we use plus which thermodynamic
    inputs are still needed before the engine can predict solubility.
    """
    out = []
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            r = {k.strip().lower(): (v or "").strip() for k, v in row.items()}
            name = r.get("common name") or r.get("name") or ""
            known = ALL_SOLUTES.get(name.lower()) or EXTRA_COMPOUNDS.get(name.lower())
            out.append({
                "name": name, "pubchem_cid": r.get("pubchem id"),
                "M": float(r["molecular weight"]) if r.get("molecular weight") else None,
                "flavor_profile": [p.strip() for p in r.get("flavor profile", "").replace("@", ",").split(",") if p.strip()],
                "in_project_library": known is not None,
                "still_needed": [] if known is not None and known.Tm else ["Tm", "groups", "dHfus (optional)"],
            })
    return out
