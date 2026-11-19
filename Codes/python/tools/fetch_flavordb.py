"""
Download the FlavorDB2 records we need (Member 1 - Data).

Run this on a laptop with normal internet access (the IIIT-Delhi site refuses
automated access from cloud servers). The data files in data/ already hold the
scan of 2026-11-09; rerun only to refresh them:

    python tools/fetch_flavordb.py              # scans entity ids 1..1000
    python tools/fetch_flavordb.py --max-id 1200 --delay 0.5

What it does
1. Reads https://cosylab.iiitd.edu.in/flavordb2/entities_json?id=N for each N.
2. Keeps the entities whose name matches one of our spices (turmeric, black
   pepper, pepper, vanilla, ajwain/carom) and writes them, with the target
   compounds they contain, to data/flavordb_spice_entities.csv.
3. For each target compound (vanillin, curcumin, piperine, thymol) found in any
   entity, copies its molecule record into data/flavordb_records.csv
   (status = retrieved). Compounds not found keep their previous row.

Only the Python standard library is used. Be polite to the server: keep the
delay, and run it once - the CSV files are what the code reads afterwards.
"""
import argparse
import csv
import json
import os
import sys
import time
import urllib.request

BASE = "https://cosylab.iiitd.edu.in/flavordb2/entities_json?id={}"
TARGETS = {1183: "vanillin", 969516: "curcumin", 638024: "piperine", 6989: "thymol"}
SPICE_WORDS = ("turmeric", "pepper", "vanilla", "ajwain", "carom", "fenugreek")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
RECORDS = os.path.join(DATA, "flavordb_records.csv")
ENTITIES = os.path.join(DATA, "flavordb_spice_entities.csv")

FIELDS = ["compound", "pubchem_cid", "status", "flavordb_source", "cas", "fema_number",
          "fema_flavor_profile", "fooddb_flavor_profile", "odor", "taste", "functional_groups",
          "molecular_weight", "xlogp", "hbd", "hba", "fooddb_id", "note"]


def get_json(url, timeout=30):
    req = urllib.request.Request(url, headers={"User-Agent": "CL202-course-project/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8", "replace"))


def pick(rec, *names):
    """Return the first field whose key contains one of `names` (FlavorDB key names vary)."""
    for n in names:
        for k, v in rec.items():
            if n == k.lower() or n in k.lower():
                if isinstance(v, list):
                    v = "; ".join(str(x) for x in v)
                return "" if v is None else str(v)
    return ""


def molecule_row(name, m, entity_id, entity_name):
    return {
        "compound": name,
        "pubchem_cid": pick(m, "pubchem_id"),
        "status": "retrieved",
        "flavordb_source": f"FlavorDB2 entities_json?id={entity_id} ({entity_name})",
        "cas": pick(m, "cas"),
        "fema_number": pick(m, "fema_number"),
        "fema_flavor_profile": pick(m, "fema_flavor_profile"),
        "fooddb_flavor_profile": pick(m, "fooddb_flavor_profile"),
        "odor": pick(m, "odor"),
        "taste": pick(m, "taste"),
        "functional_groups": pick(m, "functional_groups").replace("@", "; "),
        "molecular_weight": pick(m, "molecular_weight"),
        "xlogp": pick(m, "xlogp"),
        "hbd": pick(m, "hbd", "num_hbd", "h_bond_donor"),
        "hba": pick(m, "hba", "num_hba", "h_bond_acceptor"),
        "fooddb_id": pick(m, "fooddb_id", "foodb_id"),
        "note": "Retrieved " + time.strftime("%Y-%m-%d") + " by tools/fetch_flavordb.py",
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--min-id", type=int, default=1)
    ap.add_argument("--max-id", type=int, default=1000)
    ap.add_argument("--delay", type=float, default=0.3, help="seconds between requests")
    a = ap.parse_args(argv)

    old = {}
    if os.path.exists(RECORDS):
        with open(RECORDS, newline="", encoding="utf-8") as f:
            old = {r["compound"]: r for r in csv.DictReader(f)}

    found, entities, seen_ids = {}, [], set()
    for eid in range(a.min_id, a.max_id + 1):
        try:
            e = get_json(BASE.format(eid))
        except Exception as exc:                      # network error or no such id
            print(f"id {eid}: {exc}", file=sys.stderr)
            time.sleep(a.delay)
            continue
        real_id = e.get("entity_id")
        if real_id in seen_ids:                       # server returns a default entity for unknown ids
            time.sleep(a.delay)
            continue
        seen_ids.add(real_id)
        name = str(e.get("entity_alias_readable", ""))
        mols = e.get("molecules", []) or []
        hits = []
        for m in mols:
            try:
                cid = int(m.get("pubchem_id"))
            except (TypeError, ValueError):
                continue
            if cid in TARGETS:
                hits.append(TARGETS[cid])
                found.setdefault(TARGETS[cid], molecule_row(TARGETS[cid], m, real_id, name))
        if any(w in name.lower() for w in SPICE_WORDS):
            entities.append({"entity_id": real_id, "name": name, "category": e.get("category", ""),
                             "scientific_name": e.get("natural_source_name", ""),
                             "n_molecules": len(mols), "target_compounds": "; ".join(sorted(set(hits)))})
            print(f"{real_id:>5}  {name:<25} {len(mols):>4} molecules  {hits}")
        time.sleep(a.delay)

    with open(ENTITIES, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["entity_id", "name", "category", "scientific_name",
                                          "n_molecules", "target_compounds"])
        w.writeheader()
        w.writerows(entities)

    rows = []
    for cid, name in TARGETS.items():
        rows.append(found.get(name) or old.get(name) or {"compound": name, "pubchem_cid": cid,
                                                         "status": "not_retrieved"})
    with open(RECORDS, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {ENTITIES} ({len(entities)} spice entities)")
    print(f"wrote {RECORDS}: " + ", ".join(f"{r['compound']}={r['status']}" for r in rows))


if __name__ == "__main__":
    main()
