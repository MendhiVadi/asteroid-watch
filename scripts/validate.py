"""Validate data/neo_elements.json.  Run: python -I scripts/validate.py  (exit 1 on failure)"""
import json, math, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
data = json.load(open(ROOT / "data" / "neo_elements.json"))
raw = json.load(open(ROOT / "data" / "raw" / "sbdb_neo.json"))
fails = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        fails.append(msg)


REQ = ["a", "e", "i", "om", "w", "ma", "epoch"]
NUM = REQ + ["per", "moid", "H", "diameter_km"]
KEYS = ["id", "name", "a", "e", "i", "om", "w", "ma", "epoch", "per", "moid", "H", "diameter_km", "pha",
        "in_csv", "hazardous_csv", "min_miss_distance_km_csv"]

check(isinstance(data, list) and len(data) > 0, f"array non-empty ({len(data):,})")
check(all(all(k in o for k in KEYS) for o in data), "all required keys present")
bad_nan = [o["id"] for o in data for k in NUM if o[k] is not None and not math.isfinite(o[k])]
check(not bad_nan, f"no NaN/inf numeric values ({len(bad_nan)} bad)")
bad_null = [o["id"] for o in data if any(o[k] is None for k in REQ)]
check(not bad_null, f"no null orbital elements ({len(bad_null)} bad)")
check(all(o["e"] < 1 and o["e"] >= 0 for o in data), "0 <= e < 1 for all")
check(all(o["a"] > 0 for o in data), "a > 0 for all")
check(all(0 <= o["i"] <= 180 for o in data), "0 <= i <= 180")
check(all(0 <= o[k] <= 360 for o in data for k in ("om", "w", "ma")), "om, w, ma within [0,360]")
check(all(o["epoch"] > 2400000 for o in data), "epoch is a Julian Date")
ids = [o["id"] for o in data]
check(len(ids) == len(set(ids)), "ids unique")
csv_ids = [o["csv_id"] for o in data if o["csv_id"]]
check(len(csv_ids) == len(set(csv_ids)), "csv_id unique among in_csv objects")
check(all((o["csv_id"] is not None) == o["in_csv"] for o in data), "in_csv consistent with csv_id")
check(all((o["hazardous_csv"] is not None) == o["in_csv"] for o in data), "hazardous_csv set iff in_csv")
check(all((o["min_miss_distance_km_csv"] is not None) == o["in_csv"] for o in data),
      "min_miss_distance_km_csv set iff in_csv")
nodiam = [o["id"] for o in data if o["diameter_km"] is None or o["diameter_km"] <= 0]
print(f"INFO objects without usable diameter: {len(nodiam)}")
check(len(nodiam) <= 10, "<=10 objects lack diameter (no H, no CSV)")
# independent recount of the build filter on the raw SBDB rows (null element, e>=1, a<=0 are dropped)
_fi = {k: raw["fields"].index(k) for k in ("a", "e", "i", "om", "w", "ma", "epoch")}


def _keep(row):
    try:
        v = {k: float(row[j]) for k, j in _fi.items()}
    except (TypeError, ValueError):
        return False
    return v["e"] < 1 and v["a"] > 0


check(raw["count"] == len(raw["data"]), f"SBDB count field {raw['count']:,} == rows {len(raw['data']):,}")
expected_keep = sum(1 for r in raw["data"] if _keep(r))
check(len(data) == expected_keep, f"neo_elements rows {len(data):,} == raw rows passing filter {expected_keep:,}")
dropped = raw["count"] - len(data)
print(f"INFO dropped vs SBDB: {dropped}")
n_csv = sum(o["in_csv"] for o in data)
check(27000 <= n_csv <= 27423, f"in_csv count plausible ({n_csv:,})")
print(f"INFO total={len(data):,} in_csv={n_csv:,} not_in_csv={len(data)-n_csv:,} pha={sum(o['pha'] for o in data):,}")
print("\nRESULT:", "FAILED" if fails else "ALL CHECKS PASSED")
sys.exit(1 if fails else 0)
