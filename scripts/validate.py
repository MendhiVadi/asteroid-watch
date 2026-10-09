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
# ---- derived outputs (app/public/data/asteroids.json, data/asteroids_trajectories.json) ----------
# Analysis window 2026-10-09 .. 2100-12-31; see scripts/trajectory.py and data/trajectory_report.md.
LO_ISO, HI_ISO = "2026-10-09T00:00:00Z", "2100-12-31T00:00:00Z"
pub = ROOT / "app" / "public" / "data" / "asteroids.json"
if pub.exists():
    rows = json.load(open(pub))
    mb = pub.stat().st_size / 1e6
    print(f"INFO asteroids.json {len(rows):,} rows, {mb:.2f} MB")
    check(len(rows) == len(data) and [r["id"] for r in rows] == ids, "asteroids.json rows == neo_elements rows (same ids, same order)")
    check(mb < 16.0, f"asteroids.json size {mb:.2f} MB < 16 MB")
    contract = ["id", "name", "diameter_km", "pha", "category", "min_dist_au", "min_date", "a", "e", "i", "om", "w", "ma", "epoch",
                "strict_collision", "trend_au", "trend_au_10y", "trend_au_full"]
    check(all(all(r.get(k) is not None for k in contract) for r in rows), "asteroids.json: contract + trend fields present and non-null")
    nums = ["diameter_km", "min_dist_au", "a", "e", "i", "om", "w", "ma", "epoch", "trend_au", "trend_au_10y", "trend_au_full"]
    check(all(math.isfinite(r[k]) for r in rows for k in nums), "asteroids.json: no NaN/inf")
    check(all(r["category"] in ("impact", "approaching", "receding") for r in rows), "asteroids.json: category in {impact, approaching, receding}")
    check(all(LO_ISO <= r["min_date"] <= HI_ISO for r in rows), f"asteroids.json: min_date inside {LO_ISO[:10]} .. {HI_ISO[:10]}")
    ca_ok = True
    for r in rows:
        lst = r.get("next_close_approach")
        if lst is None:
            continue
        d_ = [c["date"] for c in lst]
        ca_ok &= 1 <= len(lst) <= 3 and d_ == sorted(d_) and all(0 <= c["dist_au"] < 0.05 + 1e-6 and LO_ISO < c["date"] < HI_ISO for c in lst)
        ca_ok &= ("n_close_approaches" not in r) or (r["n_close_approaches"] > len(lst))
    check(ca_ok, "asteroids.json: next_close_approach lists well-formed (<= 3, chronological, < 0.05 AU, inside window)")
    cnt = {c: sum(r["category"] == c for r in rows) for c in ("impact", "approaching", "receding")}
    print(f"INFO categories {cnt}, strict_collision {sum(bool(r['strict_collision']) for r in rows)}, "
          f"with next_close_approach {sum('next_close_approach' in r for r in rows):,}")
    trj = ROOT / "data" / "asteroids_trajectories.json"
    if trj.exists():
        tr = json.load(open(trj))
        check(len(tr["objects"]) == 300 and all(len(o[k]) == 180 for o in tr["objects"] for k in ("helio_orbit", "geo_window", "geo_approach"))
              and len(tr["earth_orbit"]) == 180 and tr["meta"]["end_iso"] == HI_ISO,
              "asteroids_trajectories.json: 300 objects x 60-point (x3) polylines, window ends 2100-12-31")
else:
    print("INFO asteroids.json not generated yet (run scripts/trajectory.py)")


print("\nRESULT:", "FAILED" if fails else "ALL CHECKS PASSED")
sys.exit(1 if fails else 0)
