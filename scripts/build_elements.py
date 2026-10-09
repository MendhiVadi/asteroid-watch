"""Match neo_v2.csv against the JPL SBDB NEO catalogue; write data/neo_elements.json,
data/coverage_report.md and data/raw/match_log.json.  Run: python -I scripts/build_elements.py"""
import csv, json, math, re, time, urllib.parse, urllib.request
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "sbdb_neo.json"
ALIAS_CACHE = ROOT / "data" / "raw" / "sbdb_alias_resolve.json"
CSV = ROOT / "neo_v2.csv"
OUT = ROOT / "data" / "neo_elements.json"
REPORT = ROOT / "data" / "coverage_report.md"
ALBEDO = 0.14


def f(x):
    if x is None or x == "":
        return None
    try:
        v = float(x)
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def h_to_d(H, p=ALBEDO):
    return 1329.0 / math.sqrt(p) * 10 ** (-H / 5.0)


# ---------- load ----------
sb = json.load(open(RAW))
fields = sb["fields"]
recs = [dict(zip(fields, r)) for r in sb["data"]]
by_spk = {str(r["spkid"]): r for r in recs}
by_pdes = {r["pdes"].strip(): r for r in recs}
by_fulldes = {}
for r in recs:
    m = re.search(r"\(([^)]*)\)\s*$", r["full_name"])
    if m:
        by_fulldes[m.group(1).strip()] = r

csv_obj = {}  # id -> aggregate
with open(CSV, newline="") as fh:
    for row in csv.DictReader(fh):
        o = csv_obj.setdefault(row["id"], dict(name=row["name"].strip(), haz=False, mind=math.inf,
                                               dmin=[], dmax=[], n=0))
        o["haz"] |= row["hazardous"].strip().lower() == "true"
        o["mind"] = min(o["mind"], float(row["miss_distance"]))
        o["dmin"].append(float(row["est_diameter_min"]))
        o["dmax"].append(float(row["est_diameter_max"]))
        o["n"] += 1

# ---------- matching ----------
alias_cache = json.load(open(ALIAS_CACHE)) if ALIAS_CACHE.exists() else {}
alias_cache = {k: v for k, v in alias_cache.items() if not str(v.get("msg") or "").startswith("ERR")}


def alias_key(n):
    m = re.match(r"^(\d+)", n.strip())
    if m:
        return m.group(1)
    return n.strip().strip("()") if n.strip().startswith("(") else n.strip()


def resolve_alias(name):
    """ask SBDB single-object API for the current spkid of this name"""
    key = alias_key(name)
    if key not in alias_cache:
        try:
            with urllib.request.urlopen(
                    "https://ssd-api.jpl.nasa.gov/sbdb.api?sstr=" + urllib.parse.quote(key), timeout=60) as resp:
                r = json.load(resp)
            o = r.get("object", {})
            alias_cache[key] = dict(spkid=o.get("spkid"), neo=o.get("neo"), fullname=o.get("fullname"),
                                    msg=r.get("message"))
        except Exception as e:
            alias_cache[key] = dict(spkid=None, neo=None, fullname=None, msg="ERR " + repr(e))
        time.sleep(0.2)
    return alias_cache[key]


match = {}  # csv id -> sbdb record
how = defaultdict(int)
unmatched = []
for cid, o in csv_obj.items():
    n = o["name"]
    rec = None
    h = None
    if cid in by_spk:
        rec, h = by_spk[cid], "spkid"
    else:
        mm = re.match(r"^(\d+)\b", n)
        paren = re.findall(r"\(([^)]*)\)", n)
        d = paren[-1].strip() if paren else None
        if mm and mm.group(1) in by_pdes:
            rec, h = by_pdes[mm.group(1)], "number/pdes"
        elif d and d in by_pdes:
            rec, h = by_pdes[d], "designation/pdes"
        elif d and d in by_fulldes:
            rec, h = by_fulldes[d], "designation/full_name"
    if rec is None:
        a = resolve_alias(n)
        if a["spkid"] and str(a["spkid"]) in by_spk:
            rec, h = by_spk[str(a["spkid"])], "SBDB alias lookup"
    if rec is None:
        unmatched.append(cid)
        how["unmatched"] += 1
    else:
        match[cid] = rec
        how[h] += 1
with open(ALIAS_CACHE, "w") as fh:
    json.dump(alias_cache, fh, indent=1)

unm_info = []
for cid in unmatched:
    a = alias_cache.get(alias_key(csv_obj[cid]["name"]), {})
    unm_info.append(dict(csv_id=cid, name=csv_obj[cid]["name"], sbdb_neo=a.get("neo"),
                         sbdb_fullname=a.get("fullname")))
n_not_neo = sum(1 for u in unm_info if u["sbdb_neo"] is False)

grp = defaultdict(list)
for cid, rec in match.items():
    grp[str(rec["spkid"])].append(cid)
multi = {k: v for k, v in grp.items() if len(v) > 1}

# ---------- build elements ----------
drops = defaultdict(int)
diam_src = defaultdict(int)
out = []
for rec in recs:
    sid = str(rec["spkid"])
    el = {k: f(rec[k]) for k in ("a", "e", "i", "om", "w", "ma", "epoch", "per", "moid", "H", "diameter", "albedo")}
    if any(el[k] is None for k in ("a", "e", "i", "om", "w", "ma", "epoch")):
        drops["null_element"] += 1
        continue
    if el["e"] >= 1:
        drops["e>=1 (hyperbolic/parabolic)"] += 1
        continue
    if el["a"] <= 0:
        drops["a<=0"] += 1
        continue
    cids = grp.get(sid, [])
    D = el["diameter"]
    if D is not None and D > 0:
        src = "sbdb"
    elif el["H"] is not None:
        D, src = h_to_d(el["H"]), "H_albedo0.14"
    elif cids:
        D = sum((sum(csv_obj[c]["dmin"]) / len(csv_obj[c]["dmin"]) + sum(csv_obj[c]["dmax"]) / len(csv_obj[c]["dmax"])) / 2
                for c in cids) / len(cids)
        src = "csv_mean"
    else:
        D, src = None, "none"
    diam_src[src] += 1
    out.append(dict(
        id=sid, name=rec["full_name"].strip(), csv_id=cids[0] if cids else None,
        a=el["a"], e=el["e"], i=el["i"], om=el["om"], w=el["w"], ma=el["ma"], epoch=el["epoch"],
        per=el["per"], moid=el["moid"], H=el["H"],
        diameter_km=D, diameter_src=src, pha=rec["pha"] == "Y",
        in_csv=bool(cids),
        hazardous_csv=(any(csv_obj[c]["haz"] for c in cids) if cids else None),
        min_miss_distance_km_csv=(min(csv_obj[c]["mind"] for c in cids) if cids else None)))
with open(OUT, "w") as fh:
    json.dump(out, fh, separators=(",", ":"))

# ---------- report ----------
n_csv = len(csv_obj)
n_sb = len(recs)
matched_csv = len(match)
sb_matched = len(grp)
sb_missing = n_sb - sb_matched
pha_null = sum(1 for r in recs if r["pha"] is None)
miss_in_csv_pha = sum(1 for r in recs if str(r["spkid"]) not in grp and r["pha"] == "Y")
dropped_total = sum(drops.values())
in_csv_out = sum(1 for o in out if o["in_csv"])
lines = [
    "# NEO coverage report", "",
    f"Source catalogue: JPL SBDB Query API (`sb-kind=a`, `sb-group=neo`), fetched {time.strftime('%Y-%m-%d')}; raw: `data/raw/sbdb_neo.json`.", "",
    "## Headline", "",
    f"- SBDB total known NEO asteroids: **{n_sb:,}**",
    f"- CSV (`neo_v2.csv`) rows: 90,836; unique (id,name) objects: **{n_csv:,}**",
    f"- CSV objects matched to an SBDB NEO: **{matched_csv:,}** ({100*matched_csv/n_csv:.1f}%) -> {sb_matched:,} distinct SBDB objects "
    f"({len(multi)} SBDB objects appear under >1 CSV id, i.e. duplicated/renumbered designations)",
    f"- CSV objects NOT found in SBDB NEO list: **{len(unmatched)}** ({n_not_neo} confirmed by SBDB as no longer NEO-class, "
    f"{len(unmatched)-n_not_neo} not resolvable)",
    f"- SBDB NEOs missing from CSV: **{sb_missing:,}** ({100*sb_missing/n_sb:.1f}% of all NEOs); of these {miss_in_csv_pha} are flagged PHA",
    f"- **The CSV does NOT contain all asteroids**: it covers {100*sb_matched/n_sb:.1f}% of today's NEO catalogue "
    "(it is a ~2014-2022 close-approach set; later discoveries and objects without a listed approach are absent).", "",
    "## Matching method (CSV id -> SBDB)", "",
    "CSV ids are old-style SPK ids (2,000,000+number) while SBDB uses 20,000,000+number, so direct spkid matching only hits some. Order tried:", ""]
lines += [f"- {k}: {v:,}" for k, v in sorted(how.items(), key=lambda kv: -kv[1])]
lines += ["", "## Element table (`data/neo_elements.json`)", "",
          f"- SBDB NEOs in: {n_sb:,}",
          f"- Dropped: {dropped_total} " + ("(" + "; ".join(f"{k}: {v}" for k, v in drops.items()) + ")" if drops else ""),
          f"- Kept: **{len(out):,}** ({in_csv_out:,} in CSV, {len(out)-in_csv_out:,} not in CSV)",
          f"- PHA flag null in SBDB (treated as false): {pha_null}",
          "- Diameter source: " + ", ".join(f"{k}={v:,}" for k, v in diam_src.items()),
          "- Units: a in AU, angles in degrees, epoch as JD (TDB), per in days, moid in AU, diameter_km in km.", "",
          "## CSV objects unmatched (first 30)", ""]
lines += [f"- {u['csv_id']} {u['name']} (SBDB neo flag: {u['sbdb_neo']}, SBDB: {u['sbdb_fullname']})" for u in unm_info[:30]]
REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
with open(ROOT / "data" / "raw" / "match_log.json", "w") as fh:
    json.dump(dict(unmatched=unm_info, multi_csv_ids=multi, drops=drops), fh, indent=1)
print("\n".join(lines[:32]))
