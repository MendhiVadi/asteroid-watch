"""Fetch planetary / lunar ephemerides from JPL Horizons for the N-body engine.

Output (app/public/data/ephemeris/):
  planets.bin   Float64, layout [body][step][x,y,z,vx,vy,vz]  (barycentric, ecliptic J2000), 1-day steps
                bodies: sun, mercury, venus, earth, mars, jupiter, saturn, uranus, neptune
  moon.bin      Float64, layout [step][x,y,z,vx,vy,vz]        (GEOCENTRIC, ecliptic J2000), 12-hour steps
  manifest.json grids, GM values (AU^3/day^2 and km^3/s^2), radii, provenance, validation numbers

Coverage: 2020-01-01 .. 2100-12-31 (JD 2458849.5 .. 2488433.5).

Moon step: cubic-Hermite interpolation error of the geocentric Moon vs. Horizons was measured for steps of
2/4/6/8/12/16/24 h (worst case over 2020-2036: 4 m / 16 m / 63 m / 319 m / 1.0 km / 5.1 km for 4/6/8/12/16/24 h);
12 h is the coarsest divisor of one day under 1 km, and makes moon.bin smaller than the old 2-hour table.

Units: AU and AU/day, time = JD (TDB).  Frame = ICRF axes rotated to the J2000
ecliptic (what Horizons calls REF_PLANE=ECLIPTIC, REF_SYSTEM=ICRF), which is the
same frame the SBDB orbital elements are expressed in.

Masses come from the JPL DE440 header (ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de440/header.440).
Mars, Jupiter, Saturn, Uranus, Neptune use their system barycentres (GM includes moons).

If Horizons or the GM source is unreachable the script exits non-zero and writes
nothing; it never fabricates data.  Requests run in a thread pool (--workers).

Usage:
  python -I scripts/fetch_ephemeris.py                      full fetch, 2020-01-01 .. 2100-12-31
  python -I scripts/fetch_ephemeris.py --extend --validate  fetch ONLY the span missing from the existing
                                                            tables, append, re-grid the Moon, validate
  python -I scripts/fetch_ephemeris.py --validate-only      measure Hermite error of the files on disk
"""
import argparse
import json
import math
import os
import random
import re
import sys
import time
import urllib.parse
import urllib.request
from array import array
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "app" / "public" / "data" / "ephemeris"
HORIZONS = "https://ssd.jpl.nasa.gov/api/horizons.api"
DE440_HEADER = "https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de440/header.440"

# name -> Horizons COMMAND (barycentre codes for Mars+ systems)
MAIN_BODIES = [
    ("sun", "10"),
    ("mercury", "199"),
    ("venus", "299"),
    ("earth", "399"),
    ("mars", "4"),
    ("jupiter", "5"),
    ("saturn", "6"),
    ("uranus", "7"),
    ("neptune", "8"),
]

# 2020-01-01 00:00 TDB .. 2100-12-31 00:00 TDB inclusive = 29585 daily rows
DEFAULT_START_JD = 2458849.5
DEFAULT_DAYS = 29585      # number of daily points (JD 2488433.5 = 2100-12-31)
DEFAULT_MOON_HOURS = 12   # Moon node spacing

EARTH_RADIUS_KM = 6371.0
MOON_RADIUS_KM = 1737.4
SUN_RADIUS_KM = 695700.0
AU_KM_CONST = 149597870.7

PLANET_CHUNK = 3100   # rows per daily request
MOON_CHUNK = 4400     # rows per Moon request


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def http_get(url, tries=10, timeout=180):
    last = None
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001
            last = e
            # Horizons answers 503 when it is busy / throttling: back off exponentially (capped), with jitter
            time.sleep(min(60.0, 3.0 * 1.6 ** k) + random.random() * 3)
    raise RuntimeError(f"NETWORK_FAILURE: {url[:120]}... -> {last!r}")


CACHE_DIR = None  # set by --cache-dir: completed request chunks are stored here so a re-run resumes cheaply


def _cache_path(command, center, start_jd, n, step_days):
    if CACHE_DIR is None:
        return None
    key = f"{command}_{center.replace('@', 'at')}_{start_jd:.6f}_{n}_{step_days:.9f}.bin"
    return Path(CACHE_DIR) / key


def horizons_vectors(command, center, start_jd, n, step_days):
    """Return list of (jd, x,y,z,vx,vy,vz) for n points start_jd + i*step_days (cached on disk if --cache-dir)."""
    cp = _cache_path(command, center, start_jd, n, step_days)
    if cp is not None and cp.exists():
        a = array("d")
        a.frombytes(cp.read_bytes())
        if len(a) == 7 * n:
            return [tuple(a[7 * i:7 * i + 7]) for i in range(n)]
    rows = _horizons_vectors_net(command, center, start_jd, n, step_days)
    if cp is not None:
        a = array("d")
        for r in rows:
            a.extend(r)
        cp.parent.mkdir(parents=True, exist_ok=True)
        write_atomic(cp, a.tobytes())
    return rows


def _horizons_vectors_net(command, center, start_jd, n, step_days):
    stop_jd = start_jd + (n - 1) * step_days
    step_min = int(round(step_days * 1440))
    q = {
        "format": "text",
        "COMMAND": f"'{command}'",
        "OBJ_DATA": "'NO'",
        "MAKE_EPHEM": "'YES'",
        "EPHEM_TYPE": "'VECTORS'",
        "CENTER": f"'{center}'",
        "START_TIME": f"'JD {start_jd:.9f}'",
        "STOP_TIME": f"'JD {stop_jd:.9f}'",
        "STEP_SIZE": f"'{step_min} min'",
        "OUT_UNITS": "'AU-D'",
        "REF_PLANE": "'ECLIPTIC'",
        "REF_SYSTEM": "'ICRF'",
        "VEC_TABLE": "'2'",
        "CSV_FORMAT": "'YES'",
    }
    url = HORIZONS + "?" + urllib.parse.urlencode(q, safe="'@ ").replace(" ", "%20")
    txt = ""
    for attempt in range(4):  # an overloaded API can answer 200 with an error page: retry a few times
        txt = http_get(url)
        if "$$SOE" in txt:
            break
        time.sleep(3 + 4 * attempt + random.random() * 2)
    if "$$SOE" not in txt:
        raise RuntimeError(f"Horizons error for {command}@{center}: {txt[:600]}")
    body = txt.split("$$SOE")[1].split("$$EOE")[0]
    rows = []
    for line in body.strip().splitlines():
        p = [s.strip() for s in line.split(",")]
        rows.append((float(p[0]), *[float(v) for v in p[2:8]]))
    if len(rows) != n:
        raise RuntimeError(f"{command}@{center}: expected {n} rows, got {len(rows)}")
    for i, r in enumerate(rows):
        if abs(r[0] - (start_jd + i * step_days)) > 1e-6:
            raise RuntimeError(f"{command}: grid mismatch at row {i}: {r[0]} vs {start_jd + i * step_days}")
    return rows


def series_tasks(name, command, center, start_jd, n, step_days, chunk):
    """Split one series into <= chunk-row requests."""
    tasks = []
    done = 0
    while done < n:
        m = min(chunk, n - done)
        tasks.append({"series": name, "idx": len(tasks), "command": command, "center": center,
                      "start_jd": start_jd + done * step_days, "n": m, "step": step_days})
        done += m
    return tasks


def run_tasks(tasks, workers):
    """Run requests in a thread pool; returns {series: rows concatenated in time order}."""
    t0 = time.time()
    total_rows = sum(t["n"] for t in tasks)
    done_rows = 0
    parts = {}
    ex = ThreadPoolExecutor(max_workers=workers)
    try:
        futs = {ex.submit(horizons_vectors, t["command"], t["center"], t["start_jd"], t["n"], t["step"]): t for t in tasks}
        for k, f in enumerate(as_completed(futs), 1):
            t = futs[f]
            parts[(t["series"], t["idx"])] = f.result()
            done_rows += t["n"]
            el = time.time() - t0
            log(f"  {k}/{len(tasks)} requests, {done_rows}/{total_rows} rows ({100 * done_rows / total_rows:.0f}%), "
                f"{el:.0f}s elapsed  [last: {t['series']} chunk {t['idx']}]")
    except BaseException:
        ex.shutdown(wait=False, cancel_futures=True)
        raise
    ex.shutdown(wait=True)
    series = {}
    for (s, i) in sorted(parts):
        series.setdefault(s, []).extend(parts[(s, i)])
    return series


def parse_de440_header():
    txt = http_get(DE440_HEADER)
    sec = txt.split("GROUP   1041")[1].split("GROUP   1050")[0]
    vals = [float(m.replace("D", "E")) for m in re.findall(r"[-+]?\d\.\d+D[-+]\d+", sec)]
    # names in GROUP 1040 order (see header): DENUM LENUM TDATEF TDATEB JDEPOC CENTER CLIGHT BETA GAMMA AU
    # EMRAT GM1 GM2 GMB GM4 GM5 GM6 GM7 GM8 GM9 GMS ...
    names = ["DENUM", "LENUM", "TDATEF", "TDATEB", "JDEPOC", "CENTER", "CLIGHT", "BETA", "GAMMA", "AU",
             "EMRAT", "GM1", "GM2", "GMB", "GM4", "GM5", "GM6", "GM7", "GM8", "GM9", "GMS"]
    c = dict(zip(names, vals))
    if abs(c["DENUM"] - 440) > 1e-9 or abs(c["AU"] - 149597870.7) > 1e-3 or abs(c["CLIGHT"] - 299792.458) > 1e-6:
        raise RuntimeError(f"Unexpected DE440 header parse: {c}")
    return c


def write_atomic(path, data):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_bytes(data)
    os.replace(tmp, path)


def flatten(rows):
    a = array("d")
    for r in rows:
        a.extend(r[1:7])
    return a


def max_diff(rows_first, ref6):
    """Largest |difference| of the first fetched row vs. an existing 6-vector (pos AU, vel AU/d)."""
    r = rows_first
    dp = max(abs(r[1 + d] - ref6[d]) for d in range(3))
    dv = max(abs(r[4 + d] - ref6[3 + d]) for d in range(3))
    return dp, dv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start-jd", type=float, default=DEFAULT_START_JD)
    ap.add_argument("--days", type=int, default=DEFAULT_DAYS, help="total number of daily points")
    ap.add_argument("--moon-hours", type=float, default=DEFAULT_MOON_HOURS, help="Moon node spacing in hours")
    ap.add_argument("--workers", type=int, default=4, help="parallel Horizons requests (Horizons throttles with 503 above ~4-6)")
    ap.add_argument("--cache-dir", default=None, help="store completed request chunks here; a re-run after a failure resumes from them")
    ap.add_argument("--extend", action="store_true",
                    help="keep the existing tables, fetch only the missing span (and re-grid the Moon if its step changed)")
    ap.add_argument("--moon-only", action="store_true", help="re-fetch only the Moon table, keep existing planets.bin")
    ap.add_argument("--validate", action="store_true", help="after writing, measure cubic-Hermite error against fresh Horizons data")
    ap.add_argument("--validate-only", action="store_true", help="only validate the files on disk (writes nothing except the manifest 'validation' entry)")
    args = ap.parse_args()
    global CACHE_DIR
    CACHE_DIR = args.cache_dir

    OUT.mkdir(parents=True, exist_ok=True)

    if args.validate_only:
        man = json.loads((OUT / "manifest.json").read_text())
        planets = array("d")
        planets.frombytes((OUT / "planets.bin").read_bytes())
        moon = array("d")
        moon.frombytes((OUT / "moon.bin").read_bytes())
        res = validate(man["main"]["jd0"], man["main"]["count"], planets, moon, man["moon"]["stepDays"],
                       man.get("extension", {}).get("seamIndex"), args.workers)
        man["validation"] = res
        write_atomic(OUT / "manifest.json", json.dumps(man, indent=1).encode())
        return

    per_day = 24.0 / args.moon_hours
    if abs(per_day - round(per_day)) > 1e-9:
        raise RuntimeError("--moon-hours must divide 24")
    per_day = int(round(per_day))
    moon_step = 1.0 / per_day
    n = args.days
    start = args.start_jd
    n_moon = (n - 1) * per_day + 1
    workers = args.workers

    if args.extend:
        man = json.loads((OUT / "manifest.json").read_text())
        if abs(man["main"]["jd0"] - start) > 1e-9:
            raise RuntimeError("--extend: manifest jd0 differs from --start-jd")
        old_n = man["main"]["count"]
        old_planets = array("d")
        old_planets.frombytes((OUT / "planets.bin").read_bytes())
        if len(old_planets) != len(MAIN_BODIES) * old_n * 6:
            raise RuntimeError("--extend: planets.bin does not match manifest")
        old_moon = array("d")
        old_moon.frombytes((OUT / "moon.bin").read_bytes())
        old_moon_step = man["moon"]["stepDays"]
        old_moon_n = man["moon"]["count"]
        if len(old_moon) != old_moon_n * 6:
            raise RuntimeError("--extend: moon.bin does not match manifest")
        if old_n >= n:
            raise RuntimeError(f"--extend: tables already cover {old_n} days (target {n}); nothing to fetch")
        k = round(moon_step / old_moon_step)
        if k < 1 or abs(k * old_moon_step - moon_step) > 1e-9 or (old_moon_n - 1) % k != 0:
            raise RuntimeError(f"--extend: cannot re-grid Moon from step {old_moon_step} d to {moon_step} d")
        # Re-grid the existing Moon table (pure subsampling: 12 h nodes are every 6th 2-hour node)
        kept_moon = array("d")
        for i in range(0, old_moon_n, k):
            kept_moon.extend(old_moon[i * 6:i * 6 + 6])
        kept_moon_n = len(kept_moon) // 6
        old_end_jd = start + (old_n - 1)
        # Fetch from the LAST existing row (an overlap row that must match what is on disk), then drop it.
        n_new_planet = n - old_n + 1
        n_new_moon = n_moon - kept_moon_n + 1
        log(f"extend: planets {old_n} -> {n} daily rows (fetch {n_new_planet - 1} new + 1 overlap per body); "
            f"moon {old_moon_n} rows @{old_moon_step * 24:g}h -> {kept_moon_n} kept @{moon_step * 24:g}h + {n_new_moon - 1} new")
        tasks = []
        for name, cmd in MAIN_BODIES:
            tasks += series_tasks(name, cmd, "500@0", old_end_jd, n_new_planet, 1.0, PLANET_CHUNK)
        tasks += series_tasks("moon", "301", "500@399", old_end_jd, n_new_moon, moon_step, MOON_CHUNK)
        log(f"{len(tasks)} Horizons requests, {workers} workers")
        res = run_tasks(tasks, workers)

        # overlap check + build planets [body][step]
        planets = array("d")
        worst_p = worst_v = 0.0
        for bi, (name, _) in enumerate(MAIN_BODIES):
            rows = res[name]
            ref = old_planets[(bi * old_n + old_n - 1) * 6:(bi * old_n + old_n) * 6]
            dp, dv = max_diff(rows[0], ref)
            worst_p, worst_v = max(worst_p, dp), max(worst_v, dv)
            if dp > 1e-9 or dv > 1e-10:
                raise RuntimeError(f"{name}: new data does not match existing table at overlap JD {old_end_jd}: dpos {dp:.3e} AU, dvel {dv:.3e} AU/d")
            planets.extend(old_planets[bi * old_n * 6:(bi + 1) * old_n * 6])
            planets.extend(flatten(rows[1:]))
        mrows = res["moon"]
        dp, dv = max_diff(mrows[0], kept_moon[-6:])
        worst_p, worst_v = max(worst_p, dp), max(worst_v, dv)
        if dp > 1e-9 or dv > 1e-10:
            raise RuntimeError(f"moon: new data does not match existing table at overlap: dpos {dp:.3e} AU, dvel {dv:.3e} AU/d")
        moon = kept_moon
        moon.extend(flatten(mrows[1:]))
        log(f"seam check OK: first fetched row equals last existing row for all 10 series (max dpos {worst_p:.2e} AU, dvel {worst_v:.2e} AU/d)")
        base_manifest = man
        seam_index = old_n - 1
    else:
        log("fetching DE440 constants ...")
        c = parse_de440_header()
        base_manifest = None
        seam_index = None
        tasks = []
        planets = array("d")
        if args.moon_only:
            planets.frombytes((OUT / "planets.bin").read_bytes())
            if len(planets) != len(MAIN_BODIES) * n * 6:
                raise RuntimeError("--moon-only: existing planets.bin does not match --start-jd/--days")
        else:
            for name, cmd in MAIN_BODIES:
                tasks += series_tasks(name, cmd, "500@0", start, n, 1.0, PLANET_CHUNK)
        tasks += series_tasks("moon", "301", "500@399", start, n_moon, moon_step, MOON_CHUNK)
        log(f"{len(tasks)} Horizons requests, {workers} workers")
        res = run_tasks(tasks, workers)
        if not args.moon_only:
            for name, _ in MAIN_BODIES:
                planets.extend(flatten(res[name]))
        moon = flatten(res["moon"])

    if len(planets) != len(MAIN_BODIES) * n * 6 or len(moon) != n_moon * 6:
        raise RuntimeError("internal size mismatch")

    write_atomic(OUT / "planets.bin", planets.tobytes())
    write_atomic(OUT / "moon.bin", moon.tobytes())

    if base_manifest is not None:
        manifest = base_manifest
        manifest["extension"] = {
            "extendedUtc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "baseEndJd": start + seam_index,
            "seamIndex": seam_index,
            "note": "rows after seamIndex fetched in a second pass; first fetched row verified equal to the last existing row",
        }
    else:
        au_km = c["AU"]
        day_s = 86400.0
        k3 = au_km ** 3 / day_s ** 2  # AU^3/d^2 -> km^3/s^2
        emrat = c["EMRAT"]
        gm_emb = c["GMB"]
        gm = {
            "sun": c["GMS"],
            "mercury": c["GM1"],
            "venus": c["GM2"],
            "earth": gm_emb * emrat / (1 + emrat),
            "moon": gm_emb / (1 + emrat),
            "mars": c["GM4"],
            "jupiter": c["GM5"],
            "saturn": c["GM6"],
            "uranus": c["GM7"],
            "neptune": c["GM8"],
        }
        manifest = {
            "version": 1,
            "source": "JPL Horizons VECTORS (DE441); constants from JPL DE440 header",
            "fetchedUtc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "frame": "ICRF axes, J2000 ecliptic plane; barycentric (SSB) unless noted",
            "units": {"length": "AU", "velocity": "AU/day", "time": "JD (TDB)", "gm": "AU^3/day^2"},
            "auKm": au_km,
            "cKmS": c["CLIGHT"],
            "cAuDay": c["CLIGHT"] * day_s / au_km,
            "emrat": emrat,
            "main": {
                "file": "planets.bin",
                "dtype": "float64",
                "layout": "[body][step][x,y,z,vx,vy,vz]",
                "jd0": start,
                "stepDays": 1.0,
                "count": n,
                "bodies": [b for b, _ in MAIN_BODIES],
            },
            "moon": {},
            "gm": {k: {"au3d2": v, "km3s2": v * k3} for k, v in gm.items()},
            "radiiKm": {"earth": EARTH_RADIUS_KM, "moon": MOON_RADIUS_KM, "sun": SUN_RADIUS_KM},
        }

    manifest["main"]["count"] = n
    manifest["moon"] = {
        "file": "moon.bin",
        "dtype": "float64",
        "layout": "[step][x,y,z,vx,vy,vz]",
        "relativeTo": "earth",
        "jd0": start,
        "stepDays": moon_step,
        "count": n_moon,
    }
    manifest["jdStart"] = start
    manifest["jdEnd"] = start + (n - 1)
    manifest["notes"] = [
        "mars/jupiter/saturn/uranus/neptune are system barycentres; their GM includes satellites.",
        f"earth and moon are separate bodies; moon.bin is stored relative to earth, Float64, {args.moon_hours:g}-hour nodes "
        "(cubic Hermite on position+velocity; step chosen as the coarsest divisor of 1 day with Moon position error < 1 km, "
        "see 'validation'; Float32 gave interpolation jerk that limits step size near the Moon).",
        "No minor-planet perturbers (Ceres, Vesta, Pallas ...) and no Pluto are included.",
    ]
    write_atomic(OUT / "manifest.json", json.dumps(manifest, indent=1).encode())
    log(f"wrote {OUT}: planets.bin {len(planets) * 8} bytes; moon.bin {len(moon) * 8} bytes; "
        f"coverage JD {manifest['jdStart']} .. {manifest['jdEnd']}")

    if args.validate:
        res = validate(start, n, planets, moon, moon_step, seam_index, workers)
        manifest["validation"] = res
        write_atomic(OUT / "manifest.json", json.dumps(manifest, indent=1).encode())
        log("manifest updated with validation numbers")


def hermite(p0, v0, p1, v1, h, s):
    s2, s3 = s * s, s * s * s
    return ((2 * s3 - 3 * s2 + 1) * p0 + (s3 - 2 * s2 + s) * h * v0
            + (-2 * s3 + 3 * s2) * p1 + (s3 - s2) * h * v1)


def _err(table, base_idx, base1_idx, s, h, row):
    e2 = 0.0
    for d in range(3):
        val = hermite(table[base_idx + d], table[base_idx + 3 + d], table[base1_idx + d], table[base1_idx + 3 + d], h, s)
        e2 += (val - row[1 + d]) ** 2
    return e2 ** 0.5 * AU_KM_CONST


def validate(start, n, planets, moon, moon_step, seam_index, workers):
    """Cubic-Hermite error of the stored tables vs. fresh fine-grained Horizons data (parallel requests).

    Earth and Mercury (worst planet) : 1-day nodes vs. 6-hourly truth, 240-day windows.
    Moon (geocentric)                : stored nodes vs. 2-hourly truth, 60-day windows.
    Windows: early span, around the seam of the extension (if any), middle and end of the coverage.
    """
    win_p, win_m = 240, 60
    p_offs = {}
    m_offs = {}
    p_offs["start"] = 2000
    m_offs["start"] = 2000
    if seam_index:
        p_offs["seam"] = seam_index - win_p // 2
        m_offs["seam"] = seam_index - win_m // 2
    p_offs["mid"] = int(n * 0.6)
    m_offs["mid"] = int(n * 0.5)
    m_offs["late"] = int(n * 0.8)
    p_offs["end"] = n - 1 - win_p - 1
    m_offs["end"] = n - 1 - win_m - 1

    tasks = []
    idx_of = {name: i for i, (name, _) in enumerate(MAIN_BODIES)}
    cmd_of = dict(MAIN_BODIES)
    for body in ("earth", "mercury"):
        for w, off in p_offs.items():
            tasks.append({"series": f"{body}/{w}", "idx": 0, "command": cmd_of[body], "center": "500@0",
                          "start_jd": start + off, "n": win_p * 4, "step": 0.25})
    for w, off in m_offs.items():
        tasks.append({"series": f"moon/{w}", "idx": 0, "command": "301", "center": "500@399",
                      "start_jd": start + off, "n": win_m * 12, "step": 1.0 / 12})
    log(f"validation: {len(tasks)} requests ...")
    res = run_tasks(tasks, workers)

    out = {"method": "cubic Hermite on stored nodes (position+velocity) vs fresh Horizons samples; error = |dr| in km",
           "planets": {}, "moon": {}}
    for body in ("earth", "mercury"):
        bi = idx_of[body]
        worst = 0.0
        per = {}
        for w in p_offs:
            rows = res[f"{body}/{w}"]
            wmax = 0.0
            for r in rows:
                t = r[0] - start
                i = int(math.floor(t + 1e-9))
                s = t - i
                b0 = (bi * n + i) * 6
                b1 = (bi * n + i + 1) * 6
                wmax = max(wmax, _err(planets, b0, b1, s, 1.0, r))
            per[w] = round(wmax, 4)
            worst = max(worst, wmax)
        out["planets"][body] = {"nodeStepDays": 1.0, "worstKm": round(worst, 4), "perWindowKm": per}
        log(f"  {body} 1-day Hermite: worst {worst * 1000:.1f} m  per window (km): {per}")
    mworst = 0.0
    mrms_acc = 0.0
    mcount = 0
    per = {}
    for w in m_offs:
        rows = res[f"moon/{w}"]
        wmax = 0.0
        for r in rows:
            u = (r[0] - start) / moon_step
            i = int(math.floor(u + 1e-7))
            s = u - i
            e = _err(moon, i * 6, (i + 1) * 6, s, moon_step, r)
            wmax = max(wmax, e)
            mrms_acc += e * e
            mcount += 1
        per[w] = round(wmax, 4)
        mworst = max(mworst, wmax)
    out["moon"] = {"nodeStepDays": moon_step, "worstKm": round(mworst, 4),
                   "rmsKm": round((mrms_acc / mcount) ** 0.5, 4), "perWindowKm": per,
                   "samples": mcount}
    log(f"  moon {moon_step * 24:g}-hour Hermite: worst {mworst * 1000:.1f} m, rms {out['moon']['rmsKm'] * 1000:.1f} m  per window (km): {per}")
    return out


if __name__ == "__main__":
    try:
        main()
    except RuntimeError as e:
        print(str(e), file=sys.stderr, flush=True)
        # hard exit: do not wait for in-flight worker threads that are still retrying against a throttled API
        os._exit(2)
