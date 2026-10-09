"""Fetch planetary / lunar ephemerides from JPL Horizons for the N-body engine.

Output (app/public/data/ephemeris/):
  planets.bin   Float64, layout [body][step][x,y,z,vx,vy,vz]  (barycentric, ecliptic J2000)
                bodies: sun, mercury, venus, earth, mars, jupiter, saturn, uranus, neptune
  moon.bin      Float64, layout [step][x,y,z,vx,vy,vz]        (GEOCENTRIC, ecliptic J2000)
  manifest.json grids, GM values (AU^3/day^2 and km^3/s^2), radii, provenance

Units: AU and AU/day, time = JD (TDB).  Frame = ICRF axes rotated to the J2000
ecliptic (what Horizons calls REF_PLANE=ECLIPTIC, REF_SYSTEM=ICRF), which is the
same frame the SBDB orbital elements are expressed in.

Masses come from the JPL DE440 header (ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de440/header.440).
Mars, Jupiter, Saturn, Uranus, Neptune use their system barycentres (GM includes moons).

If Horizons or the GM source is unreachable the script exits non-zero and writes
nothing; it never fabricates data.

Usage:  python -I scripts/fetch_ephemeris.py [--start-jd 2458849.5] [--days 6210] [--validate]
"""
import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from array import array
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

# 2020-01-01 00:00 TDB .. 2036-12-31 00:00 TDB inclusive = 6210 daily rows
DEFAULT_START_JD = 2458849.5
DEFAULT_DAYS = 6210  # number of daily points
MOON_PER_DAY = 12    # 2-hour steps

EARTH_RADIUS_KM = 6371.0
MOON_RADIUS_KM = 1737.4
SUN_RADIUS_KM = 695700.0


def http_get(url, tries=5, timeout=180):
    last = None
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 + 3 * k)
    raise RuntimeError(f"NETWORK_FAILURE: {url[:120]}... -> {last!r}")


def horizons_vectors(command, center, start_jd, n, step_days):
    """Return list of (jd, x,y,z,vx,vy,vz) for n points start_jd + i*step_days."""
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
    txt = http_get(url)
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


def fetch_series(command, center, start_jd, n, step_days, chunk):
    rows = []
    done = 0
    while done < n:
        m = min(chunk, n - done)
        rows += horizons_vectors(command, center, start_jd + done * step_days, m, step_days)
        done += m
        print(f"  {command}@{center}: {done}/{n}", flush=True)
    return rows


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start-jd", type=float, default=DEFAULT_START_JD)
    ap.add_argument("--days", type=int, default=DEFAULT_DAYS)
    ap.add_argument("--moon-only", action="store_true", help="re-fetch only the Moon table, keep existing planets.bin")
    ap.add_argument("--validate", action="store_true", help="also measure cubic-Hermite interpolation error (extra requests)")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    print("fetching DE440 constants ...")
    c = parse_de440_header()
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

    n = args.days
    start = args.start_jd
    planets = array("d")
    if args.moon_only:
        planets.frombytes((OUT / "planets.bin").read_bytes())
        if len(planets) != len(MAIN_BODIES) * n * 6:
            raise RuntimeError("--moon-only: existing planets.bin does not match --start-jd/--days")
    else:
        for name, cmd in MAIN_BODIES:
            print(f"{name} ({cmd}) daily ...", flush=True)
            rows = fetch_series(cmd, "500@0", start, n, 1.0, 3100)
            for r in rows:
                planets.extend(r[1:7])

    n_moon = (n - 1) * MOON_PER_DAY + 1
    print("moon geocentric 2-hourly ...", flush=True)
    mrows = fetch_series("301", "500@399", start, n_moon, 1.0 / MOON_PER_DAY, 4400)
    moon = array("d")
    for r in mrows:
        moon.extend(r[1:7])

    if not args.moon_only:
        (OUT / "planets.bin").write_bytes(planets.tobytes())
    (OUT / "moon.bin").write_bytes(moon.tobytes())

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
        "moon": {
            "file": "moon.bin",
            "dtype": "float64",
            "layout": "[step][x,y,z,vx,vy,vz]",
            "relativeTo": "earth",
            "jd0": start,
            "stepDays": 1.0 / MOON_PER_DAY,
            "count": n_moon,
        },
        "gm": {k: {"au3d2": v, "km3s2": v * k3} for k, v in gm.items()},
        "radiiKm": {"earth": EARTH_RADIUS_KM, "moon": MOON_RADIUS_KM, "sun": SUN_RADIUS_KM},
        "notes": [
            "mars/jupiter/saturn/uranus/neptune are system barycentres; their GM includes satellites.",
            "earth and moon are separate bodies; moon.bin is stored relative to earth (Float64; Float32 gave interpolation jerk that limits step size near the Moon).",
            "No minor-planet perturbers (Ceres, Vesta, Pallas ...) and no Pluto are included.",
        ],
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=1))
    print("wrote", OUT, "planets.bin", len(planets) * 8, "bytes; moon.bin", len(moon) * 8, "bytes")

    if args.validate:
        validate_hermite(start, n, planets)


def hermite(p0, v0, p1, v1, h, s):
    s2, s3 = s * s, s * s * s
    return ((2 * s3 - 3 * s2 + 1) * p0 + (s3 - 2 * s2 + s) * h * v0
            + (-2 * s3 + 3 * s2) * p1 + (s3 - s2) * h * v1)


def validate_hermite(start, n, planets):
    """Compare cubic-Hermite midpoints of Earth against Horizons 6-hourly data for 2 years."""
    idx = [b for b, _ in MAIN_BODIES].index("earth")
    yrs = 730
    rows = fetch_series("399", "500@0", start + 2000.0, yrs * 4, 0.25, 2920)
    worst = 0.0
    for r in rows:
        t = r[0] - start
        i = int(t)
        s = t - i
        base = (idx * n + i) * 6
        base1 = (idx * n + i + 1) * 6
        err = 0.0
        for d in range(3):
            val = hermite(planets[base + d], planets[base + 3 + d], planets[base1 + d], planets[base1 + 3 + d], 1.0, s)
            err += (val - r[1 + d]) ** 2
        worst = max(worst, err ** 0.5)
    print(f"Earth cubic-Hermite (1-day nodes) worst error over {yrs} d: {worst:.3e} AU = {worst * 149597870.7:.3f} km")


if __name__ == "__main__":
    try:
        main()
    except RuntimeError as e:
        print(str(e), file=sys.stderr)
        sys.exit(2)
