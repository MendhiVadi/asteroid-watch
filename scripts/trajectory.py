"""Trajectory computation + categorisation of all NEOs.

Run:  python -I scripts/trajectory.py          (reads data/neo_elements.json)

Method
------
* Two-body (Keplerian) propagation of every asteroid from its catalogue osculating
  elements (ecliptic J2000; a [AU], e, i/Omega/omega/M [deg], epoch JD TDB).  Mean motion
  is n = k / a^1.5 (Gaussian gravitational constant k = 0.01720209895 rad/day); the
  catalogue `per` field is NOT used (SBDB rounds it to 3 significant digits).
* Kepler's equation M = E - e sin E solved with vectorised Newton-Raphson
  (Danby starting guess, tolerance 1e-13 rad).
* Earth (preferred): JPL Horizons DE441 state vectors (Earth centre minus Sun, both
  barycentric, ecliptic J2000) read from app/public/data/ephemeris/planets.bin (daily rows with
  velocities, 2020-01-01..2036-12-31), cubic-Hermite interpolated (error < 1e-9 AU).  This is
  the true geocentre, so no Earth-Moon barycentre offset.  Used for every JD inside the file's
  range, which covers the whole analysis window.
* Earth (fallback, only if the ephemeris files are absent / JD outside their range): JPL/SSD "Keplerian Elements for Approximate Positions of the Major Planets"
  (E.M. Standish, https://ssd.jpl.nasa.gov/planets/approx_pos.html), Table 1
  (valid 1800-2050), Earth-Moon barycentre row:
        a      = 1.00000261  + 0.00000562 T      [AU]
        e      = 0.01671123  - 0.00004392 T
        I      = -0.00001531 - 0.01294668 T      [deg]
        L      = 100.46457166 + 35999.37244981 T [deg]  (mean longitude)
        varpi  = 102.93768193 + 0.32327364 T     [deg]  (longitude of perihelion)
        Omega  = 0.0                             [deg]
  T = Julian centuries (36525 d) of TDB since J2000 (JD 2451545.0); omega = varpi - Omega,
  M = L - varpi.  This fallback is the EM barycentre with ~1e-4 AU formal accuracy (secular
  model plus the ~3e-5 AU barycentre offset); it is documented, not hidden.
* Window: 10 x 365 = 3650 days from START_JD (2026-10-09 00:00 TDB, JD 2461322.5), sampled
  every 5 days (731 samples).  Every local minimum of d(t) below 0.3 AU (plus the global
  minimum) is refined with a 0.05-day grid and then a 0.001-day grid, so the reported
  minimum distance/date is accurate to the integration model, not the 5-day sampling.

Limitations (documented, not hidden)
------------------------------------
* Pure two-body motion: no planetary perturbations, no Earth gravitational focusing, no
  non-gravitational forces.  Over 10 years this is a screening-quality prediction, not an
  ephemeris; real close-approach distances can differ.
* Elements are fetched with `full-prec=true` (full double precision), so the remaining error
  for close approaches is the two-body model itself, not catalogue rounding.
"""
from __future__ import annotations

import json
import math
import sys
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent

K_GAUSS = 0.01720209895          # rad/day, Gaussian gravitational constant
J2000 = 2451545.0
START_JD = 2461322.5             # 2026-10-09 00:00 TDB
WINDOW_DAYS = 3650.0             # 10 x 365 d
STEP_DAYS = 5.0
SHORT_DAYS = 730.0               # "short window" for the slope (first 2 years)
IMPACT_AU = 0.05                 # impact-capable threshold (MOID and min geocentric distance)
STRICT_AU = 0.001                # strict collision-course threshold
REFINE_CAND_AU = 0.3             # coarse local minima below this are refined
AU_KM = 149597870.7
JD_UNIX = 2440587.5


# ----------------------------------------------------------------------------------------
# Kepler / orbital geometry
# ----------------------------------------------------------------------------------------
def solve_kepler(M, e, tol=1e-13, max_iter=60):
    """Solve M = E - e sin E (radians) by Newton-Raphson, vectorised (arrays broadcast)."""
    M = np.asarray(M, dtype=float)
    e = np.asarray(e, dtype=float)
    M = np.mod(M + np.pi, 2 * np.pi) - np.pi                    # wrap to [-pi, pi]
    E = M + 0.85 * e * np.sign(np.sin(M) + (np.sin(M) == 0))     # Danby start
    for _ in range(max_iter):
        f = E - e * np.sin(E) - M
        dE = f / (1.0 - e * np.cos(E))
        E = E - dE
        if np.max(np.abs(dE)) < tol:
            break
    return E


def orientation_vectors(i_deg, om_deg, w_deg):
    """Unit vectors P (towards perihelion) and Q (90 deg ahead in orbital plane), ecliptic frame."""
    i, om, w = np.radians(i_deg), np.radians(om_deg), np.radians(w_deg)
    cO, sO, cw, sw, ci, si = np.cos(om), np.sin(om), np.cos(w), np.sin(w), np.cos(i), np.sin(i)
    P = np.stack([cw * cO - sw * sO * ci, cw * sO + sw * cO * ci, sw * si], axis=-1)
    Q = np.stack([-sw * cO - cw * sO * ci, -sw * sO + cw * cO * ci, cw * si], axis=-1)
    return P, Q


def kepler_position(a, e, P, Q, M_rad):
    """Heliocentric ecliptic position (..., 3) given a, e, P, Q (broadcast against M_rad)."""
    E = solve_kepler(M_rad, e)
    xp = a * (np.cos(E) - e)
    yp = a * np.sqrt(1.0 - e * e) * np.sin(E)
    return xp[..., None] * P + yp[..., None] * Q


def asteroid_positions(el, jd):
    """el: dict of arrays shape (n,) ; jd: (n,L) or (L,) -> positions (n,L,3)."""
    a = el["a"][:, None]
    e = el["e"][:, None]
    n = K_GAUSS / el["a"] ** 1.5                                  # rad/day
    jd = np.asarray(jd, dtype=float)
    if jd.ndim == 1:
        jd = jd[None, :]
    M = np.radians(el["ma"])[:, None] + n[:, None] * (jd - el["epoch"][:, None])
    return kepler_position(a, e, el["P"][:, None, :], el["Q"][:, None, :], M)


def earth_elements(jd):
    """Standish mean elements of the EM barycentre at JD (scalar/array) -> a,e,I,L,varpi,Omega."""
    T = (np.asarray(jd, dtype=float) - J2000) / 36525.0
    a = 1.00000261 + 0.00000562 * T
    e = 0.01671123 - 0.00004392 * T
    I = -0.00001531 - 0.01294668 * T
    L = 100.46457166 + 35999.37244981 * T
    varpi = 102.93768193 + 0.32327364 * T
    Om = 0.0 * T
    return a, e, I, L, varpi, Om


def earth_position_standish(jd):
    """Heliocentric ecliptic J2000 position of the EM barycentre from Standish mean elements."""
    a, e, I, L, varpi, Om = earth_elements(jd)
    w = varpi - Om
    M = np.radians(L - varpi)
    P, Q = orientation_vectors(I, Om, w)
    return kepler_position(a, e, P, Q, M)


EPH_DIR = ROOT / "app" / "public" / "data" / "ephemeris"
_EPH = {}


def _load_ephemeris():
    """Load Earth and Sun barycentric states from the DE441 file, or None if unavailable."""
    if "v" in _EPH:
        return _EPH["v"]
    v = None
    try:
        man = json.load(open(EPH_DIR / "manifest.json"))
        m = man["main"]
        bodies = m["bodies"]
        raw = np.fromfile(EPH_DIR / m["file"], dtype=m["dtype"])
        arr = raw.reshape(len(bodies), m["count"], 6)
        v = {"earth": arr[bodies.index("earth")].copy(), "sun": arr[bodies.index("sun")].copy(),
             "jd0": float(m["jd0"]), "step": float(m["stepDays"]), "count": int(m["count"])}
    except Exception:                                              # missing / malformed -> fallback
        v = None
    _EPH["v"] = v
    return v


def earth_source():
    return "JPL Horizons DE441 (Earth centre)" if _load_ephemeris() is not None else "Standish EM barycentre"


def _hermite(state, jd0, h, jd):
    """Cubic Hermite interpolation of [x,y,z,vx,vy,vz] rows; jd array (any shape) -> (...,3)."""
    n = state.shape[0]
    u = (jd - jd0) / h
    i = np.clip(np.floor(u).astype(int), 0, n - 2)
    s = (u - i)[..., None]
    p0, p1 = state[i, :3], state[i + 1, :3]
    m0, m1 = state[i, 3:] * h, state[i + 1, 3:] * h
    s2, s3 = s * s, s * s * s
    return ((2 * s3 - 3 * s2 + 1) * p0 + (s3 - 2 * s2 + s) * m0 + (-2 * s3 + 3 * s2) * p1 + (s3 - s2) * m1)


def earth_position(jd):
    """Heliocentric ecliptic J2000 position of Earth, shape jd.shape + (3,).
    DE441 (Earth centre - Sun) where the ephemeris file covers the date, else Standish."""
    jd = np.asarray(jd, dtype=float)
    eph = _load_ephemeris()
    if eph is None:
        return earth_position_standish(jd)
    jd_hi = eph["jd0"] + (eph["count"] - 1) * eph["step"]
    inside = (jd >= eph["jd0"]) & (jd <= jd_hi)
    if inside.all():
        return _hermite(eph["earth"], eph["jd0"], eph["step"], jd) - _hermite(eph["sun"], eph["jd0"], eph["step"], jd)
    out = earth_position_standish(jd)
    if inside.any():
        jj = jd[inside]
        out[inside] = _hermite(eph["earth"], eph["jd0"], eph["step"], jj) - _hermite(eph["sun"], eph["jd0"], eph["step"], jj)
    return out


# ----------------------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------------------
def prepare_elements(objs):
    el = {k: np.array([o[k] for o in objs], dtype=float) for k in ("a", "e", "i", "om", "w", "ma", "epoch")}
    el["P"], el["Q"] = orientation_vectors(el["i"], el["om"], el["w"])
    return el


def sub_elements(el, idx):
    return {k: v[idx] for k, v in el.items()}


def jd_to_iso(jd):
    dt = datetime(1970, 1, 1, tzinfo=timezone.utc) + timedelta(days=float(jd) - JD_UNIX)
    return dt.strftime("%Y-%m-%dT%H:%M:%SZ")


def estimate_diameter(H):
    """Diameter [km] from absolute magnitude, albedo 0.14: D = 1329/sqrt(p) * 10^(-H/5)."""
    return 1329.0 / math.sqrt(0.14) * 10 ** (-H / 5.0)


def geo_distance(el, jd_rows):
    """Geocentric distance for asteroids (n,) at per-row times (n,L). Returns (n,L)."""
    pos = asteroid_positions(el, jd_rows)
    ep = earth_position(jd_rows)
    return np.linalg.norm(pos - ep, axis=-1)


def refine_minimum(el, idx_obj, t_center, lo_jd, hi_jd):
    """Refine closest approach for candidate rows. Two-stage dense grid search.
    el: full element dict; idx_obj: (m,) object indices; t_center: (m,) JD of coarse minimum."""
    sel = sub_elements(el, idx_obj)
    t = t_center
    for half, npts in ((STEP_DAYS, 201), (0.05, 101)):
        offs = np.linspace(-half, half, npts)
        grid = np.clip(t[:, None] + offs[None, :], lo_jd, hi_jd)
        d = geo_distance(sel, grid)
        j = np.argmin(d, axis=1)
        t = grid[np.arange(len(t)), j]
        dmin = d[np.arange(len(t)), j]
    return dmin, t


def orbit_moid(a, e, i, om, w, jd_earth=START_JD, n0=1500, rounds=4):
    """Minimum orbit intersection distance [AU] between an asteroid's osculating ellipse and
    the Earth's mean orbit at jd_earth (Standish elements frozen at that date). Grid search
    over both eccentric anomalies followed by successive local refinement."""
    ea, ee, eI, eL, evp, eOm = [float(np.asarray(v)) for v in earth_elements(jd_earth)]
    Pa, Qa = orientation_vectors(i, om, w)
    Pe, Qe = orientation_vectors(eI, eOm, evp - eOm)

    def pts(aa, ee_, P, Q, E):
        return (aa * (np.cos(E) - ee_))[:, None] * P + (aa * math.sqrt(1 - ee_ ** 2) * np.sin(E))[:, None] * Q

    u = np.linspace(0, 2 * np.pi, n0, endpoint=False)
    A = pts(a, e, Pa, Qa, u)
    B = pts(ea, ee, Pe, Qe, u)
    D = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=-1)
    ia, ib = np.unravel_index(np.argmin(D), D.shape)
    ua, ub, h = u[ia], u[ib], 2 * np.pi / n0
    best = D[ia, ib]
    for _ in range(rounds):
        ga = ua + np.linspace(-2 * h, 2 * h, 41)
        gb = ub + np.linspace(-2 * h, 2 * h, 41)
        D = np.linalg.norm(pts(a, e, Pa, Qa, ga)[:, None, :] - pts(ea, ee, Pe, Qe, gb)[None, :, :], axis=-1)
        ja, jb = np.unravel_index(np.argmin(D), D.shape)
        ua, ub, best = ga[ja], gb[jb], D[ja, jb]
        h = 4 * h / 40
    return float(best)


# ----------------------------------------------------------------------------------------
# Main analysis
# ----------------------------------------------------------------------------------------
def analyse(objs, chunk=1500, verbose=True):
    n = len(objs)
    el = prepare_elements(objs)
    hi_jd = START_JD + WINDOW_DAYS
    tgrid = START_JD + np.arange(0.0, WINDOW_DAYS + 1e-9, STEP_DAYS)
    K = len(tgrid)
    k_short = int(round(SHORT_DAYS / STEP_DAYS)) + 1
    earth = earth_position(tgrid)                                   # (K,3)
    ts = (tgrid[:k_short] - tgrid[0]) / 365.25                        # years
    tsc = ts - ts.mean()

    out = {
        "d_start": np.empty(n), "d_end": np.empty(n), "slope": np.empty(n),
        "min_dist": np.empty(n), "min_jd": np.empty(n), "coarse_min": np.empty(n),
    }
    t0 = time.time()
    for s in range(0, n, chunk):
        e_ = min(n, s + chunk)
        idx = np.arange(s, e_)
        sel = sub_elements(el, idx)
        pos = asteroid_positions(sel, tgrid)                         # (c,K,3)
        d = np.linalg.norm(pos - earth[None], axis=-1)                # (c,K)
        out["d_start"][s:e_] = d[:, 0]
        out["d_end"][s:e_] = d[:, -1]
        y = d[:, :k_short]
        out["slope"][s:e_] = (y - y.mean(1, keepdims=True)) @ tsc / (tsc @ tsc)   # AU / yr
        out["coarse_min"][s:e_] = d.min(1)

        # candidate minima: interior local minima < REFINE_CAND_AU plus global argmin
        loc = np.zeros_like(d, dtype=bool)
        loc[:, 1:-1] = (d[:, 1:-1] <= d[:, :-2]) & (d[:, 1:-1] <= d[:, 2:]) & (d[:, 1:-1] < REFINE_CAND_AU)
        gi = np.argmin(d, axis=1)
        loc[np.arange(len(idx)), gi] = True
        ci, ck = np.nonzero(loc)
        dmin, tmin = refine_minimum(el, idx[ci], tgrid[ck], START_JD, hi_jd)
        best_d = np.full(len(idx), np.inf)
        best_t = np.zeros(len(idx))
        np.minimum.at(best_d, ci, dmin)                               # per-object minimum over candidates
        hit = dmin == best_d[ci]                                      # candidate rows attaining it
        best_t[ci[hit]] = tmin[hit]                                   # matching argmin time
        out["min_dist"][s:e_] = best_d
        out["min_jd"][s:e_] = best_t
        if verbose:
            print(f"  {e_:>6}/{n}  ({time.time() - t0:5.1f}s)", flush=True)
    out["trend"] = out["d_end"] - out["d_start"]
    return el, tgrid, out


def categorise(objs, out):
    """Return arrays: category (str), impact_capable, strict, moid. Every object gets exactly one.

    impact      = can crash: (MOID <= 0.05 AU or PHA or MOID unknown) AND window min distance < 0.05 AU
    receding    = never crashes: MOID > 0.05 AU and not PHA, AND trend_au = d(end)-d(start) >= 0
    approaching = everything else (diminishing but safe, or low-MOID objects not closing and not impact)
    trend_au is the value as shipped in asteroids.json (rounded to 4 dp) so labels match data.
    """
    moid = np.array([np.nan if o.get("moid") is None else o["moid"] for o in objs], dtype=float)
    pha = np.array([bool(o.get("pha")) for o in objs])
    min_dist = out["min_dist"]
    orbit_ok = np.where(np.isnan(moid), True, moid <= IMPACT_AU) | pha
    impact_capable = orbit_ok & (min_dist < IMPACT_AU)
    strict = impact_capable & (min_dist < STRICT_AU)
    trend = np.round(out["trend"], 4)
    safe = ~orbit_ok
    cat = np.where(impact_capable, "impact", np.where(safe & (trend >= 0), "receding", "approaching"))
    return cat, impact_capable, strict, moid


def downsample_polylines(el, idx, out, npts=60):
    """Per-object polylines: helio_orbit (one orbital period), geo_window (full 10 yr),
    geo_approach (+-60 d around closest approach). Rounded to 1e-4 AU (1e-5 near Earth)."""
    sel = sub_elements(el, np.asarray(idx))
    m = len(idx)
    n_rad = K_GAUSS / sel["a"] ** 1.5
    period = 2 * np.pi / n_rad
    # heliocentric: one full orbit starting at START_JD
    frac = np.linspace(0, 1, npts)
    jd_h = START_JD + period[:, None] * frac[None, :]
    helio = asteroid_positions(sel, jd_h)
    # geocentric across window
    jd_w = np.broadcast_to(START_JD + np.linspace(0, WINDOW_DAYS, npts), (m, npts))
    geo_w = asteroid_positions(sel, jd_w) - earth_position(jd_w)
    # geocentric around the closest approach
    jd_a = np.clip(out["min_jd"][np.asarray(idx)][:, None] + np.linspace(-60, 60, npts)[None, :],
                   START_JD, START_JD + WINDOW_DAYS)
    geo_a = asteroid_positions(sel, jd_a) - earth_position(jd_a)
    return helio, geo_w, geo_a


def main():
    objs = json.load(open(ROOT / "data" / "neo_elements.json"))
    print(f"{len(objs):,} NEOs; window {jd_to_iso(START_JD)} + {WINDOW_DAYS:.0f} d, step {STEP_DAYS} d")
    el, tgrid, out = analyse(objs)
    cat, impact_capable, strict, moid = categorise(objs, out)
    counts = {c: int((cat == c).sum()) for c in ("impact", "approaching", "receding")}
    assert sum(counts.values()) == len(objs)
    n_strict = int(strict.sum())
    print("categories:", counts, "strict:", n_strict)

    # ---- asteroids.json (dashboard contract) ----
    rows = []
    for k, o in enumerate(objs):
        dia = o.get("diameter_km")
        if dia is None:
            dia = estimate_diameter(o["H"]) if o.get("H") is not None else 0.05
        r = lambda v, p: round(float(v), p)
        rows.append({
            "id": o["id"], "name": o["name"], "diameter_km": r(dia, 5), "pha": bool(o["pha"]),
            "category": str(cat[k]), "min_dist_au": r(out["min_dist"][k], 6),
            "min_date": jd_to_iso(out["min_jd"][k]),
            "a": r(o["a"], 6), "e": r(o["e"], 6), "i": r(o["i"], 4), "om": r(o["om"], 4),
            "w": r(o["w"], 4), "ma": r(o["ma"], 4), "epoch": o["epoch"],
            # additive extras (not in the minimal contract)
            "strict_collision": bool(strict[k]),
            "trend_au": r(out["trend"][k], 4),
        })
    pub = ROOT / "app" / "public" / "data"
    pub.mkdir(parents=True, exist_ok=True)
    with open(pub / "asteroids.json", "w") as f:
        json.dump(rows, f, separators=(",", ":"))
    print("wrote", pub / "asteroids.json")

    # ---- trajectories for 300 closest/impact objects ----
    impact_idx = np.nonzero(cat == "impact")[0]
    impact_idx = impact_idx[np.argsort(out["min_dist"][impact_idx])]
    rest = np.nonzero(cat != "impact")[0]
    rest = rest[np.argsort(out["min_dist"][rest])]
    sel_idx = np.concatenate([impact_idx, rest])[:300]
    helio, geo_w, geo_a = downsample_polylines(el, sel_idx, out)
    ej = START_JD + np.linspace(0, 365.25, 60)
    earth_orbit = earth_position(ej)
    flat = lambda arr, p: [round(float(v), p) for v in arr.reshape(-1)]
    traj = {
        "meta": {
            "start_jd": START_JD, "start_iso": jd_to_iso(START_JD), "window_days": WINDOW_DAYS,
            "units": "AU, ecliptic J2000, flat [x,y,z,...] x 60 points",
            "helio_orbit": "heliocentric, 60 points over one full orbital period starting at start_jd",
            "geo_window": "geocentric (asteroid minus Earth), 60 points evenly over the 10-year window",
            "geo_approach": "geocentric, 60 points over +-60 days around closest approach (clipped to window)",
            "earth_orbit": "heliocentric Earth, 60 points over one year from start_jd",
            "selection": "all impact-category objects by ascending min_dist, then closest others; 300 total",
        },
        "earth_orbit": flat(earth_orbit, 5),
        "objects": [
            {"id": objs[k]["id"], "name": objs[k]["name"], "category": str(cat[k]),
             "min_dist_au": round(float(out["min_dist"][k]), 6), "min_date": jd_to_iso(out["min_jd"][k]),
             "helio_orbit": flat(helio[j], 4), "geo_window": flat(geo_w[j], 5), "geo_approach": flat(geo_a[j], 6)}
            for j, k in enumerate(sel_idx)
        ],
    }
    with open(ROOT / "data" / "asteroids_trajectories.json", "w") as f:
        json.dump(traj, f, separators=(",", ":"))
    print("wrote data/asteroids_trajectories.json", len(traj["objects"]), "objects")

    write_report(objs, out, cat, impact_capable, strict, moid, counts, n_strict)


def write_report(objs, out, cat, impact_capable, strict, moid, counts, n_strict):
    n = len(objs)
    order = np.argsort(out["min_dist"])[:20]
    rep = f"""# Trajectory report

Generated by `python -I scripts/trajectory.py` from `data/neo_elements.json` ({n:,} NEOs, JPL SBDB).

## Method
- Two-body Keplerian propagation of each asteroid from its catalogue elements, fetched with `full-prec=true` (full double precision; Newton-Raphson Kepler solver, mean motion from `a` via the Gaussian constant). Earth: **{earth_source()}** - JPL Horizons DE441 Earth-centre minus Sun state vectors from `app/public/data/ephemeris/planets.bin`, cubic-Hermite interpolated (the true geocentre, no barycentre offset). If those files are absent the code falls back to Standish (JPL SSD "Approximate Positions of the Major Planets", Table 1, 1800-2050) mean elements of the Earth-Moon barycentre, which carries a ~1e-4 AU error (secular model + ~3e-5 AU barycentre offset).
- Window: {jd_to_iso(START_JD)} (JD {START_JD}) + {WINDOW_DAYS:.0f} days (10 x 365 d), 5-day steps (731 samples); every coarse local minimum of the geocentric distance below {REFINE_CAND_AU} AU is refined with 0.05-day and 0.001-day grids, so the minimum distance/date are not limited by the 5-day sampling. Per-object minimum over candidates uses `np.minimum.at` with the matching argmin time.
- Per object: `d_start`, `d_end`, `trend_au = d(end) - d(start)` (as shipped in `asteroids.json`, rounded to 4 dp), `min_dist` and `min_date`.

## Category rules (exactly one per object)
1. **impact** (can crash) = (MOID <= {IMPACT_AU} AU, or flagged PHA, or MOID unknown) AND predicted minimum geocentric distance in the window < {IMPACT_AU} AU. Independent of the trend.
2. **receding** (never crashes) = NOT impact-capable by orbit (MOID > {IMPACT_AU} AU and not PHA) AND `trend_au >= 0`. A PHA or MOID <= {IMPACT_AU} object is never `receding`.
3. **approaching** = everything else (distance diminishing but safe, or low-MOID/PHA objects that do not come within {IMPACT_AU} AU in the window, regardless of trend).
4. `strict_collision` flag (subset of `impact`) = min_dist < {STRICT_AU} AU. Only meaningful with full-precision elements; it is still a two-body, unfocused, unperturbed prediction, not an impact probability.

## Counts
| category | count | share |
|---|---:|---:|
| impact | {counts['impact']:,} | {100*counts['impact']/n:.1f}% |
| approaching | {counts['approaching']:,} | {100*counts['approaching']/n:.1f}% |
| receding | {counts['receding']:,} | {100*counts['receding']/n:.1f}% |
| **total** | **{n:,}** | |

- strict_collision subset of `impact` (min_dist < {STRICT_AU} AU): **{n_strict:,}**
- impact-capable (= impact count): {int(impact_capable.sum()):,}
- PHA-flagged objects: {int(sum(1 for o in objs if o['pha'])):,}; PHAs in `receding`: **{int(sum(1 for k, o in enumerate(objs) if o['pha'] and cat[k] == 'receding'))}**; objects with MOID <= {IMPACT_AU} AU in `receding`: **{int(np.sum((cat == 'receding') & (np.nan_to_num(moid, nan=0.0) <= IMPACT_AU)))}**

## Top 20 closest approaches in the window
| # | name | category | min dist (AU) | min dist (km) | date | MOID (AU) | PHA |
|--:|---|---|---:|---:|---|---:|:-:|
"""
    for r, k in enumerate(order, 1):
        o = objs[k]
        mo = "n/a" if o.get("moid") is None else f"{o['moid']:g}"
        rep += (f"| {r} | {o['name'].strip()} | {cat[k]} | {out['min_dist'][k]:.5f} | "
                f"{out['min_dist'][k]*AU_KM:,.0f} | {jd_to_iso(out['min_jd'][k])[:10]} | {mo} | {'Y' if o['pha'] else ''} |\n")
    rep += f"""
## Accuracy caveats
- Pure two-body model: no planetary perturbations, no Earth gravitational focusing, no Yarkovsky. Treat the predicted miss distances as screening-grade; for deep encounters (e.g. Apophis 2029) the unfocused two-body distance is ~25% larger than the real focused one.
- Elements are full precision (`full-prec=true`), so catalogue rounding is no longer a limiting error.
- Earth position source: {earth_source()}. Fallback Standish accuracy is ~1e-4 AU.
"""
    (ROOT / "data" / "trajectory_report.md").write_text(rep, encoding="utf-8")
    print("wrote data/trajectory_report.md")


if __name__ == "__main__":
    sys.exit(main())
