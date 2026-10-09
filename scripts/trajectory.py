"""Trajectory computation + categorisation of all NEOs over 2026-10-09 .. 2100-12-31.

Run:  python -I scripts/trajectory.py          (reads data/neo_elements.json)

Method
------
* Two-body (Keplerian) propagation of every asteroid from its catalogue osculating
  elements (ecliptic J2000; a [AU], e, i/Omega/omega/M [deg], epoch JD TDB).  Mean motion
  is n = k / a^1.5 (Gaussian gravitational constant k = 0.01720209895 rad/day); the
  catalogue `per` field is NOT used (SBDB rounds it to 3 significant digits).
* Kepler's equation M = E - e sin E is solved with vectorised Halley iteration (Danby
  starting guess, tolerance 1e-13 rad) for stand-alone evaluations.  The coarse scan marches
  the 5-day grid in time and predicts E(t+dt) from E(t) with a 2nd-order Taylor step, so
  only two Halley iterations per grid point are needed (stragglers iterate to convergence);
  it evaluates the geocentric distance through |r|^2 + |R|^2 - 2 r.R without building
  (n, K, 3) arrays, and works on cache-sized blocks of objects.
* Earth, best available source, per Julian date:
    1. JPL Horizons DE441 state vectors (Earth centre minus Sun, both barycentric, ecliptic
       J2000) from app/public/data/ephemeris/planets.bin (daily rows with velocities),
       cubic-Hermite interpolated (error < 1e-9 AU).  True geocentre, no EM-barycentre offset.
       Used for every JD inside the file's range (the manifest's jd0 .. jd0 + (count-1) step).
    2. Beyond the file's last row (the file originally ended 2036-12-31 and was then extended
       to 2100-12-31 = END_JD, so this branch is only a fallback now; it is also exercised by
       the tests) the JPL/SSD "Keplerian
       Elements for Approximate Positions of the Major Planets" (E.M. Standish,
       https://ssd.jpl.nasa.gov/planets/approx_pos.html), Table 1 (valid 1800-2050),
       Earth-Moon barycentre row:
            a      = 1.00000261  + 0.00000562 T      [AU]
            e      = 0.01671123  - 0.00004392 T
            I      = -0.00001531 - 0.01294668 T      [deg]
            L      = 100.46457166 + 35999.37244981 T [deg]  (mean longitude)
            varpi  = 102.93768193 + 0.32327364 T     [deg]  (longitude of perihelion)
            Omega  = 0.0                             [deg]
       T = Julian centuries (36525 d) of TDB since J2000 (JD 2451545.0); omega = varpi - Omega,
       M = L - varpi.  EM barycentre, ~1e-4 AU formal accuracy within 1800-2050 and slowly
       growing when extrapolated to 2100 (still far below the two-body model error).
  `earth_source()` reports what was actually used; the report records it.  Standish agrees with
  DE441 to <= 1.5e-4 AU over 2020-2100 (rms 5e-5 AU), so which branch is used does not matter
  at the model's accuracy; the loader picks up any ephemeris range change with no code change.
* Window: START_JD (2026-10-09 00:00 TDB, JD 2461322.5) .. END_JD (2100-12-31 00:00 TDB,
  JD 2488433.5), 27111 days, sampled every 5 days (+ the end point; 5424 samples).  Every
  coarse local minimum of d(t) below 0.3 AU (plus the global minimum, plus the window-edge
  samples if they are minima) is refined with three 21-point grids (0.5 d, 0.06 d, 0.006 d
  spacing) and a parabola through d^2, so the reported minimum distance/date are accurate to
  the integration model, not to the 5-day sampling.  Every refined minimum with d < 0.05 AU
  strictly inside the window is a "close approach"; `next_close_approach` ships the first 3.
* Trend: `trend_au_10y` = d(start+3650 d) - d(start) is the signal used for
  approaching/receding; `trend_au_full` = d(2100-12-31) - d(start) is shipped for reference
  only (over 74 years it is dominated by orbital phase, not by "closing" or "receding").

Limitations (documented, not hidden)
------------------------------------
* Pure two-body motion: no planetary perturbations, no Earth gravitational focusing, no
  non-gravitational forces.  Over 10 years this is a screening-quality prediction; over 74
  years the osculating semi-major axis error alone (d a/a ~ 1e-3 -> d(mean motion)/n ~ 1.5e-3)
  moves a 1-year-period body by ~0.7 rad of orbital phase, so flyby *dates* beyond ~2040 are
  indicative of orbit geometry, not predictions, and objects that actually pass close to a
  planet change orbit completely (e.g. Apophis in 2029).
* Elements are fetched with `full-prec=true` (full double precision), so for the first
  decade the remaining error for close approaches is the two-body model itself.
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
END_JD = 2488433.5               # 2100-12-31 00:00 TDB
WINDOW_DAYS = END_JD - START_JD  # 27111 d (~74.2 yr)
TREND_DAYS = 3650.0              # 10 x 365 d: window of the approaching/receding trend
STEP_DAYS = 5.0
IMPACT_AU = 0.05                 # impact-capable threshold (MOID and min geocentric distance)
CA_AU = 0.05                     # a "close approach" in next_close_approach: refined minimum < this
N_CA = 3                         # next_close_approach list length
STRICT_AU = 0.001                # strict collision-course threshold
REFINE_CAND_AU = 0.3             # coarse local minima below this are refined
REFINE_HALF = (5.0, 0.6, 0.06)   # half-widths [d] of the three 21-point refinement grids
DEDUPE_DAYS = 2.0                # two refined minima of one object closer than this are the same event
AU_KM = 149597870.7
JD_UNIX = 2440587.5


# ----------------------------------------------------------------------------------------
# Kepler / orbital geometry
# ----------------------------------------------------------------------------------------
def solve_kepler(M, e, tol=1e-13, max_iter=60):
    """Solve M = E - e sin E (radians) by Halley iteration, vectorised (arrays broadcast).
    M is wrapped to [-pi, pi]; the returned E is in the same wrapped domain."""
    M = np.asarray(M, dtype=float)
    e = np.asarray(e, dtype=float)
    M = np.mod(M + np.pi, 2 * np.pi) - np.pi                    # wrap to [-pi, pi]
    E = M + 0.85 * e * np.sign(np.sin(M) + (np.sin(M) == 0))     # Danby start
    for _ in range(max_iter):
        s = np.sin(E)
        f = E - e * s - M
        f1 = 1.0 - e * np.cos(E)
        dE = f / (f1 - 0.5 * f * (e * s) / f1)
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
    """Load Earth-minus-Sun states from the DE441 file, or None if unavailable.

    The file may be rewritten by another process (range extension), so the manifest and
    binary are only accepted as a consistent pair (file size == bodies*count*6*itemsize);
    otherwise we retry for a while before falling back to Standish (with a loud warning)."""
    if "v" in _EPH:
        return _EPH["v"]
    v = None
    err = None
    for attempt in range(15):
        try:
            man = json.load(open(EPH_DIR / "manifest.json"))
            m = man["main"]
            bodies = m["bodies"]
            path = EPH_DIR / m["file"]
            expect = len(bodies) * int(m["count"]) * 6 * np.dtype(m["dtype"]).itemsize
            if path.stat().st_size != expect:
                err = f"{path.name} is {path.stat().st_size} bytes, manifest implies {expect}"
                time.sleep(2.0)
                continue
            raw = np.fromfile(path, dtype=m["dtype"])
            arr = raw.reshape(len(bodies), int(m["count"]), 6)
            rel = arr[bodies.index("earth")] - arr[bodies.index("sun")]      # Earth centre - Sun, (N,6)
            v = {"rel": np.ascontiguousarray(rel), "jd0": float(m["jd0"]), "step": float(m["stepDays"]),
                 "count": int(m["count"])}
            err = None
            break
        except Exception as ex:                                   # missing / malformed
            err = repr(ex)
            break
    if v is None and err:
        print(f"WARNING: DE441 ephemeris unusable ({err}); falling back to Standish for Earth", file=sys.stderr, flush=True)
    _EPH["v"] = v
    return v


def ephemeris_range():
    """(jd_lo, jd_hi) covered by the DE441 file, or None."""
    eph = _load_ephemeris()
    if eph is None:
        return None
    return eph["jd0"], eph["jd0"] + (eph["count"] - 1) * eph["step"]


def earth_source():
    r = ephemeris_range()
    if r is None:
        return "Standish EM barycentre (no ephemeris file) for the whole window"
    lo, hi = r
    s = f"JPL Horizons DE441 (Earth centre) {jd_to_iso(lo)[:10]}..{jd_to_iso(hi)[:10]}"
    if hi >= END_JD and lo <= START_JD:
        return s + " (covers the whole window)"
    parts = []
    if lo > START_JD:
        parts.append(f"{jd_to_iso(START_JD)[:10]}..{jd_to_iso(lo)[:10]}")
    if hi < END_JD:
        parts.append(f"{jd_to_iso(hi)[:10]}..{jd_to_iso(END_JD)[:10]}")
    return s + "; Standish EM barycentre for " + " and ".join(parts)


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
        return _hermite(eph["rel"], eph["jd0"], eph["step"], jd)
    out = earth_position_standish(jd)
    if inside.any():
        out[inside] = _hermite(eph["rel"], eph["jd0"], eph["step"], jd[inside])
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


def iso_to_jd(iso):
    """ISO 'YYYY-MM-DD[THH:MM:SSZ]' -> Julian date (UTC treated as TDB)."""
    dt = datetime.fromisoformat(iso.replace("Z", "+00:00"))
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return JD_UNIX + (dt - datetime(1970, 1, 1, tzinfo=timezone.utc)).total_seconds() / 86400.0


def estimate_diameter(H):
    """Diameter [km] from absolute magnitude, albedo 0.14: D = 1329/sqrt(p) * 10^(-H/5)."""
    return 1329.0 / math.sqrt(0.14) * 10 ** (-H / 5.0)


def make_grid():
    """Coarse time grid: START_JD + 0, 5, 10 ... days, plus the window end point."""
    t = START_JD + np.arange(0.0, WINDOW_DAYS, STEP_DAYS)
    return np.append(t, END_JD)


def geo_distance(el, jd_rows):
    """Geocentric distance for asteroids (n,) at per-row times (n,L). Returns (n,L)."""
    pos = asteroid_positions(el, jd_rows)
    ep = earth_position(jd_rows)
    return np.linalg.norm(pos - ep, axis=-1)


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
# Coarse scan (streaming over the 5-day grid, vectorised over a block of objects)
# ----------------------------------------------------------------------------------------
def scan_block(sel, tgrid, earth, k10):
    """March the coarse grid for one block of objects.

    Returns d_start, d_10y (at index k10), d_end, global-minimum value/index and the coarse
    local-minimum candidates (object index within block, grid index) below REFINE_CAND_AU,
    including the two window-edge samples when they are minima."""
    a, e = sel["a"], sel["e"]
    m, K = a.size, len(tgrid)
    n = K_GAUSS / a ** 1.5
    ae = a * e
    b = a * np.sqrt(1.0 - e * e)
    P = np.ascontiguousarray(sel["P"])
    Q = np.ascontiguousarray(sel["Q"])
    t0 = tgrid[0]
    M0 = np.radians(sel["ma"]) + n * (t0 - sel["epoch"])           # unwrapped mean anomaly at t0
    Mw = np.mod(M0 + np.pi, 2 * np.pi) - np.pi
    E = solve_kepler(Mw, e) + (M0 - Mw)                             # unwrapped eccentric anomaly
    sE, cE = np.sin(E), np.cos(E)
    R2 = np.einsum("ij,ij->i", earth, earth)

    d_start = np.empty(m)
    d_10 = np.empty(m)
    gmin = np.full(m, np.inf)
    gk = np.zeros(m, dtype=np.int64)
    cand_o, cand_k = [], []
    d1 = d2 = None                                                   # d at k-1, k-2
    for k in range(K):
        R = earth[k]
        r = a - ae * cE
        d2sq = r * r + R2[k] - 2.0 * ((a * cE - ae) * (P @ R) + b * sE * (Q @ R))
        np.maximum(d2sq, 0.0, out=d2sq)
        d = np.sqrt(d2sq)
        if k == 0:
            d_start[:] = d
        if k == k10:
            d_10[:] = d
        better = d < gmin
        np.putmask(gmin, better, d)
        np.putmask(gk, better, k)
        if k >= 1:
            if k == 1:                                               # window-start edge sample
                lm = (d1 <= d) & (d1 < REFINE_CAND_AU)
                idx = np.flatnonzero(lm)
                if idx.size:
                    cand_o.append(idx); cand_k.append(np.zeros(idx.size, dtype=np.int64))
            else:                                                    # interior sample k-1
                lm = (d1 <= d2) & (d1 <= d) & (d1 < REFINE_CAND_AU)
                idx = np.flatnonzero(lm)
                if idx.size:
                    cand_o.append(idx); cand_k.append(np.full(idx.size, k - 1, dtype=np.int64))
        d2, d1 = d1, d
        if k == K - 1:
            break
        # ---- advance E from grid point k to k+1 (Taylor predictor + 2 Halley iterations)
        M1 = M0 + n * (tgrid[k + 1] - t0)
        f1 = 1.0 - e * cE
        g = (n * (tgrid[k + 1] - tgrid[k])) / f1
        E = E + g - 0.5 * e * sE * g * g / f1
        for _ in range(2):
            s = np.sin(E)
            c = np.cos(E)
            f = E - e * s - M1
            f1 = 1.0 - e * c
            step = f / (f1 - 0.5 * f * (e * s) / f1)
            E = E - step
        dl = -step                                                   # sin/cos at the updated E from (s, c)
        sE = s + c * dl - 0.5 * s * dl * dl
        cE = c - s * dl - 0.5 * c * dl * dl
        bad = np.abs(step) > 1e-6
        if bad.any():                                                # rare: slow convergence -> iterate exactly
            ib = np.flatnonzero(bad)
            Eb, Mb, eb = E[ib], M1[ib], e[ib]
            for _ in range(60):
                s_ = np.sin(Eb)
                f = Eb - eb * s_ - Mb
                f1_ = 1.0 - eb * np.cos(Eb)
                st = f / (f1_ - 0.5 * f * (eb * s_) / f1_)
                Eb = Eb - st
                if np.max(np.abs(st)) < 1e-11:
                    break
            E[ib] = Eb
            sE[ib], cE[ib] = np.sin(Eb), np.cos(Eb)
    # window-end edge sample
    lm = (d1 <= d2) & (d1 < REFINE_CAND_AU)
    idx = np.flatnonzero(lm)
    if idx.size:
        cand_o.append(idx); cand_k.append(np.full(idx.size, K - 1, dtype=np.int64))
    cand_o.append(np.arange(m)); cand_k.append(gk)                  # the global minimum always gets refined
    return {"d_start": d_start, "d_10": d_10, "d_end": d1, "coarse_min": gmin,
            "cand_o": np.concatenate(cand_o), "cand_k": np.concatenate(cand_k)}


# ----------------------------------------------------------------------------------------
# Refinement of candidate minima
# ----------------------------------------------------------------------------------------
def refine_candidates(el, obj, t_center, lo_jd=START_JD, hi_jd=END_JD, block=3000):
    """Refine closest approaches for candidate rows: three 21-point grids of shrinking width
    around the coarse minimum, then the vertex of a parabola through d^2 on the last grid.
    el: full element dict; obj: (m,) object indices; t_center: (m,) JD of the coarse minimum.
    Returns refined JD, distance, and `inside` (False if the minimum sits on the window edge)."""
    m = len(obj)
    t_out = np.empty(m)
    d_out = np.empty(m)
    ok_out = np.empty(m, dtype=bool)
    for s in range(0, m, block):
        sl = slice(s, min(m, s + block))
        sel = sub_elements(el, obj[sl])
        t = t_center[sl].copy()
        c = t.size
        rows = np.arange(c)
        for half in REFINE_HALF:
            offs = np.linspace(-half, half, 21)
            grid = np.clip(t[:, None] + offs[None, :], lo_jd, hi_jd)
            d = geo_distance(sel, grid)
            j = np.argmin(d, axis=1)
            t = grid[rows, j]
        # parabola through d^2 around the last grid minimum (interior, unclipped, uniform spacing)
        h = REFINE_HALF[-1] / 10.0
        y = d * d
        jm = np.clip(j, 1, 19)
        y0, y1, y2 = y[rows, jm - 1], y[rows, jm], y[rows, jm + 1]
        den = y0 - 2 * y1 + y2
        with np.errstate(divide="ignore", invalid="ignore"):
            dlt = 0.5 * h * (y0 - y2) / den
        good = (j > 0) & (j < 20) & (den > 0) & (np.abs(dlt) <= h)
        t_v = np.clip(np.where(good, t + dlt, t), lo_jd, hi_jd)
        d_v = geo_distance(sel, t_v[:, None])[:, 0]
        d_g = d[rows, j]
        use_v = good & (d_v <= d_g)
        t_out[sl] = np.where(use_v, t_v, t)
        d_out[sl] = np.where(use_v, d_v, d_g)
        ok_out[sl] = (t_out[sl] > lo_jd + 1e-3) & (t_out[sl] < hi_jd - 1e-3)
    return t_out, d_out, ok_out


# ----------------------------------------------------------------------------------------
# Main analysis
# ----------------------------------------------------------------------------------------
def analyse(objs, chunk=5000, verbose=True):
    """Coarse scan + refinement for every object.

    Returns (el, tgrid, out); `out` holds per-object d_start / d_10 / d_end, trend_10y,
    trend_full, coarse_min, min_dist / min_jd (refined global minimum over the window),
    n_ca (count of close approaches < CA_AU) and the flat, object-then-time sorted list
    of all of them: fly_obj, fly_jd, fly_d."""
    n = len(objs)
    el = prepare_elements(objs)
    tgrid = make_grid()
    K = len(tgrid)
    k10 = int(round(TREND_DAYS / STEP_DAYS))
    assert abs(tgrid[k10] - (START_JD + TREND_DAYS)) < 1e-9
    earth = earth_position(tgrid)                                   # (K,3)

    d_start, d_10, d_end, coarse_min = (np.empty(n) for _ in range(4))
    co, ck = [], []
    t0 = time.time()
    for s in range(0, n, chunk):
        e_ = min(n, s + chunk)
        r = scan_block(sub_elements(el, slice(s, e_)), tgrid, earth, k10)
        d_start[s:e_], d_10[s:e_], d_end[s:e_], coarse_min[s:e_] = r["d_start"], r["d_10"], r["d_end"], r["coarse_min"]
        co.append(r["cand_o"] + s)
        ck.append(r["cand_k"])
        if verbose:
            print(f"  scan   {e_:>6}/{n}  ({time.time() - t0:5.1f}s)", flush=True)
    key = np.unique(np.concatenate(co) * K + np.concatenate(ck))   # unique (object, grid index) candidates
    c_obj, c_k = key // K, key % K
    if verbose:
        print(f"  {len(key):,} candidate minima ({len(key) / n:.1f} / object); refining ...", flush=True)
    t_ref, d_ref, inside = refine_candidates(el, c_obj, tgrid[c_k])
    if verbose:
        print(f"  refine done ({time.time() - t0:5.1f}s)", flush=True)

    # per-object global minimum over candidates (never worse than the coarse grid value)
    best_d = np.full(n, np.inf)
    np.minimum.at(best_d, c_obj, d_ref)
    hit = d_ref == best_d[c_obj]
    best_t = np.zeros(n)
    best_t[c_obj[hit]] = t_ref[hit]
    out = {"d_start": d_start, "d_10": d_10, "d_end": d_end, "coarse_min": coarse_min,
           "min_dist": np.minimum(best_d, coarse_min), "min_jd": best_t}
    out["trend_10y"] = d_10 - d_start
    out["trend_full"] = d_end - d_start

    # all close approaches (< CA_AU, refined minimum strictly inside the window), de-duplicated
    sel_ca = np.flatnonzero((d_ref < CA_AU) & inside)
    fo, fj, fd = c_obj[sel_ca], t_ref[sel_ca], d_ref[sel_ca]
    order = np.lexsort((fj, fo))
    fo, fj, fd = fo[order], fj[order], fd[order]
    if len(fo) > 1:
        dup = (fo[1:] == fo[:-1]) & (fj[1:] - fj[:-1] < DEDUPE_DAYS)       # same event found twice
        drop = np.zeros(len(fo), dtype=bool)
        for i in np.flatnonzero(dup) + 1:
            if fd[i] < fd[i - 1]:
                drop[i - 1] = True
            else:
                drop[i] = True
        fo, fj, fd = fo[~drop], fj[~drop], fd[~drop]
    out["fly_obj"], out["fly_jd"], out["fly_d"] = fo, fj, fd
    out["n_ca"] = np.bincount(fo, minlength=n)
    return el, tgrid, out


def categorise(objs, out):
    """Return arrays: category (str), impact_capable, strict, moid. Every object gets exactly one.

    impact      = can crash: (MOID <= 0.05 AU or PHA or MOID unknown) AND window min distance < 0.05 AU
    receding    = never crashes: MOID > 0.05 AU and not PHA, AND trend_au_10y >= 0
    approaching = everything else (diminishing but safe, or low-MOID objects not closing and not impact)
    The window for the minimum distance is 2026-10-09 .. 2100-12-31; the trend is the geocentric-distance
    change over the FIRST 10 YEARS only (d(+3650 d) - d(start), rounded to 4 dp as shipped in
    asteroids.json so labels match the data).  Over 74 years d(end) - d(start) is mostly orbital phase.
    """
    moid = np.array([np.nan if o.get("moid") is None else o["moid"] for o in objs], dtype=float)
    pha = np.array([bool(o.get("pha")) for o in objs])
    min_dist = out["min_dist"]
    orbit_ok = np.where(np.isnan(moid), True, moid <= IMPACT_AU) | pha
    impact_capable = orbit_ok & (min_dist < IMPACT_AU)
    strict = impact_capable & (min_dist < STRICT_AU)
    trend = np.round(out["trend_10y"], 4)
    safe = ~orbit_ok
    cat = np.where(impact_capable, "impact", np.where(safe & (trend >= 0), "receding", "approaching"))
    return cat, impact_capable, strict, moid


def downsample_polylines(el, idx, out, npts=60):
    """Per-object polylines: helio_orbit (one orbital period), geo_window (full 74 yr),
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
                   START_JD, END_JD)
    geo_a = asteroid_positions(sel, jd_a) - earth_position(jd_a)
    return helio, geo_w, geo_a


def ca_lists(out, n):
    """CSR pointers into the flat close-approach arrays: object k owns [ptr[k], ptr[k+1])."""
    ptr = np.zeros(n + 1, dtype=np.int64)
    np.cumsum(out["n_ca"], out=ptr[1:])
    return ptr


def main():
    objs = json.load(open(ROOT / "data" / "neo_elements.json"))
    n = len(objs)
    print(f"{n:,} NEOs; window {jd_to_iso(START_JD)} .. {jd_to_iso(END_JD)} ({WINDOW_DAYS:.0f} d), step {STEP_DAYS} d", flush=True)
    print("Earth:", earth_source(), flush=True)
    t_all = time.time()
    el, tgrid, out = analyse(objs)
    cat, impact_capable, strict, moid = categorise(objs, out)
    counts = {c: int((cat == c).sum()) for c in ("impact", "approaching", "receding")}
    assert sum(counts.values()) == n
    n_strict = int(strict.sum())
    print("categories:", counts, "strict:", n_strict, f"({time.time() - t_all:.0f}s)", flush=True)
    ptr = ca_lists(out, n)

    # ---- asteroids.json (dashboard contract) ----
    rows = []
    r = lambda v, p: round(float(v), p)
    for k, o in enumerate(objs):
        dia = o.get("diameter_km")
        if dia is None:
            dia = estimate_diameter(o["H"]) if o.get("H") is not None else 0.05
        row = {
            "id": o["id"], "name": o["name"], "diameter_km": r(dia, 5), "pha": bool(o["pha"]),
            "category": str(cat[k]), "min_dist_au": r(out["min_dist"][k], 6),
            "min_date": jd_to_iso(out["min_jd"][k]),
            "a": r(o["a"], 6), "e": r(o["e"], 6), "i": r(o["i"], 4), "om": r(o["om"], 4),
            "w": r(o["w"], 4), "ma": r(o["ma"], 4), "epoch": o["epoch"],
            # additive extras (not in the minimal contract)
            "strict_collision": bool(strict[k]),
            "trend_au": r(out["trend_10y"][k], 4),            # alias of trend_au_10y (backward compatible)
            "trend_au_10y": r(out["trend_10y"][k], 4),
            "trend_au_full": r(out["trend_full"][k], 4),
        }
        c0, c1 = ptr[k], ptr[k + 1]
        if c1 > c0:                                            # first N_CA close approaches, soonest first
            row["next_close_approach"] = [
                {"date": jd_to_iso(out["fly_jd"][q]), "dist_au": r(out["fly_d"][q], 6)}
                for q in range(c0, min(c1, c0 + N_CA))]
            if c1 - c0 > N_CA:                                  # list truncated -> ship the full count
                row["n_close_approaches"] = int(c1 - c0)
        rows.append(row)
    pub = ROOT / "app" / "public" / "data"
    pub.mkdir(parents=True, exist_ok=True)
    with open(pub / "asteroids.json", "w") as f:
        json.dump(rows, f, separators=(",", ":"))
    size_mb = (pub / "asteroids.json").stat().st_size / 1e6
    print(f"wrote {pub / 'asteroids.json'} ({size_mb:.2f} MB)", flush=True)

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
            "start_jd": START_JD, "start_iso": jd_to_iso(START_JD), "end_jd": END_JD, "end_iso": jd_to_iso(END_JD),
            "window_days": WINDOW_DAYS,
            "units": "AU, ecliptic J2000, flat [x,y,z,...] x 60 points",
            "helio_orbit": "heliocentric, 60 points over one full orbital period starting at start_jd",
            "geo_window": "geocentric (asteroid minus Earth), 60 points evenly over the 74-year window (coarse: ~1.2 yr apart)",
            "geo_approach": "geocentric, 60 points over +-60 days around closest approach (clipped to window)",
            "earth_orbit": "heliocentric Earth, 60 points over one year from start_jd",
            "selection": "all impact-category objects by ascending min_dist, then closest others; 300 total",
            "earth_source": earth_source(),
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
    print("wrote data/asteroids_trajectories.json", len(traj["objects"]), "objects",
          f"({(ROOT / 'data' / 'asteroids_trajectories.json').stat().st_size / 1e6:.2f} MB)", flush=True)

    write_report(objs, out, cat, impact_capable, strict, moid, counts, n_strict, size_mb, time.time() - t_all)
    print(f"total {time.time() - t_all:.0f}s", flush=True)


def _period_table(out):
    """Close-approach statistics by calendar span."""
    spans = [("2026-10-09 .. 2036-10-07 (first 10 yr)", START_JD, START_JD + TREND_DAYS),
             ("2036-10-07 .. 2050-01-01", START_JD + TREND_DAYS, 2469807.5),
             ("2050-01-01 .. 2075-01-01", 2469807.5, 2478938.5),
             ("2075-01-01 .. 2100-12-31", 2478938.5, END_JD + 1)]
    fj, fd, fo = out["fly_jd"], out["fly_d"], out["fly_obj"]
    lines = ["| span | flybys < 0.05 AU | < 0.01 AU | < 0.002 AU | < 0.001 AU | distinct objects < 0.05 AU |",
             "|---|---:|---:|---:|---:|---:|"]
    for name, a, b in spans:
        m = (fj >= a) & (fj < b)
        lines.append(f"| {name} | {int(m.sum()):,} | {int((m & (fd < 0.01)).sum()):,} | {int((m & (fd < 0.002)).sum()):,} | "
                     f"{int((m & (fd < 0.001)).sum()):,} | {len(np.unique(fo[m])):,} |")
    return "\n".join(lines)


def write_report(objs, out, cat, impact_capable, strict, moid, counts, n_strict, size_mb=None, runtime_s=None):
    n = len(objs)
    pha = np.array([bool(o["pha"]) for o in objs])
    order = np.argsort(out["min_dist"])[:20]
    fo, fj, fd = out["fly_obj"], out["fly_jd"], out["fly_d"]
    t10_end = START_JD + TREND_DAYS

    # ---- rule re-check numbers
    orbit_ok = np.where(np.isnan(moid), True, moid <= IMPACT_AU) | pha
    n_ok = int(orbit_ok.sum())
    n_imp = int(impact_capable.sum())
    outside_hit = int(((~orbit_ok) & (out["min_dist"] < IMPACT_AU)).sum())
    # old 10-year window: minimum over [d_start, d(+10y), close approaches within 10 y]
    min10 = np.minimum(out["d_start"], out["d_10"])
    m10 = fj <= t10_end
    np.minimum.at(min10, fo[m10], fd[m10])
    imp10 = int((orbit_ok & (min10 < IMPACT_AU)).sum())
    strict10 = int((orbit_ok & (min10 < STRICT_AU)).sum())
    t10 = np.round(out["trend_10y"], 4)
    tf = np.round(out["trend_full"], 4)
    rec = cat == "receding"
    sign_flip = float(np.mean((t10 >= 0) != (tf >= 0)))
    rec_full_neg = int((rec & (tf < 0)).sum())
    appr_safe_full_pos = int(((cat == "approaching") & ~orbit_ok & (tf >= 0)).sum())
    n_flyby = len(fo)
    n_obj_fly = int((out["n_ca"] > 0).sum())
    n_more = int((out["n_ca"] > N_CA).sum())
    strict_early = int((strict & (out["min_jd"] <= t10_end)).sum())
    at_start = (out["min_jd"] <= START_JD + 1e-6) & (out["min_dist"] < CA_AU)
    at_end = (out["min_jd"] >= END_JD - 1e-6) & (out["min_dist"] < CA_AU)
    n_edge = int((at_start | at_end).sum())
    n_edge_imp = int((impact_capable & (at_start | at_end)).sum())
    n_edge_nolist = int(((at_start | at_end) & (out["n_ca"] == 0)).sum())

    # ---- top 20 individual approaches in 2037..2100 (best per object)
    late = fj >= 2465059.5                                           # 2037-01-01
    li = np.flatnonzero(late)
    best = {}
    for q in li:
        k = int(fo[q])
        if k not in best or fd[q] < fd[best[k]]:
            best[k] = q
    late_top = sorted(best.values(), key=lambda q: fd[q])[:20]
    n_late_obj = len(best)
    apo = next((k for k, o in enumerate(objs) if "Apophis" in o["name"]), None)
    apo_line, apo68 = "", ""
    if apo is not None:
        apo_line = (f"Check case: Apophis (99942) minimum in this model {out['min_dist'][apo]:.6f} AU on "
                    f"{jd_to_iso(out['min_jd'][apo])[:10]} versus the real {0.000254:.6f} AU on 2029-04-13 (JPL, perturbed N-body).")
        t68, d68, _ = refine_candidates(prepare_elements([objs[apo]]), np.array([0]), np.array([iso_to_jd("2068-04-12")]))
        apo68 = f"The unperturbed model has its own 2068 pass on {jd_to_iso(t68[0])[:10]} at {d68[0]:.3f} AU (> {CA_AU} AU, so not listed as a close approach)."
    seam = ""
    rng_ = ephemeris_range()
    if rng_ is not None and rng_[1] < END_JD:
        tt = np.linspace(rng_[0], rng_[1], 4000)
        dif = np.linalg.norm(earth_position(tt) - earth_position_standish(tt), axis=-1)
        seam = (f" Over the DE441 span Standish differs from DE441 by at most {dif.max():.1e} AU (rms {np.sqrt((dif ** 2).mean()):.1e} AU; mostly the Earth-to-EM-barycentre offset), "
                f"and the hand-over at {jd_to_iso(rng_[1])[:10]} jumps by {np.linalg.norm(earth_position(rng_[1] - 1e-6) - earth_position(rng_[1] + 1e-6)):.1e} AU.")

    rep = f"""# Trajectory report

Generated by `python -I scripts/trajectory.py` from `data/neo_elements.json` ({n:,} NEOs, JPL SBDB). Runtime {runtime_s or 0:.0f} s.

## Method
- Two-body Keplerian propagation of each asteroid from its catalogue elements, fetched with `full-prec=true` (full double precision; Halley Kepler solver, mean motion from `a` via the Gaussian constant). The coarse scan marches the 5-day grid with a Taylor-predicted eccentric anomaly (2 Halley iterations per sample) and evaluates the geocentric distance as |r|^2 + |R|^2 - 2 r.R on cache-sized blocks of 5,000 objects.
- **Earth: {earth_source()}.** DE441 = JPL Horizons Earth-centre minus Sun state vectors from `app/public/data/ephemeris/planets.bin`, cubic-Hermite interpolated (true geocentre, no barycentre offset). Outside the file's range the code uses the Standish (JPL SSD "Approximate Positions of the Major Planets", Table 1, valid 1800-2050) mean elements of the Earth-Moon barycentre, ~1e-4 AU formal error inside 1800-2050 and slowly growing when extrapolated toward 2100; this is far below the two-body error and is not a limiting factor.{seam} The loader picks up an extended ephemeris file automatically (and re-checks that manifest and binary agree in size before using them).
- **Window: {jd_to_iso(START_JD)} (JD {START_JD}) .. {jd_to_iso(END_JD)} (JD {END_JD}), {WINDOW_DAYS:.0f} days (~74.2 yr), 5-day steps (5,424 samples including the end point).** Every coarse local minimum of the geocentric distance below {REFINE_CAND_AU} AU (plus the global minimum and the two window-edge samples if they are minima) is refined with three 21-point grids (0.5 d, 0.06 d, 0.006 d spacing) and a parabola through d^2, so minimum distances/dates are not limited by the 5-day sampling. Per-object minimum over candidates uses `np.minimum.at` with the matching argmin time.
- **Close approaches** (`next_close_approach`, up to {N_CA}, soonest first, ISO date + distance in AU): every refined local minimum with d < {CA_AU} AU strictly inside the window (events closer than {DEDUPE_DAYS:.0f} d are merged). `n_close_approaches` is the total count for the object and is only present when the list was truncated (more than {N_CA}); otherwise the count is simply the list length. Window-edge minima (a flyby that already happened just before the start) are not listed.
- Per object: `d_start`, `d(+10 yr)`, `d_end`, `trend_au_10y = d(2036-10-07) - d(start)`, `trend_au_full = d(2100-12-31) - d(start)` (both rounded to 4 dp in `asteroids.json`; `trend_au` is kept as an alias of `trend_au_10y` for older consumers), `min_dist` and `min_date`.

## Category rules (exactly one per object) - unchanged, 74-year window
1. **impact** (can crash) = (MOID <= {IMPACT_AU} AU, or flagged PHA, or MOID unknown) AND predicted minimum geocentric distance anywhere in 2026-10-09..2100-12-31 < {IMPACT_AU} AU. Independent of the trend.
2. **receding** (never crashes) = NOT impact-capable by orbit (MOID > {IMPACT_AU} AU and not PHA) AND `trend_au_10y >= 0`. A PHA or MOID <= {IMPACT_AU} object is never `receding`.
3. **approaching** = everything else (distance diminishing but safe, or low-MOID/PHA objects that do not come within {IMPACT_AU} AU in the window, regardless of trend).
4. `strict_collision` flag (subset of `impact`) = min_dist < {STRICT_AU} AU. It is still a two-body, unfocused, unperturbed prediction, not an impact probability.

### Re-check of the rules for a 74-year window
- **Trend must stay a 10-year quantity.** `trend_au_full = d(2100-12-31) - d(start)` compares two points 74 years apart (~74 orbital phases of a 1-year body); for {100 * sign_flip:.1f}% of the objects its sign differs from `trend_au_10y`, so it says almost nothing about whether an object is closing or receding. Receding/approaching therefore uses `trend_au_10y` (the same quantity as the previous 10-year release); `trend_au_full` is shipped for reference only. Because `receding` already requires MOID > {IMPACT_AU} AU and not PHA, the receding set does not depend on the window length at all: {counts['receding']:,} objects, identical in definition to the previous release ({rec_full_neg:,} of them would have a negative full-window trend, which shows how arbitrary that signal is).
- **The MOID/PHA gate is nearly redundant, as it should be.** Objects outside the gate (MOID > {IMPACT_AU} AU, not PHA) that nevertheless reach d < {IMPACT_AU} AU in the window: **{outside_hit}** (expected ~0, since min distance >= MOID up to ~1e-4 AU of Earth-orbit variation). The gate stays as a cheap consistency guard and for objects with unknown MOID.
- **The `impact` label is now mostly an orbit-geometry statement.** {n_ok:,} objects are impact-capable by orbit (MOID <= {IMPACT_AU} AU, PHA or unknown MOID); {n_imp:,} of them ({100 * n_imp / n_ok:.1f}%) pass within {IMPACT_AU} AU at least once in the 74-year window, versus {imp10:,} ({100 * imp10 / n_ok:.1f}%) within the first 10 years only (the previous window; the previous release reported 2,309). Over 74 years roughly half of all impact-capable orbits get a conjunction near a node crossing (the per-year flyby rate is flat after the first decade, see the table below), so `impact` is better read together with `min_dist_au`, `min_date`, `strict_collision` and `next_close_approach` (the UI shows those) than as "will hit". Rules were kept as requested; if a tighter notion is wanted, thresholds on `min_dist_au` (e.g. < 0.01 AU) or `strict_collision` are the knobs.
- **Window-edge effect.** {n_edge} objects have their window minimum on the first sample (2026-10-09: a flyby already in progress or just past, {int(at_start.sum())}) or on the last one (2100-12-31: an approach that peaks after the window, {int(at_end.sum())}); {n_edge_imp} of them are `impact`; the edge minimum itself is not an interior future flyby, so it is never listed in `next_close_approach` ({n_edge_nolist} of these objects have no other close approach and therefore no list at all). This is inherited from the previous release's definition (the window includes its end samples).
- Beyond ~2040 the two-body flyby *dates* are not predictions (see caveats); the list is an unperturbed-geometry screening.

## Counts
| category | count | share |
|---|---:|---:|
| impact | {counts['impact']:,} | {100*counts['impact']/n:.1f}% |
| approaching | {counts['approaching']:,} | {100*counts['approaching']/n:.1f}% |
| receding | {counts['receding']:,} | {100*counts['receding']/n:.1f}% |
| **total** | **{n:,}** | |

- strict_collision subset of `impact` (min_dist < {STRICT_AU} AU): **{n_strict:,}** ({strict_early:,} reach their minimum in the first 10 years, {n_strict - strict_early:,} later; the previous 10-year release had 2, this run restricted to the first 10 years gives {strict10:,})
- impact-capable (= impact count): {n_imp:,}; the same rule on the first 10 years only: {imp10:,}
- PHA-flagged objects: {int(pha.sum()):,}; PHAs in `receding`: **{int(np.sum(pha & rec))}**; objects with MOID <= {IMPACT_AU} AU in `receding`: **{int(np.sum(rec & (np.nan_to_num(moid, nan=0.0) <= IMPACT_AU)))}**
- Close approaches < {CA_AU} AU inside the window: **{n_flyby:,}** events on **{n_obj_fly:,}** objects ({n_more:,} objects have more than {N_CA}, only the first {N_CA} are shipped; `n_close_approaches` has the count).

## Close approaches by period
{_period_table(out)}

## Top 20 closest approaches in the window (2026-10-09 .. 2100-12-31)
| # | name | category | min dist (AU) | min dist (km) | date | MOID (AU) | PHA |
|--:|---|---|---:|---:|---|---:|:-:|
"""
    for r_, k in enumerate(order, 1):
        o = objs[k]
        mo = "n/a" if o.get("moid") is None else f"{o['moid']:g}"
        rep += (f"| {r_} | {o['name'].strip()} | {cat[k]} | {out['min_dist'][k]:.5f} | "
                f"{out['min_dist'][k]*AU_KM:,.0f} | {jd_to_iso(out['min_jd'][k])[:10]} | {mo} | {'Y' if o['pha'] else ''} |\n")
    rep += f"""
## Top 20 closest approaches in 2037-2100 (best flyby per object; {n_late_obj:,} objects have one < {CA_AU} AU after 2036)
Unperturbed two-body geometry; dates are indicative only.

| # | name | category | dist (AU) | dist (km) | date | MOID (AU) | PHA | diameter (km) |
|--:|---|---|---:|---:|---|---:|:-:|---:|
"""
    for r_, q in enumerate(late_top, 1):
        k = int(fo[q])
        o = objs[k]
        mo = "n/a" if o.get("moid") is None else f"{o['moid']:g}"
        dia = o.get("diameter_km")
        dia = estimate_diameter(o["H"]) if dia is None and o.get("H") is not None else dia
        rep += (f"| {r_} | {o['name'].strip()} | {cat[k]} | {fd[q]:.5f} | {fd[q]*AU_KM:,.0f} | {jd_to_iso(fj[q])[:10]} | {mo} | "
                f"{'Y' if o['pha'] else ''} | {'' if dia is None else f'{dia:.3g}'} |\n")
    rep += f"""
## Accuracy caveats
- Pure two-body model: no planetary perturbations, no Earth gravitational focusing (matters below ~0.001 AU), no Yarkovsky. Treat the predicted miss distances as screening-grade. {apo_line}
- **Error growth over 74 years.** The osculating semi-major axis differs from the long-term mean by ~1e-4..1e-3 relative for typical NEOs, so the mean-motion error is ~1.5e-4..1.5e-3 and a 1-year-period body accumulates up to ~0.07..0.7 rad of orbital phase by 2100 (more for shorter periods); a close planetary encounter changes the orbit outright. The *existence* of flybys follows orbit geometry (MOID, node phasing); their *dates and distances* after ~2040 are not forecasts. Example: after Apophis' 2029-04-13 pass (real 0.000254 AU) JPL's N-body orbit (SBDB close-approach API, fetched 2026-10-09) has a ~1.1 AU semi-major axis and no approach < 2 AU between 2066-09-17 and 2073-02-14, while the unperturbed pre-2029 orbit keeps its 0.92 AU orbit. {apo68} Its sub-0.05 AU passes (2060, 2091) also differ from JPL's (2051-04-20 at 0.041 AU, 2066-09-16 at 0.069 AU); see `scripts/test_trajectory.py`.
- Elements are full precision (`full-prec=true`), so catalogue rounding is no longer a limiting error.
- Earth position source: {earth_source()}. Fallback Standish accuracy is ~1e-4 AU (EM barycentre, includes the ~3e-5 AU barycentre-to-geocentre offset).
- Output size: `app/public/data/asteroids.json` {'' if size_mb is None else f'{size_mb:.1f} MB'}.
"""
    (ROOT / "data" / "trajectory_report.md").write_text(rep, encoding="utf-8")
    print("wrote data/trajectory_report.md", flush=True)


if __name__ == "__main__":
    sys.exit(main())
