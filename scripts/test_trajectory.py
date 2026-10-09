"""Tests for scripts/trajectory.py.   Run:  python -I scripts/test_trajectory.py   (exit 1 on failure)

Window under test: 2026-10-09 .. 2100-12-31 (27111 d), 5-day coarse steps + refinement."""
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))      # needed because `python -I` drops the script dir from sys.path
import numpy as np  # noqa: E402

import trajectory as T  # noqa: E402

ROOT = HERE.parent
fails = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        fails.append(msg)


def closest_approach(obj):
    """Run the production pipeline (coarse scan + refinement) on one object."""
    el, tgrid, out = T.analyse([obj], verbose=False)
    return out["min_dist"][0], out["min_jd"][0]


def flybys(out, k):
    """(jd array, dist array) of the close approaches (< 0.05 AU) of object k in an analyse() result."""
    s, e = np.searchsorted(out["fly_obj"], [k, k + 1])
    return out["fly_jd"][s:e], out["fly_d"][s:e]


# 0. window constants -------------------------------------------------------------------
check(T.jd_to_iso(T.START_JD) == "2026-10-09T00:00:00Z", "window starts 2026-10-09")
check(T.jd_to_iso(T.END_JD) == "2100-12-31T00:00:00Z", "window ends 2100-12-31")
check(T.WINDOW_DAYS == 27111.0 and T.STEP_DAYS == 5.0, "window 27111 d, step 5 d")
tg = T.make_grid()
check(len(tg) == 5424 and tg[0] == T.START_JD and tg[-1] == T.END_JD and np.all(np.diff(tg) <= 5.0 + 1e-9) and abs(tg[730] - T.START_JD - 3650) < 1e-9,
      "coarse grid: 5424 samples, 5-day steps, last = END_JD, index 730 = +3650 d")
check(abs(T.iso_to_jd("2029-04-13T21:46:00Z") - 2462240.4069) < 1e-3, "iso_to_jd round-trips a known date")

# 1. Kepler solver residual -------------------------------------------------------------
rng = np.random.default_rng(1)
e = np.concatenate([rng.uniform(0, 0.99, 200000), np.full(1000, 0.999), np.zeros(1000)])
M = rng.uniform(-10, 10, e.size)
E = T.solve_kepler(M, e)
res = np.abs(np.mod(E - e * np.sin(E) - M + np.pi, 2 * np.pi) - np.pi)
check(res.max() < 1e-11, f"Kepler residual max {res.max():.2e} < 1e-11 (200k random, e up to 0.999)")

# 2. circular-orbit sanity --------------------------------------------------------------
o = {"a": 1.5, "e": 0.0, "i": 30.0, "om": 40.0, "w": 10.0, "ma": 25.0, "epoch": T.START_JD}
el = T.prepare_elements([o])
period = 2 * math.pi / (T.K_GAUSS / 1.5 ** 1.5)
jd = T.START_JD + np.linspace(0, period, 721)
pos = T.asteroid_positions(el, jd)[0]
r = np.linalg.norm(pos, axis=1)
check(np.ptp(r) < 1e-12 and abs(r[0] - 1.5) < 1e-12, f"circular orbit radius constant = a (spread {np.ptp(r):.1e})")
check(np.linalg.norm(pos[-1] - pos[0]) < 1e-9, "circular orbit closes after one Gaussian period")
check(abs(np.abs(pos[:, 2]).max() - 1.5 * math.sin(math.radians(30))) < 2e-3, "max |z| = a sin(i) for i=30 deg")
ang = np.arccos(np.clip((pos[0] @ pos[180]) / (r[0] * r[180]), -1, 1))
check(abs(ang - math.pi / 2) < 1e-9, "quarter period -> 90 deg of arc")
P, Q = el["P"][0], el["Q"][0]
check(abs(P @ Q) < 1e-14 and abs(np.linalg.norm(P) - 1) < 1e-14, "orientation vectors orthonormal")

# 3. Earth ------------------------------------------------------------------------------
tj = np.linspace(T.J2000, T.J2000 + 400, 400001)
ep = T.earth_position(tj)
lon = np.unwrap(np.arctan2(ep[:, 1], ep[:, 0]))
t_period = np.interp(lon[0] + 2 * math.pi, lon, tj) - tj[0]
check(abs(t_period - 365.256) < 0.05, f"Earth sidereal period {t_period:.3f} d (expect ~365.256)")
rr = np.linalg.norm(T.earth_position(T.J2000 + np.linspace(0, 365.25, 3000)), axis=1)
check(0.9832 < rr.min() < 0.9835 and 1.0165 < rr.max() < 1.0169, f"Earth r in [{rr.min():.4f}, {rr.max():.4f}] AU (perihelion 0.9833, aphelion 1.0167)")
e0 = T.earth_position(T.J2000)
lon0 = math.degrees(math.atan2(e0[1], e0[0])) % 360
check(abs(lon0 - 100.38) < 0.2, f"Earth heliocentric longitude at J2000 = {lon0:.2f} deg (expect ~100.4)")
check(abs(T.earth_position(T.J2000)[2]) < 1e-6, "Earth stays in the ecliptic (|z| < 1e-6 AU)")
sp = T.earth_position(2459580.5)  # 2022-01-04 ~ perihelion
check(np.linalg.norm(sp) < 0.9834, "Earth near perihelion at early January 2022")
print("     Earth source:", T.earth_source())

# 3b. DE441 Earth vs Standish fallback, hand-over and the 2037..2100 tail ----------------
rng_eph = T.ephemeris_range()
tw_all = T.START_JD + np.linspace(0, T.WINDOW_DAYS, 3000)
rr = np.linalg.norm(T.earth_position(tw_all), axis=1)
check(0.982 < rr.min() and rr.max() < 1.018, f"Earth distance from Sun over the whole window in [{rr.min():.4f}, {rr.max():.4f}] AU")
if rng_eph is not None:
    lo, hi = rng_eph
    tw = np.linspace(max(lo, T.START_JD), min(hi, T.END_JD), 800)
    diff = np.linalg.norm(T.earth_position(tw) - T.earth_position_standish(tw), axis=-1)
    check(diff.max() < 3e-4, f"DE441 Earth vs Standish EMB over the DE441 span: max diff {diff.max():.2e} AU (< 3e-4)")
    t0 = T.START_JD + 100.3
    check(np.linalg.norm(T.earth_position(t0) - T.earth_position(t0 + 1e-6)) < 1e-6, "DE441 interpolation is continuous")
    full = T._load_ephemeris()
    if hi < T.END_JD:
        short = full
    else:                                   # file covers the window: emulate the old 2036-12-31 end to exercise the hand-over
        short = dict(full)
        short["count"] = 6210
        short["rel"] = full["rel"][:6210]
        T._EPH["v"] = short
    try:
        hi_s = short["jd0"] + (short["count"] - 1) * short["step"]
        jump = np.linalg.norm(T.earth_position(hi_s - 1e-6) - T.earth_position(hi_s + 1e-6))
        check(jump < 1e-4, f"hand-over DE441 -> Standish at {T.jd_to_iso(hi_s)[:10]} jumps by {jump:.2e} AU (< 1e-4)")
        tail = np.linspace(hi_s + 1, T.END_JD, 400)
        check(np.allclose(T.earth_position(tail), T.earth_position_standish(tail), rtol=0, atol=1e-12), "beyond the ephemeris range the Standish elements are used")
        mixed = np.concatenate([np.linspace(hi_s - 50, hi_s, 5), np.linspace(hi_s + 1, hi_s + 50, 5)])
        check(np.linalg.norm(T.earth_position(mixed) - T.earth_position_standish(mixed), axis=-1).max() < 3e-4 and T.earth_position(mixed).shape == (10, 3),
              "arrays straddling the hand-over are evaluated piecewise")
    finally:
        T._EPH["v"] = full
    if hi >= T.END_JD:
        tail = np.linspace(T.START_JD, T.END_JD, 400)
        d_tail = np.linalg.norm(T.earth_position(tail) - T.earth_position_standish(tail), axis=-1)
        check(d_tail.max() < 5e-4, f"ephemeris covers the whole window: DE441 vs Standish max {d_tail.max():.2e} AU over 2026-2100 (< 5e-4)")
else:
    print("     (no ephemeris files: Standish fallback in use)")

# 4. coarse scan == generic path --------------------------------------------------------
els = json.load(open(ROOT / "data" / "neo_elements.json"))
rng2 = np.random.default_rng(5)
sample = [els[i] for i in rng2.choice(len(els), 250, replace=False)]
sample += sorted(els, key=lambda x: -x["e"])[:10] + sorted(els, key=lambda x: x["a"])[:10] + sorted(els, key=lambda x: x["epoch"])[:5]
elS = T.prepare_elements(sample)
K = len(tg)
earthg = T.earth_position(tg)
scan = T.scan_block(elS, tg, earthg, 730)
ks = np.array([0, 1, 2, 365, 730, 2711, K - 2, K - 1])
dg = T.geo_distance(elS, np.broadcast_to(tg[ks], (len(sample), len(ks))))
dscan = {0: scan["d_start"], 4: scan["d_10"], len(ks) - 1: scan["d_end"]}
err = max(np.abs(dscan[j] - dg[:, j]).max() for j in dscan)
check(err < 1e-9, f"streaming coarse scan matches generic Kepler/ephemeris path at d_start, d(+10y), d_end (max |diff| {err:.1e} AU, {len(sample)} objects incl. extreme e, a, epoch)")
idx = rng2.choice(len(sample), 12, replace=False)
dfull = np.stack([T.geo_distance(T.sub_elements(elS, np.array([i])), tg[None, :])[0] for i in idx])
check(np.abs(dfull.min(1) - scan["coarse_min"][idx]).max() < 1e-9, "coarse global minimum equals brute-force minimum over the 5424-point grid")

# 5. Apophis ----------------------------------------------------------------------------
# (a) full-precision JPL SBDB osculating elements for 99942 Apophis, epoch JD 2461200.5
#     (fetched from ssd-api.jpl.nasa.gov/sbdb.api?sstr=99942&full-prec=1 on 2026-10-09)
apo_full = {"a": 0.9223592206975018, "e": 0.1911492279663492, "i": 3.340996879880978,
            "om": 203.8936514240762, "w": 126.6795706895841, "ma": 175.3304026592739, "epoch": 2461200.5}
d, jd_min = closest_approach(apo_full)
date = T.jd_to_iso(jd_min)
print(f"     Apophis (full precision): {d:.6f} AU = {d * T.AU_KM:,.0f} km at {date}")
# Two-body (no Earth focusing) closest approach is a few 1e-4 AU, vs the real (focused, perturbed)
# 0.000254 AU.  Two-body propagation from a 2026 epoch to 2029 also loses accuracy to planetary
# perturbations, so assert a sensible physical band (well inside the Moon distance of 0.00257 AU,
# above Earth radius 4.3e-5 AU) rather than an exact value; the date is checked to +-3 d.
check(date.startswith("2029-04"), f"Apophis closest approach date {date[:10]} is in 2029-04 (still the global minimum over 2026-2100)")
check(abs(jd_min - 2462240.4) < 3.0, f"Apophis closest approach within 3 d of 2029-04-13 21:46 UT (got {jd_min - 2462240.4:+.2f} d)")
check(0.0001 < d < 0.0008, f"Apophis min distance {d:.6f} AU in [0.0001, 0.0008] (real 0.000254; two-body from a 2026 epoch is perturbation-limited, got ~0.00014)")
apo = next(x for x in els if "Apophis" in x["name"])
d2, jd2 = closest_approach(apo)
print(f"     Apophis (catalogue elements): {d2:.6f} AU at {T.jd_to_iso(jd2)}")
check(T.jd_to_iso(jd2).startswith("2029-04"), "Apophis (catalogue elements) in 2029-04")
check(0.0001 < d2 < 0.0008, f"Apophis (catalogue elements) min {d2:.6f} AU in [0.0001, 0.0008]")
el1, tg1, out1 = T.analyse([apo_full], verbose=False)
check(out1["min_dist"][0] <= out1["coarse_min"][0] + 1e-12 and out1["coarse_min"][0] - out1["min_dist"][0] > 1e-4,
      f"refinement improves coarse 5-day minimum ({out1['coarse_min'][0]:.5f} -> {out1['min_dist'][0]:.5f} AU)")
# (c) all close approaches through 2100 and the 2068 flyby.
#     JPL SBDB close-approach API (orbit solution 220, queried 2026-10-09): Apophis 2029-04-13 21:46 at
#     0.000254 AU, then 2044-08-24 0.080, 2051-04-20 0.0415, 2066-09-16 0.0694 AU, and NO approach
#     < 2 AU between 2066-09-17 and 2073-02-14 (i.e. none in 2068): the 2029 encounter changes the
#     orbit (a: 0.92 -> ~1.1 AU).  An unperturbed two-body model cannot know that, so it keeps the
#     pre-encounter 0.886-yr orbit and reproduces the *old* 2068-04-12 geometry.
fj, fd = flybys(out1, 0)
print("     Apophis model close approaches < 0.05 AU:", [(T.jd_to_iso(a)[:10], round(float(b), 4)) for a, b in zip(fj, fd)])
check(len(fj) >= 1 and T.jd_to_iso(fj[0]).startswith("2029-04") and abs(fd[0] - d) < 1e-9, "Apophis: first listed close approach is the 2029-04 pass and equals the global minimum")
check(bool(np.all(np.diff(fj) > 0)) and bool(np.all(fd < T.CA_AU)) and bool(np.all((fj > T.START_JD) & (fj < T.END_JD))), "Apophis flybys: chronological, all < 0.05 AU, inside the window")
t68, d68, ok68 = T.refine_candidates(el1, np.array([0]), np.array([T.iso_to_jd("2068-04-12")]))
tt = T.iso_to_jd("2067-06-01") + np.arange(0, 730, 0.05)
dd68 = T.geo_distance(el1, tt[None, :])[0]
print(f"     Apophis 2068: model {d68[0]:.4f} AU on {T.jd_to_iso(t68[0])[:10]}  |  JPL: nothing within 2 AU (2066-09-17..2073-02-14)")
check(abs(d68[0] - dd68.min()) < 1e-6 and abs(t68[0] - tt[dd68.argmin()]) < 0.1, "pipeline refinement of the 2068 pass equals a brute-force 0.05-day scan")
check(T.jd_to_iso(t68[0]).startswith("2068-04") and 0.05 < d68[0] < 0.3,
      f"model 2068 pass ({d68[0]:.3f} AU, {T.jd_to_iso(t68[0])[:10]}) is an unperturbed-geometry artefact: > 0.05 AU so not listed, JPL has no 2068 flyby < 2 AU")
check(not any(T.jd_to_iso(a).startswith("2068") for a in fj), "Apophis: no 2068 entry in the close-approach list")
check(all(abs(np.min(np.abs(fj - T.iso_to_jd(s))) - 0) > 30 for s in ("2044-08-24", "2051-04-20", "2066-09-16")),
      "model flybys after 2029 do not coincide with JPL's (2044-08-24, 2051-04-20, 2066-09-16): post-encounter orbit is not modelled")

# 6. close-approach detection: synthetic synodic geometry + brute-force completeness ----
# (a) circular, coplanar orbit with a = 1.03 AU: conjunctions every synodic period (~23 yr), each closer
#     than 0.05 AU because |a - r_Earth| in [0.0133, 0.0467].
circ = {"a": 1.03, "e": 0.0, "i": 0.0, "om": 0.0, "w": 0.0, "ma": 10.0, "epoch": T.START_JD}
elc, tgc, outc = T.analyse([circ], verbose=False)
fjc, fdc = flybys(outc, 0)
n_ = T.K_GAUSS / 1.03 ** 1.5
td = T.START_JD + np.arange(0, T.WINDOW_DAYS, 1.0)
lon_a = np.radians(10.0) + n_ * (td - T.START_JD)
pe = T.earth_position(td)
dl = np.mod(lon_a - np.arctan2(pe[:, 1], pe[:, 0]) + np.pi, 2 * np.pi) - np.pi
conj = td[:-1][(dl[:-1] < 0) & (dl[1:] >= 0)] if n_ > 2 * math.pi / 365.256 else td[:-1][(dl[:-1] > 0) & (dl[1:] <= 0)]
print(f"     synodic test: equal-longitude dates {[T.jd_to_iso(c)[:10] for c in conj]}, pipeline minima {[T.jd_to_iso(c)[:10] for c in fjc]}")
check(len(fjc) == len(conj) >= 3, f"circular a=1.03 AU orbit: {len(fjc)} close approaches == {len(conj)} geometric conjunctions in 74 yr")
# independent analytic reference: circular orbit position a*(cos L, sin L, 0); distance minimum searched at 0.01 d around each conjunction.
# (The distance minimum trails equal longitude by ~2-3 weeks: Earth's radial velocity shifts it, so compare against the analytic d(t), not the longitude zero.)
ref_t, ref_d = [], []
for c in conj:
    tt = c + np.arange(-60.0, 60.0, 0.01)
    L = np.radians(10.0) + n_ * (tt - T.START_JD)
    dd_ = np.linalg.norm(1.03 * np.stack([np.cos(L), np.sin(L), 0 * L], axis=-1) - T.earth_position(tt), axis=-1)
    ref_t.append(tt[dd_.argmin()]); ref_d.append(dd_.min())
ref_t, ref_d = np.array(ref_t), np.array(ref_d)
check(len(fjc) == len(conj) and bool(np.all(np.abs(fjc - ref_t) < 0.02)), "flyby dates equal the analytic circular-orbit distance minima to < 0.02 d (and trail the equal-longitude instant by weeks)")
check(len(fjc) == len(conj) and bool(np.all(np.abs(fdc - ref_d) < 1e-6)), f"flyby distances equal the analytic circular-orbit minima to < 1e-6 AU ({', '.join(f'{x:.4f}' for x in fdc)})")
# (b) brute-force 0.5-day scan vs pipeline for 30 low-MOID objects (tens of flybys in total)
low = [x for x in els if x["moid"] is not None and x["moid"] < 0.02]
rng3 = np.random.default_rng(11)
pick = [low[i] for i in rng3.choice(len(low), 30, replace=False)]
elp, tgp, outp = T.analyse(pick, verbose=False)
step = 0.5
tdn = T.START_JD + np.arange(0, T.WINDOW_DAYS + step, step)
tdn = tdn[tdn <= T.END_JD]
n_dense = missed = extra = 0
for k in range(len(pick)):
    dk = T.geo_distance(T.sub_elements(elp, np.array([k])), tdn[None, :])[0]
    lm = np.flatnonzero((dk[1:-1] <= dk[:-2]) & (dk[1:-1] <= dk[2:]) & (dk[1:-1] < T.CA_AU)) + 1
    pt, pdist = flybys(outp, k)
    n_dense += len(lm)
    for t_ in tdn[lm]:
        if not len(pt) or np.min(np.abs(pt - t_)) > 1.5:
            missed += 1
    for t_ in pt:
        if not len(lm) or np.min(np.abs(tdn[lm] - t_)) > 1.5:
            extra += 1
check(n_dense >= 20 and missed == 0 and extra == 0,
      f"pipeline close approaches == brute-force 0.5-day scan on 30 low-MOID objects ({n_dense} events, {missed} missed, {extra} extra)")

# 7. category rules (unit test of categorise on hand-made cases) ----------------------
cases = [  # (moid, pha, min_dist, trend_10y, expected)
    (0.01, False, 0.01, 0.5, "impact"),          # low MOID and comes within 0.05
    (0.01, False, 0.20, 0.5, "approaching"),     # low MOID never receding even with trend >= 0
    (0.20, False, 0.30, 0.5, "receding"),        # safe orbit, distance not decreasing over 10 yr
    (0.20, False, 0.30, -0.5, "approaching"),    # safe orbit but closing
    (0.20, False, 0.30, 0.0, "receding"),        # trend exactly 0 counts as receding
    (0.20, True, 0.02, 0.5, "impact"),           # PHA flag overrides MOID gate
    (0.20, True, 0.30, 0.5, "approaching"),      # PHA never receding
    (None, False, 0.04, 0.5, "impact"),          # unknown MOID counts as impact-capable
    (None, False, 0.30, 0.5, "approaching"),     # ... and is never receding
    (0.01, False, 0.05, 0.5, "approaching"),     # min_dist must be strictly < 0.05
    (0.05, False, 0.049, -1.0, "impact"),        # MOID == 0.05 is inside the gate
    (0.0004, False, 0.0005, 0.0, "impact"),      # strict candidate
]
fake = [{"moid": c[0], "pha": c[1]} for c in cases]
fo_ = {"min_dist": np.array([c[2] for c in cases]), "trend_10y": np.array([c[3] for c in cases]),
       "trend_full": -np.array([c[3] for c in cases]) - 1.0}     # full trend deliberately contradicts: must be ignored
cat_, imp_, strict_, _ = T.categorise(fake, fo_)
check(list(cat_) == [c[4] for c in cases], f"categorise() matches the documented rules on {len(cases)} hand-made cases (trend_full ignored)")
check(bool(strict_[-1]) and int(strict_.sum()) == 1, "strict_collision only for impact objects with min_dist < 0.001 AU")

# 8. MOID -------------------------------------------------------------------------------
bennu = next(x for x in els if "Bennu" in x["name"])
moid = T.orbit_moid(bennu["a"], bennu["e"], bennu["i"], bennu["om"], bennu["w"], jd_earth=T.J2000 + 30 * 365.25)
print(f"     Bennu MOID computed {moid:.5f} AU (JPL: {bennu['moid']})")
check(abs(moid - 0.00322) < 0.0005 and abs(moid - bennu["moid"]) < 0.0005, f"Bennu MOID {moid:.5f} AU ~ JPL 0.00322 (tol 5e-4)")
moid_apo = T.orbit_moid(apo_full["a"], apo_full["e"], apo_full["i"], apo_full["om"], apo_full["w"], jd_earth=T.START_JD)
check(abs(moid_apo - 0.000108) < 0.0001, f"Apophis MOID {moid_apo:.6f} AU ~ JPL 0.000108 (tol 1e-4)")

# 9. Output files -----------------------------------------------------------------------
pub = ROOT / "app" / "public" / "data" / "asteroids.json"
rows = json.load(open(pub))
size_mb = pub.stat().st_size / 1e6
need = ["id", "name", "diameter_km", "pha", "category", "min_dist_au", "min_date", "a", "e", "i", "om", "w", "ma", "epoch"]
extra_f = ["strict_collision", "trend_au", "trend_au_10y", "trend_au_full"]
check(len(rows) == len(els), f"asteroids.json has {len(rows):,} rows == neo_elements ({len(els):,})")
check(size_mb < 16.0, f"asteroids.json is {size_mb:.2f} MB (< 16 MB)")
check(all(all(k in r and r[k] is not None for k in need + extra_f) for r in rows), "all contract + trend/strict fields present and non-null")
check(all(r["category"] in ("impact", "approaching", "receding") for r in rows), "every row has one of exactly 3 categories")
check(len({r["id"] for r in rows}) == len(rows), "ids unique")
check(all(r["diameter_km"] > 0 for r in rows), "all diameters positive (null estimated)")
check(all(r["min_date"][:2] in ("20", "21") and r["min_date"].endswith("Z") for r in rows), "min_date ISO strings")
lo_iso, hi_iso = T.jd_to_iso(T.START_JD), T.jd_to_iso(T.END_JD)
check(all(lo_iso <= r["min_date"] <= hi_iso for r in rows), "min_date within 2026-10-09 .. 2100-12-31")
check(any(r["min_date"] > "2037" for r in rows) and any(r["min_date"] > "2090" for r in rows), "min_dates extend beyond the old 10-year window and into the 2090s")
check(all(0 <= r["min_dist_au"] for r in rows), "min_dist_au >= 0")
check(all(r["trend_au"] == r["trend_au_10y"] for r in rows), "trend_au is an alias of trend_au_10y")
imp = [r for r in rows if r["category"] == "impact"]
check(all(r["min_dist_au"] < T.IMPACT_AU + 1e-6 for r in imp), "all impact objects have min_dist < 0.05 AU")
moid_of = {x["id"]: x["moid"] for x in els}
pha_of = {x["id"]: x["pha"] for x in els}
rec = [r for r in rows if r["category"] == "receding"]
check(not any(pha_of[r["id"]] for r in rec), "no PHA is labelled receding")
check(all(moid_of[r["id"]] is not None and moid_of[r["id"]] > T.IMPACT_AU for r in rec), "receding => MOID > 0.05 AU (known)")
check(all(r["trend_au_10y"] >= 0 for r in rec), "receding => trend_au_10y >= 0 (10-year trend, not the 74-year one)")
check(sum(1 for r in rows if r["trend_au_10y"] < 0 <= r["trend_au_full"]) > 100 and sum(1 for r in rows if r["trend_au_10y"] >= 0 > r["trend_au_full"]) > 100,
      "trend_au_10y and trend_au_full genuinely differ in sign for many objects (the full-window trend is a weak signal)")


def expected_category(r):
    m = moid_of[r["id"]]
    orbit_ok = (m is None or m <= T.IMPACT_AU) or pha_of[r["id"]]
    if orbit_ok and r["min_dist_au"] < T.IMPACT_AU:
        return "impact"
    return "receding" if (not orbit_ok and r["trend_au_10y"] >= 0) else "approaching"


edge = lambda r: abs(r["min_dist_au"] - T.IMPACT_AU) < 2e-6
check(all(r["category"] == expected_category(r) for r in rows if not edge(r)),
      "every row's category equals the documented rule applied to its own shipped fields")
check(all(r["strict_collision"] == (r["category"] == "impact" and r["min_dist_au"] < T.STRICT_AU) for r in rows if abs(r["min_dist_au"] - T.STRICT_AU) > 2e-6),
      "strict_collision <=> impact and min_dist < 0.001")
n_nca = 0
bad = []
for r in rows:
    lst = r.get("next_close_approach")
    if lst is None:
        continue
    n_nca += 1
    dts = [c["date"] for c in lst]
    ok = (1 <= len(lst) <= 3 and dts == sorted(dts) and len(set(dts)) == len(dts)
          and all(c["dist_au"] < T.CA_AU + 1e-6 and lo_iso < c["date"] < hi_iso and c["date"].endswith("Z") for c in lst))
    if "n_close_approaches" in r:
        ok = ok and r["n_close_approaches"] > 3 and len(lst) == 3
    if not ok:
        bad.append(r["id"])
check(not bad, f"next_close_approach well-formed on {n_nca:,} rows (<=3, soonest first, < 0.05 AU, inside the window)  [{len(bad)} bad]")
miss_list = [r["id"] for r in rows if r["min_dist_au"] < T.CA_AU - 1e-6 and r["min_date"] not in (lo_iso, hi_iso) and "next_close_approach" not in r]
check(len(miss_list) == 0, f"every object with min_dist < 0.05 AU (interior minimum) has a close-approach list  [{len(miss_list)} without]")
mism = []
for r in rows:
    lst = r.get("next_close_approach")
    if lst and "n_close_approaches" not in r and r["min_date"] not in (lo_iso, hi_iso):   # untruncated list contains the (interior) global minimum
        best = min(lst, key=lambda c: c["dist_au"])
        if abs(best["dist_au"] - r["min_dist_au"]) > 2e-6 or best["date"] != r["min_date"]:
            mism.append(r["id"])
check(not mism, f"untruncated lists contain the object's global minimum (same date and distance)  [{len(mism)} mismatches]")
edge_min = [r for r in rows if r["min_dist_au"] < T.CA_AU and r["min_date"] in (lo_iso, hi_iso)]
no_list_impact = [r for r in rows if r["category"] == "impact" and "next_close_approach" not in r]
check(all(r["min_date"] in (lo_iso, hi_iso) for r in no_list_impact),
      f"the only impact objects without a close-approach list have their minimum on a window edge ({len(no_list_impact)} of {len(edge_min)} edge minima < 0.05 AU: flyby in progress at the start / peaking after 2100)")
apo_row = next(r for r in rows if "Apophis" in r["name"])
check(apo_row["category"] == "impact" and apo_row["min_date"].startswith("2029-04"), f"asteroids.json: Apophis = impact, {apo_row['min_date']}")
check(apo_row["next_close_approach"][0]["date"].startswith("2029-04") and len(apo_row["next_close_approach"]) == 3,
      f"asteroids.json: Apophis list {[(c['date'][:10], c['dist_au']) for c in apo_row['next_close_approach']]}")
tr = json.load(open(ROOT / "data" / "asteroids_trajectories.json"))
check(len(tr["objects"]) == 300 and all(len(o[k]) == 180 for o in tr["objects"] for k in ("helio_orbit", "geo_window", "geo_approach")),
      "trajectories: 300 objects x 60 points (x3 coords) per polyline")
check(len(tr["earth_orbit"]) == 180 and tr["meta"]["end_iso"] == "2100-12-31T00:00:00Z" and tr["meta"]["window_days"] == T.WINDOW_DAYS,
      "trajectories meta: 74-year window, earth_orbit 60 points")
by_id = {r["id"]: r for r in rows}
check(all(by_id[o["id"]]["min_dist_au"] == o["min_dist_au"] and by_id[o["id"]]["category"] == o["category"] for o in tr["objects"]),
      "trajectory objects agree with asteroids.json (category, min_dist)")
md = [o["min_dist_au"] for o in tr["objects"]]
check(all(o["category"] == "impact" for o in tr["objects"]) and md == sorted(md), "trajectories: the 300 closest objects, ascending (all impact-category with the 74-year window)")
rep = (ROOT / "data" / "trajectory_report.md").read_text(encoding="utf-8")
check("2100-12-31" in rep and "trend_au_full" in rep and "Top 20 closest approaches" in rep, "trajectory_report.md describes the 74-year window and both trends")
if T.earth_source() not in rep:
    print("NOTE the report was generated with a different Earth source than is available now - re-run scripts/trajectory.py")

print("\n%d failure(s)" % len(fails))
sys.exit(1 if fails else 0)
