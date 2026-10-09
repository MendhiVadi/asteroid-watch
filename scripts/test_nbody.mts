// N-body engine tests.  Run:  node scripts/test_nbody.mts      (Node >= 22.18 / 24, no dependencies)
// Needs app/public/data/ephemeris/ (python -I scripts/fetch_ephemeris.py).
//
// Reference data embedded below was fetched from JPL (SBDB full-precision Apophis elements,
// Horizons barycentric vectors / geocentric minimum) - see comments.

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const root = path.resolve(here, '..');
const lib = (f: string) => pathToFileURL(path.join(root, 'app/src/lib/nbody', f)).href;

const { Ephemeris, IDX_EARTH, IDX_MOON, IDX_SUN, BODY_ORDER, N_BODIES } = await import(lib('ephemeris.ts'));
const { IAS15, StepPoly } = await import(lib('ias15.ts'));
const { makeForce } = await import(lib('dynamics.ts'));
const { elementsToState, stateToElements, propagateElementsTwoBody } = await import(lib('elements.ts'));
const { simulate, simulateStream, sampleStateAt, NBodyRangeError, SAMPLE_STRIDE, S } = await import(lib('sim.ts'));
const { AUD_TO_KMS, AU_KM } = await import(lib('constants.ts'));

// ------------------------------------------------------------------------------------------
const D = path.join(root, 'app/public/data/ephemeris');
const manifest = JSON.parse(fs.readFileSync(path.join(D, 'manifest.json'), 'utf8'));
const buf = (f: string) => {
  const b = fs.readFileSync(path.join(D, f));
  return b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength) as ArrayBuffer;
};
const eph = new Ephemeris(manifest, buf('planets.bin'), buf('moon.bin'));

let pass = 0;
let fail = 0;
const failures: string[] = [];
function check(name: string, ok: boolean, detail = ''): void {
  if (ok) pass++;
  else {
    fail++;
    failures.push(name);
  }
  console.log(`${ok ? 'PASS' : 'FAIL'}  ${name}${detail ? '  -- ' + detail : ''}`);
}
const fmt = (x: number) => x.toExponential(3);

/** Brute-force version of the look-ahead verdict on a full sample array: after jdFrom the object must
 *  first get beyond 0.05 AU of Earth and then never re-enter it.  Returns true if it 'comes back'
 *  (or never leaves). */
function bruteReturns(r: any, jdFrom: number, SS: any, STRIDE: number): boolean {
  let outside = false;
  for (let i = 0; i < r.nSamples; i++) {
    const o = i * STRIDE;
    if (r.samples[o] <= jdFrom) continue;
    const d = r.samples[o + SS.DE];
    if (d >= 0.05) outside = true;
    else if (outside) return true;
  }
  return !outside;
}

// ------------------------------------------------------------------------------------------
// 1. Two-body energy / angular momentum / return-to-start / reversibility
// ------------------------------------------------------------------------------------------
{
  const mu = 0.01720209895 ** 2;
  const acc = (_t: number, x: Float64Array, _v: Float64Array, a: Float64Array) => {
    const r = Math.hypot(x[0], x[1], x[2]);
    const f = -mu / (r * r * r);
    a[0] = f * x[0];
    a[1] = f * x[1];
    a[2] = f * x[2];
  };
  const e = 0.5;
  const x0 = [1 - e, 0, 0.1];
  const v0 = [0, Math.sqrt((mu * (1 + e)) / (1 - e)) * 0.99, 0.01];
  const energy = (x: ArrayLike<number>, v: ArrayLike<number>) =>
    0.5 * (v[0] ** 2 + v[1] ** 2 + v[2] ** 2) - mu / Math.hypot(x[0], x[1], x[2]);
  const L = (x: ArrayLike<number>, v: ArrayLike<number>) => [
    x[1] * v[2] - x[2] * v[1],
    x[2] * v[0] - x[0] * v[2],
    x[0] * v[1] - x[1] * v[0],
  ];
  const E0 = energy(x0, v0);
  const L0 = L(x0, v0);
  const ias = new IAS15(acc, 0, x0, v0, { dtInit: 0.05 });
  const a0 = -mu / (2 * E0);
  const T = 2 * Math.PI * Math.sqrt((a0 * a0 * a0) / mu);
  const N = 100;
  const tEnd = N * T;
  let maxdE = 0;
  let maxdL = 0;
  while (ias.t < tEnd - 1e-9) {
    ias.step(tEnd);
    maxdE = Math.max(maxdE, Math.abs((energy(ias.x, ias.v) - E0) / E0));
    const l = L(ias.x, ias.v);
    maxdL = Math.max(maxdL, Math.hypot(l[0] - L0[0], l[1] - L0[1], l[2] - L0[2]) / Math.hypot(L0[0], L0[1], L0[2]));
  }
  const posErr = Math.hypot(ias.x[0] - x0[0], ias.x[1] - x0[1], ias.x[2] - x0[2]);
  check('2-body: |dE/E| over 100 orbits < 1e-13', maxdE < 1e-13, `max ${fmt(maxdE)}`);
  check('2-body: |dL/L| over 100 orbits < 1e-13', maxdL < 1e-13, `max ${fmt(maxdL)}`);
  check('2-body: back at start after 100 periods (<1e-10 AU)', posErr < 1e-10, `err ${fmt(posErr)} AU, ${ias.steps} steps`);
  // reversibility
  const xf = Float64Array.from(ias.x);
  const vf = Float64Array.from(ias.v);
  const back = new IAS15(acc, tEnd, xf, vf, { dtInit: -0.05 });
  while (back.t > 1e-9) back.step(0);
  const rev = Math.hypot(back.x[0] - x0[0], back.x[1] - x0[1], back.x[2] - x0[2]);
  check('2-body: forward+backward returns to start (<1e-10 AU)', rev < 1e-10, `err ${fmt(rev)} AU`);
  // elements <-> state round trip
  const el = { a: 1.7, e: 0.43, i: 12.3, om: 77.7, w: 201.2, ma: 333.3, epoch: 0 };
  const st = elementsToState(el, mu);
  const back2 = stateToElements(st.x, st.v, mu);
  check(
    'elements -> state -> elements round trip',
    Math.abs(back2.a - el.a) < 1e-12 && Math.abs(back2.e - el.e) < 1e-12 && Math.abs(back2.i - el.i) < 1e-9 && Math.abs(back2.om - el.om) < 1e-9 && Math.abs(back2.w - el.w) < 1e-9,
  );
  // hyperbolic elements: energy and eccentricity
  const hy = { a: -2.0, e: 1.6, i: 20, om: 10, w: 30, ma: 40, epoch: 0 };
  const sh = elementsToState(hy, mu);
  const eh = stateToElements(sh.x, sh.v, mu);
  check('hyperbolic elements -> state consistent (a<0, e=1.6)', Math.abs(eh.e - 1.6) < 1e-10 && Math.abs(eh.a + 2.0) < 1e-10, `e=${eh.e.toFixed(10)} a=${eh.a.toFixed(10)}`);
}

// ------------------------------------------------------------------------------------------
// 2. Ephemeris / Earth-Moon dynamics sanity
// ------------------------------------------------------------------------------------------
{
  // Moon's geocentric acceleration from finite-differencing the ephemeris velocity vs Newton
  const gmEM = eph.gm[IDX_EARTH] + eph.gm[IDX_MOON];
  const h = 0.02;
  const t = eph.toT(2462100.3);
  const m1 = new Float64Array(6);
  const m0 = new Float64Array(6);
  const m2 = new Float64Array(6);
  eph.moonRelState(t - h, m1);
  eph.moonRelState(t, m0);
  eph.moonRelState(t + h, m2);
  const acc = [0, 1, 2].map((c) => (m2[3 + c] - m1[3 + c]) / (2 * h));
  const r = Math.hypot(m0[0], m0[1], m0[2]);
  const newton = [0, 1, 2].map((c) => (-gmEM * m0[c]) / r ** 3);
  const diff = Math.hypot(acc[0] - newton[0], acc[1] - newton[1], acc[2] - newton[2]) / Math.hypot(newton[0], newton[1], newton[2]);
  check('Moon geocentric acceleration ~ -(GM_E+GM_M) r/r^3 (solar tide < 3%)', diff < 0.03, `rel diff ${diff.toFixed(4)}`);
  check('Moon distance in 356000-407000 km', r * AU_KM > 356000 && r * AU_KM < 407000, `${(r * AU_KM).toFixed(0)} km`);

  const force = makeForce(eph);
  const E = new Float64Array(6);
  const Mo = new Float64Array(6);
  const circular = (centre: Float64Array, gm: number, rKm: number, tag: string, days: number, tol: number) => {
    const r0 = rKm / AU_KM;
    const t0 = eph.toT(2462100.3);
    eph.state(centre === E ? IDX_EARTH : IDX_MOON, t0, centre);
    // orbit in a plane tilted w.r.t. the ecliptic
    const ux = [1, 0, 0];
    const uy = [0, Math.cos(0.4), Math.sin(0.4)];
    const vc = Math.sqrt(gm / r0);
    const x = [0, 1, 2].map((c) => centre[c] + r0 * ux[c]);
    const v = [0, 1, 2].map((c) => centre[3 + c] + vc * uy[c]);
    const ias = new IAS15(force, t0, x, v, { dtInit: 0.001 });
    const tEnd = t0 + days;
    let rmin = Infinity;
    let rmax = 0;
    const B = new Float64Array(6);
    let steps = 0;
    while (ias.t < tEnd - 1e-9) {
      ias.step(tEnd, 0.02);
      eph.state(centre === E ? IDX_EARTH : IDX_MOON, ias.t, B);
      const rr = Math.hypot(ias.x[0] - B[0], ias.x[1] - B[1], ias.x[2] - B[2]);
      rmin = Math.min(rmin, rr);
      rmax = Math.max(rmax, rr);
      steps++;
    }
    const dev = Math.max(Math.abs(rmin / r0 - 1), Math.abs(rmax / r0 - 1));
    check(`${tag}: circular orbit radius stays within ${(tol * 100).toFixed(1)}%`, dev < tol, `r range [${(rmin * AU_KM).toFixed(1)}, ${(rmax * AU_KM).toFixed(1)}] km, ${steps} steps`);
  };
  circular(Mo, eph.gm[IDX_MOON], 1837.4, 'lunar orbiter (100 km alt), 2 days', 2, 0.01);
  circular(E, eph.gm[IDX_EARTH], 20000, 'Earth orbiter (20000 km), 10 days', 10, 0.01);
}

// ------------------------------------------------------------------------------------------
// 3. Apophis 2029 (SBDB full-precision elements, orbit solution 220, epoch JD 2461200.5)
// ------------------------------------------------------------------------------------------
const APOPHIS = {
  a: 0.9223592206975018,
  e: 0.1911492279663492,
  i: 3.340996879880978,
  om: 203.8936514240762,
  w: 126.6795706895841,
  ma: 175.3304026592739,
  epoch: 2461200.5,
};
// Horizons (DE441) barycentric ecliptic vectors for Apophis (AU, AU/d): [JD, x,y,z,vx,vy,vz]
const HZ = [
  [2461322.5, -0.1105618051078034, -0.8263611092945964, 0.04134264285504785, 0.01980526809405685, 0.0006577317728774277, 0.0004330524225866625],
  [2461700.5, 0.7354111389415391, -0.1591551216834056, 0.02572567193294531, 0.005520828338010832, 0.02086754755948199, -0.0009833253396573098],
  [2462239.5, -0.9246841060995419, -0.3936217653088041, -0.0008062719440310929, 0.008889261946540581, -0.01368236205380207, 0.0009635479566724925],
  [2462240.5, -0.9153932758697816, -0.4072715321456429, 0.0001890733649595224, 0.01035802614649523, -0.01462512795686298, 0.0008486476385556299],
  [2462600.5, -1.035341342453986, 0.5884754080599729, -0.03690185182331981, -0.005105245555734841, -0.01421935690976864, 0.0004266550004488913],
  [2463500.5, -1.00084953517144, -0.2617779637203684, -0.00599483333765123, 0.007371724384955634, -0.01585040141007828, 0.0006777492138790244],
  [2465000.5, 0.551028261419338, 1.135953756321344, -0.03212890503254839, -0.01187625226035125, 0.007502588374447981, -0.0004510724190142465],
];
const HZ_MIN_DIST_AU = 0.00025409077441363115; // Horizons geocentric 1-minute scan
const HZ_MIN_JD = 2462240.406944444;

{
  // (a) pure integrator vs Horizons, starting from the Horizons vector (JPL's model also has
  //     Yarkovsky A1/A2 and 16 massive asteroids which we do not model)
  const force = makeForce(eph);
  const run = (eps: number) => {
    const ias = new IAS15(force, eph.toT(HZ[0][0]), HZ[0].slice(1, 4), HZ[0].slice(4, 7), { epsilon: eps, dtInit: 0.1, maxDt: 20 });
    const errs: number[] = [];
    const posAt: number[][] = [];
    for (let k = 1; k < HZ.length; k++) {
      const tt = eph.toT(HZ[k][0]);
      while (ias.t < tt - 1e-10) ias.step(tt, 1e9);
      posAt.push(Array.from(ias.x));
      errs.push(Math.hypot(ias.x[0] - HZ[k][1], ias.x[1] - HZ[k][2], ias.x[2] - HZ[k][3]));
    }
    return { errs, steps: ias.steps, xs: Array.from(ias.x), pos: posAt };
  };
  const r9 = run(1e-9);
  const r12 = run(1e-12);
  // errs[] = 2461700 (+0.4 yr), 2462239.5 (0.9 d before), 2462240.5 (0.1 d after), 2462600.5 (+1 yr), ...
  check('Apophis vs Horizons 0.9 d BEFORE flyby (3 years of integration) < 1e-7 AU', r9.errs[1] < 1e-7, `err ${fmt(r9.errs[1])} AU (${(r9.errs[1] * AU_KM).toFixed(1)} km)`);
  check('Apophis vs Horizons 0.1 d AFTER flyby < 1e-7 AU', r9.errs[2] < 1e-7, `err ${fmt(r9.errs[2])} AU (${(r9.errs[2] * AU_KM).toFixed(1)} km)`);
  check('Apophis vs Horizons 1 yr / 10 yr after flyby < 3e-4 / 1e-3 AU (model-limited, see sensitivity)', r9.errs[3] < 3e-4 && r9.errs[5] < 1e-3, `errs ${r9.errs.map(fmt).join(' ')}`);
  const conv0 = Math.hypot(r9.pos[2][0] - r12.pos[2][0], r9.pos[2][1] - r12.pos[2][1], r9.pos[2][2] - r12.pos[2][2]);
  const conv = Math.hypot(r9.xs[0] - r12.xs[0], r9.xs[1] - r12.xs[1], r9.xs[2] - r12.xs[2]);
  check('integrator self-convergence eps 1e-9 vs 1e-12 right after the flyby < 1e-9 AU', conv0 < 1e-9, `diff ${fmt(conv0)} AU (${(conv0 * AU_KM * 1000).toFixed(1)} m), steps ${r9.steps} vs ${r12.steps}`);
  console.log(`INFO  self-convergence difference 10 years later (flyby-amplified): ${fmt(conv)} AU = ${(conv * AU_KM).toFixed(2)} km`);
  {
    // sensitivity: shift the start state by 1 km -> how big is the difference 1 yr after the flyby?
    const ias = new IAS15(force, eph.toT(HZ[0][0]), [HZ[0][1] + 1 / AU_KM, HZ[0][2], HZ[0][3]], HZ[0].slice(4, 7), { dtInit: 0.1, maxDt: 20 });
    const tt = eph.toT(HZ[4][0]);
    while (ias.t < tt - 1e-10) ias.step(tt, 1e9);
    const sens = Math.hypot(ias.x[0] - r9.pos[3][0], ias.x[1] - r9.pos[3][1], ias.x[2] - r9.pos[3][2]);
    console.log(`INFO  flyby amplification: 1 km start offset in 2026 -> ${(sens * AU_KM).toFixed(0)} km difference in 2030 (explains the ${(r9.errs[3] * AU_KM).toFixed(0)} km model-limited gap from Horizons)`);
  }

  // (b) full pipeline from elements
  const t0 = Date.now();
  const res = simulate(APOPHIS, 2461322.5, { ephemeris: eph });
  const ms = Date.now() - t0;
  const dMin = Math.abs(res.minEarthDistAu - 0.000254);
  check('Apophis 2029 min distance within 1e-4 AU of 0.000254 AU (spec)', dMin < 1e-4, `${res.minEarthDistAu.toFixed(9)} AU = ${(res.minEarthDistAu * AU_KM).toFixed(0)} km`);
  check(
    'Apophis 2029 min distance vs Horizons (0.000254091 AU) within 1e-7 AU',
    Math.abs(res.minEarthDistAu - HZ_MIN_DIST_AU) < 1e-7,
    `diff ${fmt(res.minEarthDistAu - HZ_MIN_DIST_AU)} AU (${((res.minEarthDistAu - HZ_MIN_DIST_AU) * AU_KM).toFixed(2)} km)`,
  );
  const dtMin = (res.minEarthJd - HZ_MIN_JD) * 1440;
  check('Apophis 2029 closest-approach time within 10 min of Horizons (2029-04-13 21:46 UTC)', Math.abs(dtMin) < 10, `${dtMin.toFixed(2)} min`);
  check('Apophis run: window_end or escaped, not impact', res.status === 'escaped' || res.status === 'window_end', `status ${res.status}, ${res.nSamples} samples, ${res.stats.steps} steps, ${ms} ms`);
  const ca = res.events.find((e: any) => e.type === 'close_approach' && e.body === 'earth');
  check('Apophis close_approach event recorded with ~7.4 km/s relative speed', !!ca && Math.abs(ca.relSpeedKmS - 7.42) < 0.1, ca ? `${ca.relSpeedKmS.toFixed(3)} km/s` : 'none');

  // sample cadence: minute spacing inside the Hill sphere
  const n = res.nSamples;
  let minuteOk = 0;
  let inHill = 0;
  for (let i = 1; i < n; i++) {
    const o = i * SAMPLE_STRIDE;
    if (res.samples[o + S.DE] < 0.009 && res.samples[o - SAMPLE_STRIDE + S.DE] < 0.009) {
      inHill++;
      const dtMinutes = (res.samples[o] - res.samples[o - SAMPLE_STRIDE]) * 1440;
      if (Math.abs(dtMinutes - 1) < 0.01) minuteOk++;
    }
  }
  check('1-minute sampling inside 0.01 AU, daily outside', inHill > 1000 && minuteOk / inHill > 0.99, `${minuteOk}/${inHill} minute-spaced`);
  const dayGap = (res.samples[SAMPLE_STRIDE] - res.samples[0]);
  check('daily sampling far from Earth', Math.abs(dayGap - 1) < 1e-9, `first gap ${dayGap}`);

  // (c) escape classification on Apophis: escaped after the flyby, consistent with a brute-force
  //     run with escape disabled (min Earth distance after Hill exit must stay > 0.05 AU)
  const hillExit = res.events.find((e: any) => e.type === 'hill_exit');
  const full = simulate(APOPHIS, 2461322.5, { ephemeris: eph, escape: { enabled: false } });
  const bf = hillExit ? bruteReturns(full, hillExit.jd, S, SAMPLE_STRIDE) : true;
  check('Apophis: escape verdict agrees with brute force to 2036', (res.status === 'escaped') === !bf, `status ${res.status}, brute-force returns=${bf}`);

  // (d) backward direction from 2030 reproduces the 2029 approach
  const back = simulate(APOPHIS, 2462600.5, { ephemeris: eph, direction: -1, jdEnd: 2461322.5, escape: { enabled: false } });
  const backCa = back.events.find((e: any) => e.type === 'close_approach' && e.body === 'earth');
  check('backward integration finds the same 2029 approach (<1e-6 AU)', !!backCa && Math.abs(backCa.distAu - res.minEarthDistAu) < 1e-6, backCa ? `${backCa.distAu.toFixed(9)} AU` : 'none');

  // (e) streaming == non-streaming
  let nChunks = 0;
  let total = 0;
  const gen = simulateStream(APOPHIS, 2461322.5, { ephemeris: eph, chunkSamples: 1000 });
  let sumr: any;
  for (;;) {
    const r = gen.next();
    if (r.done) {
      sumr = r.value;
      break;
    }
    nChunks++;
    total += r.value.samples.length / SAMPLE_STRIDE;
  }
  check('streaming yields several chunks and the same sample count', nChunks > 3 && total === res.nSamples && sumr.stopJd === res.stopJd, `${nChunks} chunks, ${total} samples`);

  // (f) Hermite sample interpolation: daily run vs half-day run (nodes at midpoints)
  const daily = simulate(APOPHIS, 2461322.5, { ephemeris: eph, escape: { enabled: false }, jdEnd: 2461400.5, sampling: [{ maxDistAu: Infinity, stepDays: 1 }] });
  const half = simulate(APOPHIS, 2461322.5, { ephemeris: eph, escape: { enabled: false }, jdEnd: 2461400.5, sampling: [{ maxDistAu: Infinity, stepDays: 0.5 }] });
  const out = new Float64Array(6);
  let worst = 0;
  for (let i = 1; i < half.nSamples; i += 2) {
    const o = i * SAMPLE_STRIDE;
    sampleStateAt(daily.samples, daily.nSamples, half.samples[o], out);
    worst = Math.max(worst, Math.hypot(out[0] - half.samples[o + 1], out[1] - half.samples[o + 2], out[2] - half.samples[o + 3]));
  }
  check('sampleStateAt (cubic Hermite on daily samples) error < 1e-6 AU', worst < 1e-6, `worst ${fmt(worst)} AU`);
}

// ------------------------------------------------------------------------------------------
// 4. Synthetic cases: Earth impact, Moon impact, hyperbolic escape, out-of-coverage handling
// ------------------------------------------------------------------------------------------
{
  const jd = 2462100.5;
  const t = eph.toT(jd);
  const E = new Float64Array(6);
  const Mo = new Float64Array(6);
  eph.state(IDX_EARTH, t, E);
  eph.state(IDX_MOON, t, Mo);
  const k = 1 / AUD_TO_KMS;

  // Earth impactor: 0.008 AU out, aimed at Earth's centre, 8 km/s relative
  {
    const u = [0.6, 0.8, 0];
    const st = {
      jd,
      x: [E[0] + 0.008 * u[0], E[1] + 0.008 * u[1], E[2] + 0.008 * u[2]],
      v: [E[3] - 8 * k * u[0], E[4] - 8 * k * u[1], E[5] - 8 * k * u[2]],
    };
    const r = simulate(st, jd, { ephemeris: eph });
    const imp = r.events.find((e: any) => e.type === 'impact');
    check('synthetic impactor -> impact_earth', r.status === 'impact_earth' && !!imp && imp.body === 'earth', `status ${r.status}`);
    if (imp) {
      const last = r.samples.subarray((r.nSamples - 1) * SAMPLE_STRIDE);
      const dEnd = last[S.DE] * AU_KM;
      check('impact: final sample at Earth radius (6371 km) and speed > 11.2 km/s', Math.abs(dEnd - 6371) < 1 && imp.relSpeedKmS > 11.2, `d=${dEnd.toFixed(3)} km, v=${imp.relSpeedKmS.toFixed(2)} km/s, t=${((imp.jd - jd) * 24).toFixed(2)} h after start`);
      check('impact: stopJd equals impact time', Math.abs(r.stopJd - imp.jd) < 1e-6, `${fmt(r.stopJd - imp.jd)} d`);
    }
  }
  // Earth grazer just missing: impact parameter 1.5 R -> no impact
  {
    // aim at point offset by 20000 km perpendicular
    const u = [0.6, 0.8, 0];
    const p = [-0.8, 0.6, 0];
    const off = 20000 / AU_KM;
    const st = {
      jd,
      x: [E[0] + 0.008 * u[0] + off * p[0], E[1] + 0.008 * u[1] + off * p[1], E[2] + 0.008 * u[2]],
      v: [E[3] - 8 * k * u[0], E[4] - 8 * k * u[1], E[5] - 8 * k * u[2]],
    };
    const r = simulate(st, jd, { ephemeris: eph, jdEnd: jd + 10 });
    const ca = r.events.find((e: any) => e.type === 'close_approach' && e.body === 'earth');
    check('non-impacting flyby is not flagged as impact and records close approach', r.status !== 'impact_earth' && !!ca && ca.distKm > 6371, ca ? `perigee ${ca.distKm.toFixed(0)} km, status ${r.status}` : r.status);
  }
  // Moon impactor
  {
    const away = [Mo[0] - E[0], Mo[1] - E[1], Mo[2] - E[2]];
    const n = Math.hypot(away[0], away[1], away[2]);
    const u = away.map((c) => c / n);
    const d0 = 0.0003;
    const st = {
      jd,
      x: [Mo[0] + d0 * u[0], Mo[1] + d0 * u[1], Mo[2] + d0 * u[2]],
      v: [Mo[3] - 6 * k * u[0], Mo[4] - 6 * k * u[1], Mo[5] - 6 * k * u[2]],
    };
    const r = simulate(st, jd, { ephemeris: eph, jdEnd: jd + 2 });
    check('synthetic impactor aimed at the Moon -> impact_moon', r.status === 'impact_moon', `status ${r.status}`);
  }
  // Hyperbolic outbound (heliocentrically unbound) starting inside the Hill sphere
  {
    const u = [E[3], E[4], E[5]];
    const un = Math.hypot(u[0], u[1], u[2]);
    const dir = u.map((c) => c / un);
    const st = {
      jd,
      x: [E[0] + 0.004 * dir[0], E[1] + 0.004 * dir[1], E[2] + 0.004 * dir[2]],
      v: [E[3] + 25 * k * dir[0], E[4] + 25 * k * dir[1], E[5] + 25 * k * dir[2]],
    };
    const r = simulate(st, jd, { ephemeris: eph });
    const chk = r.events.find((e: any) => e.type === 'escape_check');
    check('hyperbolic outbound -> escaped', r.status === 'escaped' && !!chk && chk.result === 'escaped', `status ${r.status}, stop ${(r.stopJd - jd).toFixed(2)} d after start, escape at +${r.escapeJd ? (r.escapeJd - jd).toFixed(3) : '?'} d`);
    check('escaped run stops postEscapeDays (10 d) after Hill exit', !!r.escapeJd && Math.abs(r.stopJd - r.escapeJd - 10) < 1e-6);
  }
  // Slow exit: energy barely positive -> classification must agree with brute force
  {
    const u = [E[3], E[4], E[5]];
    const un = Math.hypot(u[0], u[1], u[2]);
    const dir = u.map((c) => c / un);
    const st = {
      jd,
      x: [E[0] + 0.004 * dir[0], E[1] + 0.004 * dir[1], E[2] + 0.004 * dir[2]],
      v: [E[3] + 1.5 * k * dir[0], E[4] + 1.5 * k * dir[1], E[5] + 1.5 * k * dir[2]],
    };
    const r = simulate(st, jd, { ephemeris: eph });
    const chk = r.events.find((e: any) => e.type === 'escape_check');
    const bf = simulate(st, jd, { ephemeris: eph, escape: { enabled: false } });
    const ex = bf.events.find((e: any) => e.type === 'hill_exit');
    const returns = ex ? bruteReturns(bf, ex.jd, S, SAMPLE_STRIDE) : true;
    check('slow-exit case: escape verdict agrees with brute-force re-entry check', !!chk && (chk.result === 'escaped') === !returns, `verdict ${chk ? chk.result : 'none'}, brute-force returns=${returns}, status ${r.status}`);
  }
  // Out-of-coverage handling
  {
    let code = '';
    try {
      simulate(APOPHIS, eph.jdStart - 100, { ephemeris: eph });
    } catch (e) {
      code = e instanceof NBodyRangeError ? e.code : String(e);
    }
    check('start before coverage -> typed NBodyRangeError(OUT_OF_COVERAGE)', code === 'OUT_OF_COVERAGE', code);
    const clamped = simulate(APOPHIS, eph.jdEnd + 500, { ephemeris: eph, clampToCoverage: true });
    check('clampToCoverage clamps instead of throwing', clamped.startJd === eph.jdEnd, `start ${clamped.startJd}`);
    // element epoch in 1987: re-epoch Apophis two-body then simulate from inside coverage
    const old = propagateElementsTwoBody(APOPHIS, 2447000.5, manifest.gm.sun.au3d2);
    const r = simulate(old, 2461322.5, { ephemeris: eph, jdEnd: 2461400.5, escape: { enabled: false } });
    const dd = Math.hypot(r.startState.x[0] - HZ[0][1], r.startState.x[1] - HZ[0][2], r.startState.x[2] - HZ[0][3]);
    check('element epoch before coverage (1987): two-body pre-propagation runs', r.twoBodyDays > 10000 && dd < 0.05, `twoBodyDays ${r.twoBodyDays.toFixed(0)}, start-state offset from Horizons ${dd.toFixed(5)} AU (sun-only drift)`);
  }
}

// ------------------------------------------------------------------------------------------
// 4b. Extended coverage 2020-01-01 .. 2100-12-31: lookups in 2029 and 2095, the 2036/2037 seam, the last
//     day, and range errors at both ends.
// ------------------------------------------------------------------------------------------
{
  // Independent JPL Horizons vectors at epochs BETWEEN table nodes (barycentric, ecliptic J2000, AU and AU/day;
  // moonRel is geocentric), fetched with horizons_vectors() from scripts/fetch_ephemeris.py.
  const REF = [
    // 2029-04-13 18:00 TDB (Apophis flyby day)
    {
      tag: '2029',
      jd: 2462240.25,
      sun: [0.001145488363936348, -0.0008552883404964813, 3.115325761124493e-05, -1.259173462994582e-06, 4.988712646574992e-06, 3.90014925368425e-09],
      earth: [-0.9172783408226306, -0.403674141790749, 6.063776518834709e-05, 0.006631885622686458, -0.01582127297872252, 1.182973034883342e-06],
      jupiter: [-5.043010050510045, -2.064019942629867, 0.121464110596568, 0.002768433443493603, -0.006633869454057104, -3.437328084284276e-05],
      neptune: [29.54165838339578, 4.284041734368702, -0.7690408568364652, -0.0004710178978599656, 0.003125279307509938, -5.350562135038597e-05],
      moonRel: [0.002510143933696842, 0.00101428915788273, 0.0002292893308936725, -0.0002124665649081209, 0.0005192026951418629, -1.277132162631147e-05],
    },
    // 2036-12-31 18:00 TDB (end of the original table = seam of the extension)
    {
      tag: 'seam 2036/2037',
      jd: 2465058.75,
      sun: [0.000412151369806379, -0.007683456469433188, 1.082290537922298e-05, 8.42155783371175e-06, -7.117632911680498e-09, -2.151509520749108e-07],
      earth: [-0.1640486467486204, 0.9617831045983397, -6.823581991059252e-05, -0.01722629881462544, -0.002943305347483311, -3.026711092551175e-07],
      jupiter: [0.687453510858268, 5.043885159109904, -0.0363552403397251, -0.00756498258562046, 0.001369725427186808, 0.0001635206751923058],
      neptune: [26.93581398956325, 12.77503800461293, -0.8838464407090906, -0.001364951169115302, 0.002855192984916657, -2.73418432920911e-05],
      moonRel: [0.0006535348355406739, 0.002390836753233106, -0.0001813591441759575, -0.0005972366098008301, 0.000130099481914453, 3.593094552579227e-05],
    },
    // 2095-07-04 06:00 TDB
    {
      tag: '2095',
      jd: 2486427.75,
      sun: [-5.973266312792907e-05, -0.006269490691810774, 1.391180764774364e-05, 7.116027349770342e-06, -3.103576834370313e-06, -1.846033399835872e-07],
      earth: [0.1989613387791475, -1.003292223981129, 0.0002232316837303462, 0.01660529725688274, 0.00330242431054286, -1.625888090089229e-06],
      jupiter: [2.782970674633181, 4.160413144222001, -0.07965382290716824, -0.006362492878198388, 0.004545143386101936, 0.0001230157721072935],
      neptune: [-27.2028606931433, 13.01808441836146, 0.3588694535622779, -0.001372523110442353, -0.002811278422014665, 8.952922988275442e-05],
      moonRel: [-0.001697329607085044, 0.001987028638656229, 0.0001902431450428128, -0.0004573617093394483, -0.0003564476273148808, 3.395310148667666e-05],
    },
    // 2100-12-30 18:00 TDB (last day of coverage)
    {
      tag: 'last day 2100-12-30',
      jd: 2488433.25,
      sun: [0.007243530431021522, 0.0038253730202079, -0.0002437805917179797, -4.865597165068992e-06, 6.823090584850088e-06, 7.42737744975605e-08],
      earth: [-0.1242004603510599, 0.9784509566663336, -0.0004636321496618933, -0.01734153868826349, -0.002358017852440775, 6.07874291759428e-07],
      jupiter: [-4.321786055926784, -3.276175326622243, 0.1102958092558852, 0.004461439041814476, -0.005666210811756518, -7.577207693529152e-05],
      neptune: [-29.34784010861532, 7.139407593488341, 0.5293741991642005, -0.0007593465312969012, -0.003029933716324338, 7.990148902734419e-05],
      moonRel: [0.000179091731122767, -0.002415622342042606, -0.0001753931230232998, 0.0006252763364614877, 2.768489377531903e-05, 3.108728504158345e-05],
    },
  ];
  const IDX = { sun: IDX_SUN, earth: IDX_EARTH, jupiter: BODY_ORDER.indexOf('jupiter'), neptune: BODY_ORDER.indexOf('neptune') };
  const st = new Float64Array(6);
  const pos = new Float64Array(3 * N_BODIES);
  const km = (a: ArrayLike<number>, b: ArrayLike<number>) => Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]) * AU_KM;

  check(
    'coverage is 2020-01-01 .. 2100-12-31 (JD 2458849.5 .. 2488433.5)',
    eph.jdStart === 2458849.5 && eph.jdEnd === 2488433.5,
    `JD ${eph.jdStart} .. ${eph.jdEnd}, Moon grid ${(manifest.moon.stepDays * 24).toFixed(2)} h`,
  );
  for (const r of REF) {
    const t = eph.toT(r.jd);
    let worstP = 0;
    let worstV = 0;
    for (const name of ['sun', 'earth', 'jupiter', 'neptune'] as const) {
      eph.state(IDX[name], t, st);
      worstP = Math.max(worstP, km(st, r[name]));
      worstV = Math.max(worstV, Math.hypot(st[3] - r[name][3], st[4] - r[name][4], st[5] - r[name][5]) * AUD_TO_KMS * 1e6); // mm/s
    }
    eph.moonRelState(t, st);
    const moonKm = km(st, r.moonRel);
    const moonV = Math.hypot(st[3] - r.moonRel[3], st[4] - r.moonRel[4], st[5] - r.moonRel[5]) * AUD_TO_KMS * 1e6; // mm/s
    // positions() (used by the force model) must agree with state()
    eph.positions(t, pos);
    eph.state(IDX_MOON, t, st);
    const consist = km([pos[3 * IDX_MOON], pos[3 * IDX_MOON + 1], pos[3 * IDX_MOON + 2]], st);
    check(
      `${r.tag}: ephemeris lookup matches Horizons (Sun/Earth/Jupiter/Neptune < 0.2 km and 20 mm/s, Moon < 1 km and 0.3 m/s, positions()==state())`,
      worstP < 0.2 && worstV < 20 && moonKm < 1 && moonV < 300 && consist < 1e-6,
      `planets worst ${(worstP * 1000).toFixed(1)} m / ${worstV.toFixed(2)} mm/s, Moon(geocentric) ${(moonKm * 1000).toFixed(1)} m / ${moonV.toFixed(1)} mm/s`,
    );
  }
  // positions() at the very last node and first node must work (inclusive coverage)
  {
    let ok = true;
    try {
      eph.positions(eph.tEnd, pos);
      eph.positions(eph.tStart, pos);
    } catch {
      ok = false;
    }
    check('first and last node of the table are both readable', ok);
  }
  // Ephemeris-level range errors
  {
    const thrown = (fn: () => void) => {
      try {
        fn();
      } catch (e) {
        return e instanceof RangeError;
      }
      return false;
    };
    check(
      'ephemeris lookups outside coverage throw RangeError (both ends, state + positions)',
      thrown(() => eph.positions(eph.tEnd + 0.01, pos)) &&
        thrown(() => eph.state(IDX_EARTH, eph.tEnd + 1, st)) &&
        thrown(() => eph.positions(eph.tStart - 0.01, pos)) &&
        thrown(() => eph.state(IDX_MOON, eph.tStart - 1, st)),
    );
  }
  // A short N-body run in 2095 and one that reaches the very end of the table
  {
    const jd95 = 2486427.75;
    eph.state(IDX_EARTH, eph.toT(jd95), st);
    const seed = {
      jd: jd95,
      x: [st[0] + 0.03, st[1], st[2]],
      v: [st[3], st[4] + 0.0004, st[5]],
    };
    const r = simulate(seed, jd95, { ephemeris: eph, jdEnd: jd95 + 120, escape: { enabled: false } });
    check(
      '2095: N-body run (120 d) completes, reaches window end, finite samples',
      r.status === 'window_end' && Math.abs(r.stopJd - (jd95 + 120)) < 1e-6 && r.nSamples > 10 && r.samples.every(Number.isFinite),
      `${r.nSamples} samples, ${r.stats.steps} steps, min Earth dist ${r.minEarthDistAu.toFixed(5)} AU`,
    );
    const last = (r.nSamples - 1) * SAMPLE_STRIDE;
    const rSun = Math.hypot(r.samples[last + S.HX], r.samples[last + S.HX + 1], r.samples[last + S.HX + 2]);
    check('2095: particle seeded next to Earth is still on a ~1 AU heliocentric orbit after 120 d', rSun > 0.9 && rSun < 1.1, `heliocentric distance ${rSun.toFixed(4)} AU`);

    const jdTail = eph.jdEnd - 25;
    eph.state(IDX_EARTH, eph.toT(jdTail), st);
    const seed2 = { jd: jdTail, x: [st[0] + 0.03, st[1], st[2]], v: [st[3], st[4] + 0.0004, st[5]] };
    const r2 = simulate(seed2, jdTail, { ephemeris: eph, escape: { enabled: false } });
    check('run started 25 d before the end of the table stops exactly at coverage end (no range error)', Math.abs(r2.stopJd - eph.jdEnd) < 1e-6 && r2.samples.every(Number.isFinite), `stop JD ${r2.stopJd}, ${r2.nSamples} samples`);
  }
  // Typed NBodyRangeError at both ends (state + element inputs)
  {
    const codeOf = (fn: () => void) => {
      try {
        fn();
      } catch (e) {
        return e instanceof NBodyRangeError ? e.code : `other:${String(e)}`;
      }
      return 'no-throw';
    };
    const seed = { jd: 2486427.75, x: [0.2, -1.0, 0], v: [0.0166, 0.0033, 0] };
    check(
      'start after coverage end (2101-01-01) -> NBodyRangeError(OUT_OF_COVERAGE)',
      codeOf(() => simulate(seed, 2488434.5, { ephemeris: eph })) === 'OUT_OF_COVERAGE' &&
        codeOf(() => simulate(APOPHIS, eph.jdEnd + 1, { ephemeris: eph })) === 'OUT_OF_COVERAGE',
    );
    check(
      'start before coverage start (2019-12-31) -> NBodyRangeError(OUT_OF_COVERAGE)',
      codeOf(() => simulate(seed, eph.jdStart - 1, { ephemeris: eph })) === 'OUT_OF_COVERAGE' &&
        codeOf(() => simulate({ ...seed, jd: eph.jdEnd + 5 }, eph.jdEnd - 5, { ephemeris: eph })) === 'OUT_OF_COVERAGE',
    );
    const clamped = simulate(APOPHIS, eph.jdEnd + 3000, { ephemeris: eph, clampToCoverage: true });
    check('clampToCoverage clamps a start far beyond 2100 to the last day', clamped.startJd === eph.jdEnd, `start ${clamped.startJd}`);
  }
}

// ------------------------------------------------------------------------------------------
// 5. Rounded-element behaviour (asteroids.json precision) - informational
// ------------------------------------------------------------------------------------------
{
  const rounded = { ...APOPHIS, a: 0.922, e: 0.1911, i: 3.34, om: 203.89, w: 126.68, ma: 175.33, epoch: 2461200.5 };
  const r = simulate(rounded, 2461322.5, { ephemeris: eph, jdEnd: 2462400.5, escape: { enabled: false } });
  console.log(`INFO  Apophis from 3-4 digit rounded elements: min Earth distance ${r.minEarthDistAu.toFixed(6)} AU on JD ${r.minEarthJd.toFixed(3)} (full precision: 0.000254 AU on 2462240.407)`);
}

console.log(`\n${pass} passed, ${fail} failed`);
if (fail) {
  console.log('FAILED: ' + failures.join('; '));
  process.exit(1);
}
void StepPoly;
void IDX_SUN;
