// Trajectory simulation for ONE asteroid (massless test particle) with event detection.
//
//   simulate(asteroid, jdStart, opts)        -> SimResult   (all samples in memory)
//   simulateStream(asteroid, jdStart, opts)  -> generator of SimChunk, returns SimSummary
//
// All public times are JD (TDB).  Positions: AU, ecliptic J2000.  Velocities: AU/day.
// Internally time is t = jd - ephemeris.jd0 (full double resolution, see ephemeris.ts).

import {
  AUD_TO_KMS,
  AU_KM,
  EARTH_RADIUS_AU,
  HILL_RADIUS_AU,
  MOON_RADIUS_AU,
  RETURN_RADIUS_AU,
} from './constants.ts';
import { makeForce } from './dynamics.ts';
import { elementsToState, isState, propagateElementsTwoBody } from './elements.ts';
import type { AsteroidElements, AsteroidState } from './elements.ts';
import { Ephemeris, IDX_EARTH, IDX_MOON, IDX_SUN } from './ephemeris.ts';
import { IAS15, StepPoly } from './ias15.ts';

// ---------------------------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------------------------

/** Floats per trajectory sample: [jd, gx,gy,gz, hx,hy,hz, gvx,gvy,gvz, dEarthAu, dMoonAu]. */
export const SAMPLE_STRIDE = 12;
export const S = { JD: 0, GX: 1, HX: 4, GVX: 7, DE: 10, DM: 11 } as const;

export type TerminalStatus = 'impact_earth' | 'impact_moon' | 'escaped' | 'window_end';

export interface SamplingTier {
  /** applies while the geocentric distance is <= maxDistAu */
  maxDistAu: number;
  /** sample spacing in days */
  stepDays: number;
}

export const DEFAULT_SAMPLING: SamplingTier[] = [
  { maxDistAu: 0.01, stepDays: 1 / 1440 }, // inside Earth's Hill sphere: every minute
  { maxDistAu: 0.05, stepDays: 1 / 24 }, // hourly
  { maxDistAu: Infinity, stepDays: 1 }, // daily
];

export interface EscapeOptions {
  enabled?: boolean; // default true
  /**
   * Only declare escape after the object has been inside the Hill sphere during this simulation
   * (or starts inside it).  Default true: otherwise virtually every NEO "escapes" at t0 because
   * it is outside the Hill sphere, recedes and never comes back within 0.05 AU.
   */
  requireEncounter?: boolean;
  hillAu?: number; // 0.01
  returnAu?: number; // 0.05
  /** keep sampling this many days after leaving the Hill sphere before stopping. Default 10. */
  postEscapeDays?: number;
}

export interface SimOptions {
  ephemeris?: Ephemeris;
  /** +1 forward (default) or -1 backward in time. */
  direction?: 1 | -1;
  /** Stop time (JD).  Default: end (or start, if backward) of ephemeris coverage. */
  jdEnd?: number;
  /** IAS15 tolerance (default 1e-9 ~ machine precision for smooth problems). */
  epsilon?: number;
  /** Sun 1PN term (default true). */
  relativity?: boolean;
  sampling?: SamplingTier[];
  /** Record Earth close approaches below this distance (AU). Default 0.05. */
  closeApproachAu?: number;
  /** Record Moon close approaches below this distance (AU). Default 0.01. */
  closeApproachMoonAu?: number;
  escape?: EscapeOptions;
  /** Samples per streamed chunk. Default 4096. */
  chunkSamples?: number;
  /** Clamp jdStart (and jdEnd) into the ephemeris coverage instead of throwing NBodyRangeError. */
  clampToCoverage?: boolean;
  /** Hard cap on stored samples (truncates with `truncated: true`). Default 600000. */
  maxSamples?: number;
  /** Hard cap on integrator steps. Default 3e6. */
  maxSteps?: number;
}

export type SimEvent =
  | { type: 'close_approach'; body: 'earth' | 'moon'; jd: number; distAu: number; distKm: number; relSpeedKmS: number }
  | { type: 'impact'; body: 'earth' | 'moon'; jd: number; relSpeedKmS: number; relPosAu: [number, number, number] }
  | { type: 'hill_entry' | 'hill_exit'; jd: number; distAu: number }
  | {
      type: 'escape_check';
      jd: number;
      result: 'escaped' | 'continues';
      /** when 'continues': first return inside returnAu found by the look-ahead */
      returnJd?: number;
      returnDistAu?: number;
    };

export interface SimChunk {
  /** n * SAMPLE_STRIDE doubles; transferable. */
  samples: Float64Array;
  /** events found since the previous chunk */
  events: SimEvent[];
  /** latest time covered (JD) */
  jdReached: number;
  /** fraction of the requested span done, 0..1 */
  progress: number;
}

export interface SimSummary {
  status: TerminalStatus;
  /** time at which the trajectory ends (impact time, end of post-escape tail, window end) */
  stopJd: number;
  startJd: number;
  direction: 1 | -1;
  events: SimEvent[];
  /** closest Earth approach seen (centre distance, AU) and when */
  minEarthDistAu: number;
  minEarthJd: number;
  nSamples: number;
  /** true if maxSamples/maxSteps cut the run short */
  truncated: boolean;
  /** JD of the escape (Hill-sphere exit) when status is 'escaped' */
  escapeJd?: number;
  /** barycentric state at jdStart (after propagating from the element epoch) */
  /** days of Sun-only two-body propagation used to bring out-of-coverage element epochs into the ephemeris range */
  twoBodyDays: number;
  startState: { jd: number; x: [number, number, number]; v: [number, number, number] };
  stats: { steps: number; rejected: number; evals: number };
}

export interface SimResult extends SimSummary {
  /** nSamples * SAMPLE_STRIDE doubles */
  samples: Float64Array;
}

/** Thrown when a requested time is outside the ephemeris coverage (typed so UIs can handle it). */
export class NBodyRangeError extends RangeError {
  readonly code = 'OUT_OF_COVERAGE' as const;
  readonly coverageJd: [number, number];
  readonly requestedJd: number;
  constructor(message: string, requestedJd: number, coverageJd: [number, number]) {
    super(message);
    this.name = 'NBodyRangeError';
    this.requestedJd = requestedJd;
    this.coverageJd = coverageJd;
  }
}

let defaultEphemeris: Ephemeris | null = null;
export function setDefaultEphemeris(e: Ephemeris): void {
  defaultEphemeris = e;
}

// ---------------------------------------------------------------------------------------------
// Internals
// ---------------------------------------------------------------------------------------------

const NSCAN = 16; // sub-intervals per step when scanning for events

interface Scratch {
  X: Float64Array;
  V: Float64Array;
  E: Float64Array;
  Mo: Float64Array;
  Sun: Float64Array;
  tmp: Float64Array;
}
function makeScratch(): Scratch {
  return {
    X: new Float64Array(3),
    V: new Float64Array(3),
    E: new Float64Array(6),
    Mo: new Float64Array(6),
    Sun: new Float64Array(6),
    tmp: new Float64Array(6),
  };
}

function evalBodies(eph: Ephemeris, t: number, s: Scratch): void {
  eph.state(IDX_EARTH, t, s.E);
  eph.moonRelState(t, s.tmp);
  for (let c = 0; c < 6; c++) s.Mo[c] = s.E[c] + s.tmp[c];
}

/** distance to Earth/Moon centre and tau-derivative-sign quantity H = dt * (rel . relv). */
interface Geo {
  dE: number;
  hE: number;
  dM: number;
  hM: number;
}

function geoAt(eph: Ephemeris, poly: StepPoly, tau: number, s: Scratch, g: Geo): void {
  poly.eval(tau, s.X, s.V);
  evalBodies(eph, poly.t0 + poly.dt * tau, s);
  let rx = s.X[0] - s.E[0], ry = s.X[1] - s.E[1], rz = s.X[2] - s.E[2];
  let vx = s.V[0] - s.E[3], vy = s.V[1] - s.E[4], vz = s.V[2] - s.E[5];
  g.dE = Math.hypot(rx, ry, rz);
  g.hE = poly.dt * (rx * vx + ry * vy + rz * vz);
  rx = s.X[0] - s.Mo[0];
  ry = s.X[1] - s.Mo[1];
  rz = s.X[2] - s.Mo[2];
  vx = s.V[0] - s.Mo[3];
  vy = s.V[1] - s.Mo[4];
  vz = s.V[2] - s.Mo[5];
  g.dM = Math.hypot(rx, ry, rz);
  g.hM = poly.dt * (rx * vx + ry * vy + rz * vz);
}

/** step-size cap so a step cannot tunnel through Earth/Moon (fraction of d/v). */
function stepCap(eph: Ephemeris, t: number, x: Float64Array, v: Float64Array, s: Scratch): number {
  evalBodies(eph, t, s);
  let cap = Infinity;
  for (let b = 0; b < 2; b++) {
    const B = b === 0 ? s.E : s.Mo;
    const rx = x[0] - B[0], ry = x[1] - B[1], rz = x[2] - B[2];
    const vx = v[0] - B[3], vy = v[1] - B[4], vz = v[2] - B[5];
    const d = Math.hypot(rx, ry, rz);
    const sp = Math.hypot(vx, vy, vz);
    // 0.2*d/v, but never below a hair above the surface crossing time
    cap = Math.min(cap, (0.2 * d) / Math.max(sp, 1e-12));
  }
  return cap;
}

function centralMuEarthMoon(eph: Ephemeris): number {
  return eph.gm[IDX_EARTH] + eph.gm[IDX_MOON];
}

/**
 * Look-ahead from a state: after the object has left the sphere of radius returnAu around Earth,
 * does it come back inside it before the end of coverage (in direction dir)?  Same force model
 * and integrator as the main run.  `returns` is also true if it never gets beyond returnAu.
 */
export function lookaheadReturn(
  eph: Ephemeris,
  t: number,
  x: ArrayLike<number>,
  v: ArrayLike<number>,
  dir: 1 | -1,
  opts: { epsilon?: number; relativity?: boolean; returnAu?: number; tEnd?: number } = {},
): { returns: boolean; minDistAu: number; minJd: number } {
  const returnAu = opts.returnAu ?? RETURN_RADIUS_AU;
  const tEnd = opts.tEnd ?? (dir > 0 ? eph.tEnd : eph.tStart);
  const force = makeForce(eph, { relativity: opts.relativity ?? true });
  const ias = new IAS15(force, t, x, v, { epsilon: opts.epsilon ?? 1e-9, dtInit: 0.05 * dir, maxDt: 20 });
  const s = makeScratch();
  const g: Geo = { dE: 0, hE: 0, dM: 0, hM: 0 };
  const poly = new StepPoly();
  let outside = false;
  let minAll = Infinity; // min distance over the whole look-ahead (reported if it never leaves)
  let minAllT = t;
  let steps = 0;
  const NS = 12;
  while ((tEnd - ias.t) * dir > 1e-9 && steps++ < 2_000_000) {
    ias.step(tEnd, stepCap(eph, ias.t, ias.x, ias.v, s), poly);
    geoAt(eph, poly, 0, s, g);
    const d0 = g.dE;
    let hPrev = g.hE;
    geoAt(eph, poly, 1, s, g);
    const d1 = g.dE;
    if (d1 < minAll) {
      minAll = d1;
      minAllT = poly.t1;
    }
    if (!outside) {
      if (d1 >= returnAu) outside = true;
      else continue;
    }
    const L = 1.5 * Math.abs(poly.dt) * (Math.hypot(poly.v0[0], poly.v0[1], poly.v0[2]) + 0.03);
    if ((d0 + d1 - L) / 2 > returnAu * 1.5) continue;
    let tauPrev = 0;
    geoAt(eph, poly, 0, s, g);
    hPrev = g.hE;
    for (let m = 1; m <= NS; m++) {
      const tau = m / NS;
      geoAt(eph, poly, tau, s, g);
      if (g.dE < returnAu) return { returns: true, minDistAu: g.dE, minJd: eph.toJd(poly.t0 + poly.dt * tau) };
      if (hPrev < 0 && g.hE >= 0) {
        let lo = tauPrev;
        let hi = tau;
        for (let it = 0; it < 50; it++) {
          const mid = 0.5 * (lo + hi);
          geoAt(eph, poly, mid, s, g);
          if (g.hE < 0) lo = mid;
          else hi = mid;
        }
        const tm = 0.5 * (lo + hi);
        geoAt(eph, poly, tm, s, g);
        if (g.dE < returnAu) return { returns: true, minDistAu: g.dE, minJd: eph.toJd(poly.t0 + poly.dt * tm) };
        geoAt(eph, poly, tau, s, g);
      }
      hPrev = g.hE;
      tauPrev = tau;
    }
  }
  return { returns: !outside, minDistAu: minAll, minJd: eph.toJd(minAllT) };
}

/**
 * Barycentric initial state (internal time t).  If the element epoch lies outside the ephemeris
 * coverage the elements are first advanced by Sun-only two-body motion to 1 day inside the
 * nearest edge of the coverage (planetary perturbations over that stretch are ignored).
 */
function initialState(
  eph: Ephemeris,
  a: AsteroidElements | AsteroidState,
): { t: number; x: Float64Array; v: Float64Array; twoBodyDays: number } {
  if (isState(a)) {
    const t = eph.toT(a.jd);
    if (!eph.covers(t)) throw new NBodyRangeError(`state epoch JD ${a.jd} outside ephemeris coverage`, a.jd, [eph.jdStart, eph.jdEnd]);
    return { t, x: Float64Array.from(a.x), v: Float64Array.from(a.v), twoBodyDays: 0 };
  }
  let el = a;
  let t = eph.toT(a.epoch);
  let twoBodyDays = 0;
  if (!eph.covers(t)) {
    const tb = t < eph.tStart ? eph.tStart + 1 : eph.tEnd - 1;
    el = propagateElementsTwoBody(a, eph.toJd(tb), eph.gm[IDX_SUN]);
    twoBodyDays = Math.abs(tb - t);
    t = tb;
  }
  const { x, v } = elementsToState(el, eph.gm[IDX_SUN]);
  const sun = new Float64Array(6);
  eph.sunState(t, sun);
  return {
    t,
    x: Float64Array.from([x[0] + sun[0], x[1] + sun[1], x[2] + sun[2]]),
    v: Float64Array.from([v[0] + sun[3], v[1] + sun[4], v[2] + sun[5]]),
    twoBodyDays,
  };
}

interface PendingEvent {
  tau: number;
  kind: 'ca' | 'impact' | 'hill';
  body: 0 | 1; // 0 earth, 1 moon
  rising?: boolean; // hill: true = exit
  dist?: number;
}

// ---------------------------------------------------------------------------------------------
// Streaming simulation
// ---------------------------------------------------------------------------------------------

export function* simulateStream(
  asteroid: AsteroidElements | AsteroidState,
  jdStart: number,
  opts: SimOptions = {},
): Generator<SimChunk, SimSummary, undefined> {
  const eph = opts.ephemeris ?? defaultEphemeris;
  if (!eph) throw new Error('simulate: no ephemeris (pass opts.ephemeris or call setDefaultEphemeris)');
  const dir: 1 | -1 = opts.direction ?? 1;
  const epsilon = opts.epsilon ?? 1e-9;
  const tiers = (opts.sampling ?? DEFAULT_SAMPLING).slice().sort((a, b) => a.maxDistAu - b.maxDistAu);
  const caEarth = opts.closeApproachAu ?? 0.05;
  const caMoon = opts.closeApproachMoonAu ?? 0.01;
  const esc = opts.escape ?? {};
  const escEnabled = esc.enabled ?? true;
  const escRequireEncounter = esc.requireEncounter ?? true;
  const hillAu = esc.hillAu ?? HILL_RADIUS_AU;
  const returnAu = esc.returnAu ?? RETURN_RADIUS_AU;
  const postEscapeDays = esc.postEscapeDays ?? 10;
  const chunkSamples = opts.chunkSamples ?? 4096;
  const maxSamples = opts.maxSamples ?? 600_000;
  const maxSteps = opts.maxSteps ?? 3_000_000;
  const muEM = centralMuEarthMoon(eph);

  const cov: [number, number] = [eph.jdStart, eph.jdEnd];
  if (opts.clampToCoverage) {
    jdStart = Math.min(Math.max(jdStart, cov[0]), cov[1]);
    if (opts.jdEnd !== undefined) opts = { ...opts, jdEnd: Math.min(Math.max(opts.jdEnd, cov[0]), cov[1]) };
  }
  const tStart = eph.toT(jdStart);
  if (!Number.isFinite(jdStart) || !eph.covers(tStart)) {
    throw new NBodyRangeError(`simulate: start JD ${jdStart} outside ephemeris coverage ${cov[0]}..${cov[1]}`, jdStart, cov);
  }
  const covEnd = dir > 0 ? eph.tEnd : eph.tStart;
  let tEndUser = opts.jdEnd !== undefined ? eph.toT(opts.jdEnd) : covEnd;
  if ((tEndUser - tStart) * dir < 0) throw new RangeError('simulate: jdEnd is behind jdStart for the chosen direction');
  if ((tEndUser - covEnd) * dir > 0) tEndUser = covEnd;
  const span = Math.abs(tEndUser - tStart) || 1;

  const force = makeForce(eph, { relativity: opts.relativity ?? true });
  const s = makeScratch();
  const g: Geo = { dE: 0, hE: 0, dM: 0, hM: 0 };

  // ---- initial state at the epoch, propagated to jdStart (no sampling / events) ----
  const init = initialState(eph, asteroid);
  let ias: IAS15;
  {
    const toStart = tStart - init.t;
    const dirProp = toStart >= 0 ? 1 : -1;
    ias = new IAS15(force, init.t, init.x, init.v, { epsilon, dtInit: 0.01 * dirProp, maxDt: 20 });
    if (Math.abs(toStart) > 1e-10) {
      const poly = new StepPoly();
      let guard = 0;
      while ((tStart - ias.t) * dirProp > 1e-10 && guard++ < maxSteps) {
        ias.step(tStart, stepCap(eph, ias.t, ias.x, ias.v, s), poly);
      }
    }
    if (dirProp !== dir) ias.reset(ias.t, ias.x, ias.v, dir * Math.max(Math.abs(ias.dt), 1e-3));
    else if (Math.abs(toStart) <= 1e-10) ias.dt = dir * 0.01;
  }
  const startState = {
    jd: jdStart,
    x: [ias.x[0], ias.x[1], ias.x[2]] as [number, number, number],
    v: [ias.v[0], ias.v[1], ias.v[2]] as [number, number, number],
  };

  // ---- output buffering ----
  let buf = new Float64Array(Math.max(chunkSamples, 256) * SAMPLE_STRIDE);
  let nBuf = 0;
  let nTotal = 0;
  let truncated = false;
  const allEvents: SimEvent[] = [];
  let pendingEvents: SimEvent[] = [];
  const pushEvent = (e: SimEvent): void => {
    allEvents.push(e);
    pendingEvents.push(e);
  };

  const stepFor = (d: number): number => {
    for (const tr of tiers) if (d <= tr.maxDistAu) return tr.stepDays;
    return tiers[tiers.length - 1].stepDays;
  };

  const writeSample = (t: number, X: Float64Array, V: Float64Array): number => {
    evalBodies(eph, t, s);
    eph.sunState(t, s.Sun);
    if ((nBuf + 1) * SAMPLE_STRIDE > buf.length) {
      const nb = new Float64Array(buf.length * 2);
      nb.set(buf);
      buf = nb;
    }
    const o = nBuf * SAMPLE_STRIDE;
    const gx = X[0] - s.E[0], gy = X[1] - s.E[1], gz = X[2] - s.E[2];
    buf[o + S.JD] = eph.toJd(t);
    buf[o + S.GX] = gx;
    buf[o + S.GX + 1] = gy;
    buf[o + S.GX + 2] = gz;
    buf[o + S.HX] = X[0] - s.Sun[0];
    buf[o + S.HX + 1] = X[1] - s.Sun[1];
    buf[o + S.HX + 2] = X[2] - s.Sun[2];
    buf[o + S.GVX] = V[0] - s.E[3];
    buf[o + S.GVX + 1] = V[1] - s.E[4];
    buf[o + S.GVX + 2] = V[2] - s.E[5];
    const dE = Math.hypot(gx, gy, gz);
    buf[o + S.DE] = dE;
    buf[o + S.DM] = Math.hypot(X[0] - s.Mo[0], X[1] - s.Mo[1], X[2] - s.Mo[2]);
    nBuf++;
    nTotal++;
    return dE;
  };

  // ---- sampling state ----
  let lastSampleT = tStart;
  let nextSampleT = tStart;
  let haveFirstSample = false;

  const emitSamples = (poly: StepPoly, tCut: number): void => {
    // emit all samples with time <= tCut (in direction) that lie inside this step
    for (;;) {
      if (nTotal >= maxSamples) {
        truncated = true;
        return;
      }
      let tn = nextSampleT;
      if ((tn - tCut) * dir > 1e-12) return;
      let atStepStart = false;
      if ((tn - poly.t0) * dir <= 0) {
        tn = poly.t0; // candidate fell in an earlier step: emit at this step's start instead
        atStepStart = true;
      }
      const tau = (tn - poly.t0) / poly.dt;
      poly.eval(tau, s.X, s.V);
      evalBodies(eph, tn, s);
      const dEn = Math.hypot(s.X[0] - s.E[0], s.X[1] - s.E[1], s.X[2] - s.E[2]);
      if (haveFirstSample && !atStepStart) {
        const spacing = Math.abs(tn - lastSampleT);
        const fine = stepFor(dEn);
        if (fine < spacing * 0.999 && spacing > fine * 1.0001) {
          nextSampleT = lastSampleT + dir * Math.max(fine, spacing / 2);
          continue;
        }
      }
      if (haveFirstSample && Math.abs(tn - lastSampleT) < 1e-12) {
        nextSampleT = tn + dir * stepFor(dEn);
        continue;
      }
      writeSample(tn, s.X, s.V);
      haveFirstSample = true;
      lastSampleT = tn;
      nextSampleT = tn + dir * stepFor(dEn);
    }
  };

  const flushChunk = (tNow: number): SimChunk => {
    const chunk: SimChunk = {
      samples: buf.slice(0, nBuf * SAMPLE_STRIDE),
      events: pendingEvents,
      jdReached: eph.toJd(tNow),
      progress: Math.min(1, Math.abs(tNow - tStart) / span),
    };
    nBuf = 0;
    pendingEvents = [];
    return chunk;
  };

  // ---- geometric helpers on a step polynomial ----
  const radius = [EARTH_RADIUS_AU, MOON_RADIUS_AU];
  const distAt = (poly: StepPoly, tau: number, body: 0 | 1): number => {
    geoAt(eph, poly, tau, s, g);
    return body === 0 ? g.dE : g.dM;
  };
  const refineMin = (poly: StepPoly, a: number, b: number, body: 0 | 1): number => {
    let lo = a;
    let hi = b;
    for (let it = 0; it < 55; it++) {
      const mid = 0.5 * (lo + hi);
      geoAt(eph, poly, mid, s, g);
      if ((body === 0 ? g.hE : g.hM) < 0) lo = mid;
      else hi = mid;
    }
    return 0.5 * (lo + hi);
  };
  /** root of D(tau) - target in [a,b] where (D(a)-target) and (D(b)-target) differ in sign */
  const refineCross = (poly: StepPoly, a: number, b: number, body: 0 | 1, target: number): number => {
    let lo = a;
    let hi = b;
    const sLo = distAt(poly, lo, body) - target > 0;
    for (let it = 0; it < 60; it++) {
      const mid = 0.5 * (lo + hi);
      const sm = distAt(poly, mid, body) - target > 0;
      if (sm === sLo) lo = mid;
      else hi = mid;
    }
    return 0.5 * (lo + hi);
  };
  const relSpeedAt = (poly: StepPoly, tau: number, body: 0 | 1): number => {
    poly.eval(tau, s.X, s.V);
    evalBodies(eph, poly.t0 + poly.dt * tau, s);
    const B = body === 0 ? s.E : s.Mo;
    return Math.hypot(s.V[0] - B[3], s.V[1] - B[4], s.V[2] - B[5]) * AUD_TO_KMS;
  };

  // ---- main loop ----
  let status: TerminalStatus = 'window_end';
  let stopT = tEndUser; // may shrink after an escape is declared
  let escapeT: number | undefined;
  let encountered = false;
  let minEarthD = Infinity;
  let minEarthT = tStart;
  const poly = new StepPoly();
  const dEs = new Float64Array(NSCAN + 1);
  const hEs = new Float64Array(NSCAN + 1);
  const dMs = new Float64Array(NSCAN + 1);
  const hMs = new Float64Array(NSCAN + 1);
  let lastHE = NaN;
  let lastHM = NaN;
  let steps = 0;
  let tNow = tStart;
  let endedInside = false;
  const hillNeeded = escEnabled;

  // state at start: inside Earth/Moon? inside Hill? escape check at start?
  {
    const p0 = new StepPoly();
    p0.t0 = tStart;
    p0.dt = dir * 1e-9;
    p0.x0.set(ias.x);
    p0.v0.set(ias.v);
    p0.a0.set(ias.a);
    geoAt(eph, p0, 0, s, g);
    minEarthD = g.dE;
    if (g.dE < hillAu) encountered = true;
    if (g.dE < EARTH_RADIUS_AU || g.dM < MOON_RADIUS_AU) {
      status = g.dE < EARTH_RADIUS_AU ? 'impact_earth' : 'impact_moon';
      const body: 0 | 1 = g.dE < EARTH_RADIUS_AU ? 0 : 1;
      pushEvent({
        type: 'impact',
        body: body === 0 ? 'earth' : 'moon',
        jd: jdStart,
        relSpeedKmS: relSpeedAt(p0, 0, body),
        relPosAu: [s.X[0] - (body === 0 ? s.E : s.Mo)[0], s.X[1] - (body === 0 ? s.E : s.Mo)[1], s.X[2] - (body === 0 ? s.E : s.Mo)[2]],
      });
      stopT = tStart;
    } else if (hillNeeded && !escRequireEncounter && g.dE > hillAu) {
      // not requiring an encounter: test escape right at the start
      const rx = s.X[0] - s.E[0], ry = s.X[1] - s.E[1], rz = s.X[2] - s.E[2];
      const vx = s.V[0] - s.E[3], vy = s.V[1] - s.E[4], vz = s.V[2] - s.E[5];
      const en = 0.5 * (vx * vx + vy * vy + vz * vz) - muEM / Math.hypot(rx, ry, rz);
      if (en > 0 && (rx * vx + ry * vy + rz * vz) * dir > 0) {
        const la = lookaheadReturn(eph, tStart, ias.x, ias.v, dir, { epsilon, relativity: opts.relativity, returnAu });
        pushEvent({
          type: 'escape_check',
          jd: jdStart,
          result: la.returns ? 'continues' : 'escaped',
          returnJd: la.returns ? la.minJd : undefined,
          returnDistAu: la.returns ? la.minDistAu : undefined,
        });
        if (!la.returns) {
          escapeT = tStart;
          stopT = tStart + dir * postEscapeDays;
          if ((stopT - tEndUser) * dir > 0) stopT = tEndUser;
        }
      }
    }
    lastHE = g.hE;
    lastHM = g.hM;
  }

  const limitFor = (): number => stopT;
  let terminal = status !== 'window_end'; // impact at start

  // first sample at the start
  {
    const p0 = new StepPoly();
    p0.t0 = tStart;
    p0.dt = dir;
    p0.x0.set(ias.x);
    p0.v0.set(ias.v);
    p0.a0.set(ias.a);
    emitSamples(p0, tStart);
  }

  while (!terminal && (limitFor() - ias.t) * dir > 1e-10) {
    if (steps++ >= maxSteps || truncated) {
      truncated = true;
      break;
    }
    const tLimit = limitFor();
    const cap = stepCap(eph, ias.t, ias.x, ias.v, s);
    ias.step(tLimit, cap, poly);
    tNow = poly.t1;

    // ---------- scan this step for events ----------
    const pend: PendingEvent[] = [];
    geoAt(eph, poly, 0, s, g);
    dEs[0] = g.dE; hEs[0] = g.hE; dMs[0] = g.dM; hMs[0] = g.hM;
    geoAt(eph, poly, 1, s, g);
    dEs[NSCAN] = g.dE; hEs[NSCAN] = g.hE; dMs[NSCAN] = g.dM; hMs[NSCAN] = g.hM;
    if (dEs[NSCAN] < minEarthD) {
      minEarthD = dEs[NSCAN];
      minEarthT = poly.t1;
    }
    {
      const sp = Math.max(
        relSpeedAt(poly, 0, 0),
        relSpeedAt(poly, 1, 0),
      ) / AUD_TO_KMS;
      const L = 1.5 * Math.abs(poly.dt) * sp;
      const thr = Math.max(caEarth, hillAu, caMoon + 0.0035) * 1.2;
      const far = (dEs[0] + dEs[NSCAN] - L) / 2 > thr;
      if (!far) {
        for (let m = 1; m < NSCAN; m++) {
          geoAt(eph, poly, m / NSCAN, s, g);
          dEs[m] = g.dE; hEs[m] = g.hE; dMs[m] = g.dM; hMs[m] = g.hM;
        }
        for (let body = 0 as 0 | 1; body < 2; body = (body + 1) as 0 | 1) {
          const D = body === 0 ? dEs : dMs;
          const H = body === 0 ? hEs : hMs;
          const R = radius[body];
          const caThr = body === 0 ? caEarth : caMoon;
          for (let m = 1; m <= NSCAN; m++) {
            const a = (m - 1) / NSCAN;
            const b = m / NSCAN;
            const hPrev = m === 1 ? (body === 0 ? lastHE : lastHM) : H[m - 1];
            let minTau = NaN;
            if (hPrev < 0 && H[m] >= 0) {
              minTau = refineMin(poly, a, b, body);
              const dmin = distAt(poly, minTau, body);
              if (dmin < R && D[m - 1] >= R) {
                const tc = refineCross(poly, a, minTau, body, R);
                pend.push({ tau: tc, kind: 'impact', body, dist: R });
              } else if (dmin < caThr) {
                pend.push({ tau: minTau, kind: 'ca', body, dist: dmin });
              }
            } else if (D[m] < R && D[m - 1] >= R) {
              const tc = refineCross(poly, a, b, body, R);
              pend.push({ tau: tc, kind: 'impact', body, dist: R });
            }
            if (body === 0 && hillNeeded && (D[m - 1] - hillAu) * (D[m] - hillAu) < 0) {
              const tc = refineCross(poly, a, b, 0, hillAu);
              pend.push({ tau: tc, kind: 'hill', body: 0, rising: D[m] > D[m - 1], dist: hillAu });
            }
          }
        }
      }
    }
    lastHE = hEs[NSCAN];
    lastHM = hMs[NSCAN];
    pend.sort((p, q) => p.tau - q.tau);

    // ---------- process events in time order ----------
    let tauStop = 1;
    let tStopStep = poly.t1;
    // a previously declared escape limits the step end
    if ((stopT - poly.t1) * dir < 0) {
      tauStop = (stopT - poly.t0) / poly.dt;
      tStopStep = stopT;
    }
    for (const ev of pend) {
      if (ev.tau > tauStop) break;
      const tev = poly.t0 + poly.dt * ev.tau;
      const jdEv = eph.toJd(tev);
      if (ev.kind === 'impact') {
        pushEvent({
          type: 'impact',
          body: ev.body === 0 ? 'earth' : 'moon',
          jd: jdEv,
          relSpeedKmS: relSpeedAt(poly, ev.tau, ev.body),
          relPosAu: [
            s.X[0] - (ev.body === 0 ? s.E : s.Mo)[0],
            s.X[1] - (ev.body === 0 ? s.E : s.Mo)[1],
            s.X[2] - (ev.body === 0 ? s.E : s.Mo)[2],
          ],
        });
        status = ev.body === 0 ? 'impact_earth' : 'impact_moon';
        terminal = true;
        tauStop = ev.tau;
        tStopStep = tev;
        break;
      } else if (ev.kind === 'ca') {
        const d = ev.dist ?? 0;
        pushEvent({
          type: 'close_approach',
          body: ev.body === 0 ? 'earth' : 'moon',
          jd: jdEv,
          distAu: d,
          distKm: d * AU_KM,
          relSpeedKmS: relSpeedAt(poly, ev.tau, ev.body),
        });
        if (ev.body === 0 && d < minEarthD) {
          minEarthD = d;
          minEarthT = tev;
        }
      } else {
        // hill crossing
        if (!ev.rising) {
          encountered = true;
          pushEvent({ type: 'hill_entry', jd: jdEv, distAu: hillAu });
        } else {
          pushEvent({ type: 'hill_exit', jd: jdEv, distAu: hillAu });
          if (escEnabled && (encountered || !escRequireEncounter) && escapeT === undefined) {
            poly.eval(ev.tau, s.X, s.V);
            evalBodies(eph, tev, s);
            const rx = s.X[0] - s.E[0], ry = s.X[1] - s.E[1], rz = s.X[2] - s.E[2];
            const vx = s.V[0] - s.E[3], vy = s.V[1] - s.E[4], vz = s.V[2] - s.E[5];
            const en = 0.5 * (vx * vx + vy * vy + vz * vz) - muEM / Math.hypot(rx, ry, rz);
            const receding = (rx * vx + ry * vy + rz * vz) * dir > 0;
            if (en > 0 && receding) {
              const la = lookaheadReturn(eph, tev, s.X, s.V, dir, { epsilon, relativity: opts.relativity, returnAu });
              pushEvent({
                type: 'escape_check',
                jd: jdEv,
                result: la.returns ? 'continues' : 'escaped',
                returnJd: la.returns ? la.minJd : undefined,
                returnDistAu: la.returns ? la.minDistAu : undefined,
              });
              if (!la.returns) {
                escapeT = tev;
                stopT = tev + dir * postEscapeDays;
                if ((stopT - tEndUser) * dir > 0) stopT = tEndUser;
                if ((stopT - poly.t1) * dir < 0) {
                  tauStop = (stopT - poly.t0) / poly.dt;
                  tStopStep = stopT;
                }
              }
            } else {
              pushEvent({ type: 'escape_check', jd: jdEv, result: 'continues' });
            }
          }
        }
      }
    }

    // ---------- samples ----------
    emitSamples(poly, tStopStep);
    if (terminal || tauStop < 1) {
      // exact final sample at the impact / at the end of the post-escape tail
      poly.eval(tauStop, s.X, s.V);
      if (Math.abs(lastSampleT - tStopStep) > 1e-12) writeSample(tStopStep, s.X, s.V);
      tNow = tStopStep;
      endedInside = true;
      break;
    }
    if (nBuf >= chunkSamples) yield flushChunk(tNow);
  }

  // ---------- finish ----------
  if (!terminal) {
    status = escapeT !== undefined && !truncated ? 'escaped' : 'window_end';
    if (!endedInside) {
      tNow = ias.t;
      if (!truncated && Math.abs(lastSampleT - ias.t) > 1e-12 && (ias.t - lastSampleT) * dir > 0) {
        writeSample(ias.t, ias.x, ias.v);
      }
    }
  }

  if (nBuf > 0 || pendingEvents.length > 0) {
    yield flushChunk(tNow);
  }

  const summary: SimSummary = {
    status,
    stopJd: eph.toJd(tNow),
    startJd: jdStart,
    direction: dir,
    events: allEvents,
    minEarthDistAu: minEarthD,
    minEarthJd: eph.toJd(minEarthT),
    nSamples: nTotal,
    truncated,
    escapeJd: escapeT !== undefined ? eph.toJd(escapeT) : undefined,
    twoBodyDays: init.twoBodyDays,
    startState,
    stats: { steps: ias.steps, rejected: ias.rejected, evals: ias.evals },
  };
  return summary;
}

/** Run the whole simulation and return all samples. */
export function simulate(asteroid: AsteroidElements | AsteroidState, jdStart: number, opts: SimOptions = {}): SimResult {
  const gen = simulateStream(asteroid, jdStart, opts);
  const parts: Float64Array[] = [];
  let total = 0;
  for (;;) {
    const r = gen.next();
    if (r.done) {
      const samples = new Float64Array(total);
      let o = 0;
      for (const p of parts) {
        samples.set(p, o);
        o += p.length;
      }
      return { ...r.value, samples };
    }
    parts.push(r.value.samples);
    total += r.value.samples.length;
  }
}

/**
 * Hermite-interpolated geocentric state [x,y,z,vx,vy,vz] at jd from a sample array
 * (uses the stored geocentric positions and velocities).  Returns false if jd is outside the samples.
 * Heliocentric position = geocentric + (Earth - Sun) from the ephemeris.
 */
export function sampleStateAt(samples: Float64Array, nSamples: number, jd: number, out: Float64Array): boolean {
  if (nSamples < 2) return false;
  const first = samples[0];
  const last = samples[(nSamples - 1) * SAMPLE_STRIDE];
  const fwd = last >= first;
  const lo = fwd ? first : last;
  const hi = fwd ? last : first;
  if (jd < lo || jd > hi) return false;
  let a = 0;
  let b = nSamples - 1;
  while (b - a > 1) {
    const m = (a + b) >> 1;
    const tm = samples[m * SAMPLE_STRIDE];
    if (fwd ? tm <= jd : tm >= jd) a = m;
    else b = m;
  }
  const oa = a * SAMPLE_STRIDE;
  const ob = b * SAMPLE_STRIDE;
  const h = samples[ob] - samples[oa];
  const sN = h === 0 ? 0 : (jd - samples[oa]) / h;
  const s2 = sN * sN;
  const s3 = s2 * sN;
  const h00 = 2 * s3 - 3 * s2 + 1;
  const h10 = (s3 - 2 * s2 + sN) * h;
  const h01 = -2 * s3 + 3 * s2;
  const h11 = (s3 - s2) * h;
  const d00 = (6 * s2 - 6 * sN) / h;
  const d10 = 3 * s2 - 4 * sN + 1;
  const d01 = (-6 * s2 + 6 * sN) / h;
  const d11 = 3 * s2 - 2 * sN;
  for (let c = 0; c < 3; c++) {
    const p0 = samples[oa + S.GX + c];
    const p1 = samples[ob + S.GX + c];
    const v0 = samples[oa + S.GVX + c];
    const v1 = samples[ob + S.GVX + c];
    out[c] = h00 * p0 + h10 * v0 + h01 * p1 + h11 * v1;
    out[3 + c] = h === 0 ? v0 : d00 * p0 + d10 * v0 + d01 * p1 + d11 * v1;
  }
  return true;
}

