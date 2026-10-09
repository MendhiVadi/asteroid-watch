// Renderer-friendly trajectory arrays built from raw simulation samples.

import { SAMPLE_STRIDE, S } from './sim.ts';
import type { SimEvent, SimSummary, TerminalStatus } from './sim.ts';
import { Ephemeris } from './ephemeris.ts';

/** Struct-of-arrays piece of a trajectory (one streamed chunk, or the whole thing). */
export interface TrajectoryArrays {
  /** sample times, JD (TDB) */
  jd: Float64Array;
  /** asteroid position relative to Earth's centre, ecliptic J2000, AU: [x0,y0,z0,x1,...] */
  earth: Float32Array;
  /** asteroid heliocentric position, ecliptic J2000, AU */
  helio: Float32Array;
  /** Moon position relative to Earth's centre at the same times (AU); only if requested */
  moon?: Float32Array;
  /** |asteroid - Earth centre| and |asteroid - Moon| (AU) */
  dEarth: Float32Array;
  dMoon: Float32Array;
}

export interface Trajectory extends TrajectoryArrays {
  status: TerminalStatus;
  startJd: number;
  /** time at which the trajectory ends (impact, end of post-escape tail, or window end) */
  stopJd: number;
  direction: 1 | -1;
  events: SimEvent[];
  minEarthDistAu: number;
  minEarthJd: number;
  escapeJd?: number;
  truncated: boolean;
  count: number;
  summary: SimSummary;
}

export function samplesToArrays(samples: Float64Array, eph: Ephemeris | null, includeMoon: boolean): TrajectoryArrays {
  const n = samples.length / SAMPLE_STRIDE;
  const jd = new Float64Array(n);
  const earth = new Float32Array(3 * n);
  const helio = new Float32Array(3 * n);
  const dEarth = new Float32Array(n);
  const dMoon = new Float32Array(n);
  const moon = includeMoon && eph ? new Float32Array(3 * n) : undefined;
  const tmp = new Float64Array(6);
  for (let i = 0; i < n; i++) {
    const o = i * SAMPLE_STRIDE;
    jd[i] = samples[o + S.JD];
    for (let c = 0; c < 3; c++) {
      earth[3 * i + c] = samples[o + S.GX + c];
      helio[3 * i + c] = samples[o + S.HX + c];
    }
    dEarth[i] = samples[o + S.DE];
    dMoon[i] = samples[o + S.DM];
    if (moon && eph) {
      eph.moonRelState(eph.toT(jd[i]), tmp);
      moon[3 * i] = tmp[0];
      moon[3 * i + 1] = tmp[1];
      moon[3 * i + 2] = tmp[2];
    }
  }
  return { jd, earth, helio, moon, dEarth, dMoon };
}

export function concatArrays(parts: TrajectoryArrays[]): TrajectoryArrays {
  const n = parts.reduce((a, p) => a + p.jd.length, 0);
  const hasMoon = parts.length > 0 && parts.every((p) => p.moon);
  const out: TrajectoryArrays = {
    jd: new Float64Array(n),
    earth: new Float32Array(3 * n),
    helio: new Float32Array(3 * n),
    moon: hasMoon ? new Float32Array(3 * n) : undefined,
    dEarth: new Float32Array(n),
    dMoon: new Float32Array(n),
  };
  let o = 0;
  for (const p of parts) {
    out.jd.set(p.jd, o);
    out.earth.set(p.earth, 3 * o);
    out.helio.set(p.helio, 3 * o);
    out.dEarth.set(p.dEarth, o);
    out.dMoon.set(p.dMoon, o);
    if (out.moon && p.moon) out.moon.set(p.moon, 3 * o);
    o += p.jd.length;
  }
  return out;
}

export function makeTrajectory(arrays: TrajectoryArrays, summary: SimSummary): Trajectory {
  return {
    ...arrays,
    status: summary.status,
    startJd: summary.startJd,
    stopJd: summary.stopJd,
    direction: summary.direction,
    events: summary.events,
    minEarthDistAu: summary.minEarthDistAu,
    minEarthJd: summary.minEarthJd,
    escapeJd: summary.escapeJd,
    truncated: summary.truncated,
    count: arrays.jd.length,
    summary,
  };
}
