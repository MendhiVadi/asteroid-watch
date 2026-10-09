// Hermite-interpolated planetary / lunar ephemeris built from the binary tables
// written by scripts/fetch_ephemeris.py.  Pure data + maths, no DOM / node APIs.
//
// Interpolation: cubic Hermite on (position, velocity) node pairs.
//   planets : 1-day nodes, Float64.  Worst-case Earth error is ~0.1 km (measured vs Horizons, see the
//             manifest 'validation' block; dominated by the monthly lunar wobble of the Earth about the EMB).
//   moon    : stored geocentric on nodes spaced manifest.moon.stepDays apart (Float64; any step is
//             handled - currently 12 h, ~0.3 km worst-case position error vs Horizons, see the
//             manifest 'validation' block); barycentric Moon = Earth(t) + MoonRel(t).
//
// Coverage is whatever the manifest says (currently 2020-01-01 .. 2100-12-31, JD 2458849.5 .. 2488433.5).

export interface EphemerisManifest {
  version: number;
  auKm: number;
  cKmS: number;
  cAuDay: number;
  emrat: number;
  main: {
    file: string;
    dtype: 'float64';
    jd0: number;
    stepDays: number;
    count: number;
    bodies: string[];
  };
  moon: {
    file: string;
    dtype: 'float64' | 'float32';
    jd0: number;
    stepDays: number;
    count: number;
    relativeTo: 'earth';
  };
  gm: Record<string, { au3d2: number; km3s2: number }>;
  radiiKm: { earth: number; moon: number; sun: number };
}

/**
 * Order of the point masses used by the force model.
 * Index constants below refer to this order.
 */
export const BODY_ORDER = [
  'sun',
  'mercury',
  'venus',
  'earth',
  'moon',
  'mars',
  'jupiter',
  'saturn',
  'uranus',
  'neptune',
] as const;
export type BodyName = (typeof BODY_ORDER)[number];
export const N_BODIES = BODY_ORDER.length;
export const IDX_SUN = 0;
export const IDX_EARTH = 3;
export const IDX_MOON = 4;

export class Ephemeris {
  readonly manifest: EphemerisManifest;
  /** Coverage in internal time t (days since jd0): [tStart, tEnd]. */
  readonly tStart: number;
  readonly tEnd: number;
  /** GM of each body in BODY_ORDER, AU^3/day^2. */
  readonly gm: Float64Array;
  readonly cAuDay: number;
  private readonly main: Float64Array;
  private readonly moon: Float64Array | Float32Array;
  private readonly nMain: number;
  private readonly nMoon: number;
  /** JD (TDB) of t = 0.  Internal time `t` is days since jd0 (keeps full double resolution). */
  readonly jd0: number;
  private readonly moonT0: number;
  private readonly moonStep: number;
  /** body -> index into the [body] dimension of the main table (-1 for moon) */
  private readonly mainIndex: Int32Array;

  constructor(manifest: EphemerisManifest, planets: ArrayBuffer, moon: ArrayBuffer) {
    this.manifest = manifest;
    this.main = new Float64Array(planets);
    this.moon = manifest.moon.dtype === 'float32' ? new Float32Array(moon) : new Float64Array(moon);
    this.nMain = manifest.main.count;
    this.nMoon = manifest.moon.count;
    this.jd0 = manifest.main.jd0;
    this.moonT0 = manifest.moon.jd0 - manifest.main.jd0;
    this.moonStep = manifest.moon.stepDays;
    if (manifest.main.stepDays !== 1) throw new Error('ephemeris: expected 1-day planetary nodes');
    if (!(this.moonStep > 0 && this.moonStep <= 1) || this.nMoon < 2 || this.nMain < 2) {
      throw new Error(`ephemeris: bad moon grid (step ${this.moonStep} d, ${this.nMoon} nodes)`);
    }
    if (this.main.length !== manifest.main.bodies.length * this.nMain * 6) {
      throw new Error('ephemeris: planets.bin size does not match manifest');
    }
    if (this.moon.length !== this.nMoon * 6) throw new Error('ephemeris: moon.bin size does not match manifest');
    this.tStart = Math.max(0, this.moonT0);
    this.tEnd = Math.min(this.nMain - 1, this.moonT0 + (this.nMoon - 1) * this.moonStep);
    this.cAuDay = manifest.cAuDay;

    this.mainIndex = new Int32Array(N_BODIES);
    this.gm = new Float64Array(N_BODIES);
    for (let k = 0; k < N_BODIES; k++) {
      const name = BODY_ORDER[k];
      this.mainIndex[k] = name === 'moon' ? -1 : manifest.main.bodies.indexOf(name);
      if (name !== 'moon' && this.mainIndex[k] < 0) throw new Error(`ephemeris: body ${name} missing`);
      const g = manifest.gm[name];
      if (!g) throw new Error(`ephemeris: GM for ${name} missing`);
      this.gm[k] = g.au3d2;
    }
  }

  /** JD (TDB) -> internal time t (days since jd0). Exact for JDs near the coverage window. */
  toT(jd: number): number {
    return jd - this.jd0;
  }

  /** Internal time t -> JD (TDB). */
  toJd(t: number): number {
    return this.jd0 + t;
  }

  get jdStart(): number {
    return this.jd0 + this.tStart;
  }

  get jdEnd(): number {
    return this.jd0 + this.tEnd;
  }

  covers(t: number): boolean {
    return t >= this.tStart && t <= this.tEnd;
  }

  /**
   * All time arguments below are internal time t (days since jd0), see toT().
   * State of body k (index into BODY_ORDER) at t: out[0..2] = position, out[3..5] = velocity
   * (barycentric, ecliptic J2000, AU and AU/day).  Velocity is the derivative of the
   * interpolating polynomial.
   */
  state(k: number, t: number, out: Float64Array): void {
    if (k === IDX_MOON) {
      this.mainState(this.mainIndex[IDX_EARTH], t, out);
      this.moonRelState(t, M_TMP);
      for (let c = 0; c < 6; c++) out[c] += M_TMP[c];
    } else {
      this.mainState(this.mainIndex[k], t, out);
    }
  }

  /** Positions only for all bodies in BODY_ORDER: out[3k..3k+2]. */
  positions(t: number, out: Float64Array): void {
    const u = t;
    if (!(u >= 0 && u <= this.nMain - 1)) throw new RangeError(`ephemeris: t=${t} outside coverage`);
    let i = Math.floor(u);
    if (i > this.nMain - 2) i = this.nMain - 2;
    const s = u - i;
    const s2 = s * s;
    const s3 = s2 * s;
    const h00 = 2 * s3 - 3 * s2 + 1;
    const h10 = s3 - 2 * s2 + s;
    const h01 = -2 * s3 + 3 * s2;
    const h11 = s3 - s2;
    const main = this.main;
    const stride = this.nMain * 6;
    for (let k = 0; k < N_BODIES; k++) {
      const mi = this.mainIndex[k];
      if (mi < 0) continue;
      const b0 = mi * stride + i * 6;
      const b1 = b0 + 6;
      const o = 3 * k;
      out[o] = h00 * main[b0] + h10 * main[b0 + 3] + h01 * main[b1] + h11 * main[b1 + 3];
      out[o + 1] = h00 * main[b0 + 1] + h10 * main[b0 + 4] + h01 * main[b1 + 1] + h11 * main[b1 + 4];
      out[o + 2] = h00 * main[b0 + 2] + h10 * main[b0 + 5] + h01 * main[b1 + 2] + h11 * main[b1 + 5];
    }
    // Moon = Earth + geocentric Moon
    this.moonRelState(t, M_TMP);
    const e = 3 * IDX_EARTH;
    const m = 3 * IDX_MOON;
    out[m] = out[e] + M_TMP[0];
    out[m + 1] = out[e + 1] + M_TMP[1];
    out[m + 2] = out[e + 2] + M_TMP[2];
  }

  /**
   * Same as positions() at t = tRef + dtOff, but split as p = base + corr with base = the table node value at
   * floor(t) and corr = the (small) interpolation increment.  Lets the force model form
   * (xRef - base) + (dx - corr) without cancellation noise when the particle is close to a body.
   * Writes base[3k..], corr[3k..].
   */
  positionsSplit(tRef: number, dtOff: number, base: Float64Array, corr: Float64Array): void {
    const u = tRef + dtOff;
    if (!(u >= 0 && u <= this.nMain - 1)) throw new RangeError(`ephemeris: t=${u} outside coverage`);
    let i = Math.floor(tRef);
    let s = tRef - i + dtOff; // keeps full precision of the offset (dtOff may be ~1e-9 d)
    const carry = Math.floor(s);
    i += carry;
    s -= carry;
    if (i > this.nMain - 2) {
      s += i - (this.nMain - 2);
      i = this.nMain - 2;
    }
    const s2 = s * s;
    const s3 = s2 * s;
    const h01 = -2 * s3 + 3 * s2;
    const h10 = s3 - 2 * s2 + s;
    const h11 = s3 - s2;
    const main = this.main;
    const stride = this.nMain * 6;
    for (let k = 0; k < N_BODIES; k++) {
      const mi = this.mainIndex[k];
      if (mi < 0) continue;
      const b0 = mi * stride + i * 6;
      const b1 = b0 + 6;
      const o = 3 * k;
      for (let c = 0; c < 3; c++) {
        base[o + c] = main[b0 + c];
        corr[o + c] = h01 * (main[b1 + c] - main[b0 + c]) + h10 * main[b0 + 3 + c] + h11 * main[b1 + 3 + c];
      }
    }
    this.moonRelState(u, M_TMP);
    const e = 3 * IDX_EARTH;
    const m = 3 * IDX_MOON;
    for (let c = 0; c < 3; c++) {
      base[m + c] = base[e + c];
      corr[m + c] = corr[e + c] + M_TMP[c];
    }
  }

  /** Sun barycentric state convenience wrapper. */
  sunState(t: number, out: Float64Array): void {
    this.mainState(this.mainIndex[IDX_SUN], t, out);
  }

  private mainState(mi: number, t: number, out: Float64Array): void {
    const u = t;
    if (!(u >= 0 && u <= this.nMain - 1)) throw new RangeError(`ephemeris: t=${t} outside coverage`);
    let i = Math.floor(u);
    if (i > this.nMain - 2) i = this.nMain - 2;
    hermite6(this.main, mi * this.nMain * 6 + i * 6, u - i, 1, out);
  }

  /** Geocentric Moon state (cubic Hermite on the manifest's moon grid, e.g. 12-hour nodes). */
  moonRelState(t: number, out: Float64Array): void {
    const u = (t - this.moonT0) / this.moonStep;
    if (!(u >= 0 && u <= this.nMoon - 1)) throw new RangeError(`ephemeris: t=${t} outside coverage`);
    let i = Math.floor(u);
    if (i > this.nMoon - 2) i = this.nMoon - 2;
    hermite6(this.moon, i * 6, u - i, this.moonStep, out);
  }
}

const M_TMP = new Float64Array(6);

/** Cubic Hermite for 6-vectors [x,y,z,vx,vy,vz] stored at base and base+6; h = node spacing (days). */
function hermite6(tab: Float64Array | Float32Array, base: number, s: number, h: number, out: Float64Array): void {
  const s2 = s * s;
  const s3 = s2 * s;
  const h00 = 2 * s3 - 3 * s2 + 1;
  const h10 = (s3 - 2 * s2 + s) * h;
  const h01 = -2 * s3 + 3 * s2;
  const h11 = (s3 - s2) * h;
  const d00 = (6 * s2 - 6 * s) / h;
  const d10 = 3 * s2 - 4 * s + 1;
  const d01 = (-6 * s2 + 6 * s) / h;
  const d11 = 3 * s2 - 2 * s;
  const b1 = base + 6;
  for (let c = 0; c < 3; c++) {
    const p0 = tab[base + c];
    const v0 = tab[base + 3 + c];
    const p1 = tab[b1 + c];
    const v1 = tab[b1 + 3 + c];
    out[c] = h00 * p0 + h10 * v0 + h01 * p1 + h11 * v1;
    out[3 + c] = d00 * p0 + d10 * v0 + d01 * p1 + d11 * v1;
  }
}

/** Load an Ephemeris from a base URL (works in the main thread and in workers). */
export async function loadEphemeris(baseUrl = '/data/ephemeris/'): Promise<Ephemeris> {
  const base = baseUrl.endsWith('/') ? baseUrl : baseUrl + '/';
  const mRes = await fetch(base + 'manifest.json', { cache: 'no-cache' });
  if (!mRes.ok) throw new Error(`ephemeris manifest: HTTP ${mRes.status}`);
  const manifest = (await mRes.json()) as EphemerisManifest;
  // The query string is ignored by static hosts but changes whenever the grids change, so a stale cached .bin
  // can never be paired with a newer manifest.
  const v = `?v=${manifest.main.count}-${manifest.moon.count}-${manifest.moon.stepDays}`;
  const [pRes, lRes] = await Promise.all([fetch(base + manifest.main.file + v), fetch(base + manifest.moon.file + v)]);
  if (!pRes.ok || !lRes.ok) throw new Error(`ephemeris binaries: HTTP ${pRes.status}/${lRes.status}`);
  const [p, l] = await Promise.all([pRes.arrayBuffer(), lRes.arrayBuffer()]);
  return new Ephemeris(manifest, p, l);
}
