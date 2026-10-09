// Keplerian orbit maths shared by the main thread and the position worker.
// Units: AU, days, radians internally. Input element angles are degrees and
// the epoch is a Julian Date (the JPL SBDB convention used by the data contract).

export const DEG = Math.PI / 180;
export const TWO_PI = Math.PI * 2;
export const J2000 = 2451545.0;
/** Gaussian gravitational constant (rad/day for a = 1 AU). */
export const GAUSS_K = 0.01720209895;
export const AU_KM = 149597870.7;
export const LUNAR_DISTANCE_AU = 0.00256955;
export const EARTH_RADIUS_SCENE = 2;

export interface OrbitalElements {
  a: number;
  e: number;
  i: number;
  om: number;
  w: number;
  ma: number;
  epoch: number;
}

export function jdFromMs(ms: number): number {
  return ms / 86400000 + 2440587.5;
}

export function msFromJd(jd: number): number {
  return (jd - 2440587.5) * 86400000;
}

/** Solve Kepler's equation M = E - e sin E by Newton-Raphson. */
export function solveKepler(meanAnomaly: number, e: number): number {
  let M = meanAnomaly % TWO_PI;
  if (M > Math.PI) M -= TWO_PI;
  else if (M < -Math.PI) M += TWO_PI;
  let E = e < 0.8 ? M + e * Math.sin(M) : M >= 0 ? Math.PI : -Math.PI;
  for (let n = 0; n < 30; n++) {
    const s = Math.sin(E);
    const c = Math.cos(E);
    const d = (E - e * s - M) / (1 - e * c);
    E -= d;
    if (Math.abs(d) < 1e-11) break;
  }
  return E;
}

// ---------------------------------------------------------------------------
// Packed element table: per-object constants precomputed once so the per-frame
// work is a Newton solve plus two multiplies, not a full trig chain.
// layout: a, e, n, M0, epoch, b, Px, Py, Pz, Qx, Qy, Qz
// ---------------------------------------------------------------------------
export const ELEM_STRIDE = 12;

export function buildElementTable(rows: readonly OrbitalElements[]): Float64Array {
  const t = new Float64Array(rows.length * ELEM_STRIDE);
  for (let k = 0; k < rows.length; k++) {
    const r = rows[k];
    const a = Math.max(r.a, 1e-3);
    const e = Math.min(Math.max(r.e, 0), 0.9999);
    const inc = r.i * DEG;
    const om = r.om * DEG;
    const w = r.w * DEG;
    const cO = Math.cos(om), sO = Math.sin(om);
    const cw = Math.cos(w), sw = Math.sin(w);
    const ci = Math.cos(inc), si = Math.sin(inc);
    const b = k * ELEM_STRIDE;
    t[b] = a;
    t[b + 1] = e;
    t[b + 2] = GAUSS_K / Math.pow(a, 1.5);
    t[b + 3] = r.ma * DEG;
    t[b + 4] = r.epoch;
    t[b + 5] = a * Math.sqrt(1 - e * e);
    t[b + 6] = cw * cO - sw * sO * ci;
    t[b + 7] = cw * sO + sw * cO * ci;
    t[b + 8] = sw * si;
    t[b + 9] = -sw * cO - cw * sO * ci;
    t[b + 10] = -sw * sO + cw * cO * ci;
    t[b + 11] = cw * si;
  }
  return t;
}

/** Heliocentric ecliptic position (AU) of table entry k at Julian Date jd. */
export function helioPosition(t: Float64Array, k: number, jd: number, out: Float64Array | number[], o = 0): void {
  const b = k * ELEM_STRIDE;
  const e = t[b + 1];
  const E = solveKepler(t[b + 3] + t[b + 2] * (jd - t[b + 4]), e);
  const xp = t[b] * (Math.cos(E) - e);
  const yp = t[b + 5] * Math.sin(E);
  out[o] = xp * t[b + 6] + yp * t[b + 9];
  out[o + 1] = xp * t[b + 7] + yp * t[b + 10];
  out[o + 2] = xp * t[b + 8] + yp * t[b + 11];
}

/**
 * Earth (strictly the Earth-Moon barycentre) from the JPL "approximate
 * positions" J2000 mean elements; ample for visualisation.
 */
export function earthHelio(jd: number, out: Float64Array | number[], o = 0): void {
  const T = (jd - J2000) / 36525;
  const a = 1.00000261 + 0.00000562 * T;
  const e = 0.01671123 - 0.00004392 * T;
  const inc = (-0.00001531 - 0.01294668 * T) * DEG;
  // Table 1 (1800-2050) EM-barycentre mean-longitude rate, same value as the Python reference.
  const L = 100.46457166 + 35999.37244981 * T;
  const varpi = 102.93768193 + 0.32327364 * T;
  const om = 0;
  const w = (varpi - om) * DEG;
  const M = (L - varpi) * DEG;
  const E = solveKepler(M, e);
  const xp = a * (Math.cos(E) - e);
  const yp = a * Math.sqrt(1 - e * e) * Math.sin(E);
  const cw = Math.cos(w), sw = Math.sin(w);
  const ci = Math.cos(inc), si = Math.sin(inc);
  // Node longitude is 0 for these elements, so cos/sin of the node are 1/0.
  out[o] = xp * cw - yp * sw;
  out[o + 1] = (xp * sw + yp * cw) * ci;
  out[o + 2] = (xp * sw + yp * cw) * si;
}

// ---------------------------------------------------------------------------
// Radial compression. Real distances span ~1e-4 .. 6+ AU while Earth is 4e-5 AU
// across, so offsets from Earth are mapped through R(r) = R0 + K log10(1 + r/S):
// monotone, so ordering by distance is preserved, close objects visibly roam
// around Earth and distant ones still spread out when zoomed out.
// ---------------------------------------------------------------------------
const R0 = 3.4;
const K = 5.6;
const S = 0.001;

export function compressRadius(rAu: number): number {
  return R0 + K * Math.log10(1 + rAu / S);
}

/** Scene radius of the Moon's orbit (consistent with the asteroid mapping). */
export const MOON_ORBIT_SCENE = compressRadius(LUNAR_DISTANCE_AU);

/** Ecliptic Earth-centred offset (AU) -> compressed scene coordinates (ecliptic = XZ plane). */
export function offsetToScene(dx: number, dy: number, dz: number, out: Float32Array | number[], o = 0): void {
  const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
  if (r < 1e-12) {
    out[o] = R0;
    out[o + 1] = 0;
    out[o + 2] = 0;
    return;
  }
  const f = compressRadius(r) / r;
  out[o] = dx * f;
  out[o + 1] = dz * f;
  out[o + 2] = -dy * f;
}

const earth = new Float64Array(3);
const hp = new Float64Array(3);

/** Batch update of all scene positions (x,y,z interleaved) for Julian Date jd. */
export function computeScenePositions(table: Float64Array, count: number, jd: number, out: Float32Array): void {
  earthHelio(jd, earth);
  const ex = earth[0], ey = earth[1], ez = earth[2];
  for (let k = 0; k < count; k++) {
    helioPosition(table, k, jd, hp);
    offsetToScene(hp[0] - ex, hp[1] - ey, hp[2] - ez, out, k * 3);
  }
}

/** True Earth-asteroid distance in AU. */
export function earthDistanceAu(table: Float64Array, k: number, jd: number): number {
  earthHelio(jd, earth);
  helioPosition(table, k, jd, hp);
  return Math.hypot(hp[0] - earth[0], hp[1] - earth[1], hp[2] - earth[2]);
}

/** Earth-relative path of table entry k, sampled into scene coordinates. */
export function samplePath(
  table: Float64Array,
  k: number,
  jdStart: number,
  jdEnd: number,
  samples: number,
  out: Float32Array,
): void {
  for (let s = 0; s < samples; s++) {
    const jd = jdStart + ((jdEnd - jdStart) * s) / (samples - 1);
    earthHelio(jd, earth);
    helioPosition(table, k, jd, hp);
    offsetToScene(hp[0] - earth[0], hp[1] - earth[1], hp[2] - earth[2], out, s * 3);
  }
}

/** Moon's mean ecliptic longitude (radians) - low-order series, fine for display. */
export function moonLongitude(jd: number): number {
  return (218.316 + 13.176396 * (jd - J2000)) * DEG;
}

// ---------------------------------------------------------------------------
// Cinematic mapping: true scale near Earth (so impacts actually reach the
// surface), then logarithmic. R(rE) = Earth radius in scene units.
// ---------------------------------------------------------------------------
export const EARTH_RADIUS_AU = 4.2635e-5;
const CIN_LINEAR_LIMIT = 6 * EARTH_RADIUS_AU;

export function cinemaRadius(rAu: number): number {
  const k = EARTH_RADIUS_SCENE / EARTH_RADIUS_AU;
  if (rAu <= CIN_LINEAR_LIMIT) return rAu * k;
  return 6 * EARTH_RADIUS_SCENE + 6 * EARTH_RADIUS_SCENE * Math.log(rAu / CIN_LINEAR_LIMIT);
}

export function cinemaOffsetToScene(dx: number, dy: number, dz: number, out: Float32Array | number[], o = 0): void {
  const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
  if (r < 1e-14) {
    out[o] = out[o + 1] = out[o + 2] = 0;
    return;
  }
  const f = cinemaRadius(r) / r;
  out[o] = dx * f;
  out[o + 1] = dz * f;
  out[o + 2] = -dy * f;
}

/** Approximate Earth-centred ecliptic Moon position (AU). */
export function moonEarthRelative(jd: number, out: Float64Array | number[], o = 0): void {
  const L = moonLongitude(jd);
  const inc = 5.145 * DEG;
  const r = LUNAR_DISTANCE_AU;
  out[o] = r * Math.cos(L);
  out[o + 1] = r * Math.sin(L) * Math.cos(inc);
  out[o + 2] = r * Math.sin(L) * Math.sin(inc);
}

/** Earth heliocentric position (ecliptic AU) - exposed for sun direction. */
export function earthHelioPosition(jd: number): [number, number, number] {
  const o = new Float64Array(3);
  earthHelio(jd, o);
  return [o[0], o[1], o[2]];
}
