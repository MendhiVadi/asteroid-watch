// Osculating orbital elements <-> Cartesian state (heliocentric, ecliptic J2000).

import { DEG, TWO_PI } from './constants.ts';

export interface AsteroidElements {
  /** semi-major axis, AU (negative for hyperbolic orbits, SBDB convention) */
  a: number;
  e: number;
  /** degrees */
  i: number;
  om: number;
  w: number;
  ma: number;
  /** epoch of the elements, JD (TDB) */
  epoch: number;
  id?: string;
  name?: string;
}

/** Barycentric cartesian state of the asteroid at jd. */
export interface AsteroidState {
  jd: number;
  x: ArrayLike<number>;
  v: ArrayLike<number>;
  id?: string;
  name?: string;
}

export function isState(a: AsteroidElements | AsteroidState): a is AsteroidState {
  return (a as AsteroidState).x !== undefined;
}

/**
 * Two-body (Sun only) re-epoching of osculating elements: same orbit, new epoch (JD).
 * Only the mean anomaly changes.  mu = GM_sun in AU^3/day^2.
 */
export function propagateElementsTwoBody(el: AsteroidElements, newEpochJd: number, mu: number): AsteroidElements {
  const A = Math.abs(el.a);
  const n = Math.sqrt(mu / (A * A * A)); // rad/day
  let ma = el.ma + (n * (newEpochJd - el.epoch)) / DEG;
  if (el.e < 1) {
    ma %= 360;
    if (ma < 0) ma += 360;
  }
  return { ...el, ma, epoch: newEpochJd };
}

export function solveKeplerElliptic(M0: number, e: number): number {
  let M = M0 % TWO_PI;
  if (M > Math.PI) M -= TWO_PI;
  else if (M < -Math.PI) M += TWO_PI;
  let E = e < 0.8 ? M + e * Math.sin(M) : M >= 0 ? Math.PI : -Math.PI;
  for (let n = 0; n < 60; n++) {
    const d = (E - e * Math.sin(E) - M) / (1 - e * Math.cos(E));
    E -= d;
    if (Math.abs(d) < 1e-15) break;
  }
  return E;
}

export function solveKeplerHyperbolic(M: number, e: number): number {
  let H = Math.asinh(M / e);
  if (Math.abs(M) > 5) H = Math.sign(M) * Math.log((2 * Math.abs(M)) / e + 1.8);
  for (let n = 0; n < 100; n++) {
    const d = (e * Math.sinh(H) - H - M) / (e * Math.cosh(H) - 1);
    H -= d;
    if (Math.abs(d) < 1e-15 * Math.max(1, Math.abs(H))) break;
  }
  return H;
}

/**
 * Heliocentric ecliptic-J2000 state at the epoch of the elements.
 * mu = GM of the Sun in AU^3/day^2.
 */
export function elementsToState(el: AsteroidElements, mu: number): { x: [number, number, number]; v: [number, number, number] } {
  const { a, e } = el;
  const inc = el.i * DEG;
  const om = el.om * DEG;
  const w = el.w * DEG;
  const M = el.ma * DEG;
  let xp: number, yp: number, vxp: number, vyp: number;
  if (e < 1) {
    const E = solveKeplerElliptic(M, e);
    const sE = Math.sin(E);
    const cE = Math.cos(E);
    const b = a * Math.sqrt(1 - e * e);
    const n = Math.sqrt(mu / (a * a * a));
    const den = 1 - e * cE;
    xp = a * (cE - e);
    yp = b * sE;
    vxp = (-a * n * sE) / den;
    vyp = (b * n * cE) / den;
  } else if (e > 1) {
    const A = Math.abs(a);
    const H = solveKeplerHyperbolic(M, e);
    const sh = Math.sinh(H);
    const ch = Math.cosh(H);
    const b = A * Math.sqrt(e * e - 1);
    const n = Math.sqrt(mu / (A * A * A));
    const den = e * ch - 1;
    xp = A * (e - ch);
    yp = b * sh;
    vxp = (-A * n * sh) / den;
    vyp = (b * n * ch) / den;
  } else {
    throw new Error('elementsToState: parabolic orbits (e = 1) are not supported');
  }
  const cO = Math.cos(om), sO = Math.sin(om);
  const cw = Math.cos(w), sw = Math.sin(w);
  const ci = Math.cos(inc), si = Math.sin(inc);
  const Px = cw * cO - sw * sO * ci;
  const Py = cw * sO + sw * cO * ci;
  const Pz = sw * si;
  const Qx = -sw * cO - cw * sO * ci;
  const Qy = -sw * sO + cw * cO * ci;
  const Qz = cw * si;
  return {
    x: [xp * Px + yp * Qx, xp * Py + yp * Qy, xp * Pz + yp * Qz],
    v: [vxp * Px + vyp * Qx, vxp * Py + vyp * Qy, vxp * Pz + vyp * Qz],
  };
}

/** Classical elements of a state relative to a central mass mu (angles in degrees). */
export function stateToElements(
  x: ArrayLike<number>,
  v: ArrayLike<number>,
  mu: number,
): { a: number; e: number; i: number; om: number; w: number; nu: number; energy: number; h: number } {
  const r = Math.hypot(x[0], x[1], x[2]);
  const v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  const hx = x[1] * v[2] - x[2] * v[1];
  const hy = x[2] * v[0] - x[0] * v[2];
  const hz = x[0] * v[1] - x[1] * v[0];
  const h = Math.hypot(hx, hy, hz);
  const rv = x[0] * v[0] + x[1] * v[1] + x[2] * v[2];
  const ex = ((v2 - mu / r) * x[0] - rv * v[0]) / mu;
  const ey = ((v2 - mu / r) * x[1] - rv * v[1]) / mu;
  const ez = ((v2 - mu / r) * x[2] - rv * v[2]) / mu;
  const e = Math.hypot(ex, ey, ez);
  const energy = v2 / 2 - mu / r;
  const a = -mu / (2 * energy);
  const inc = Math.acos(hz / h);
  const nx = -hy;
  const ny = hx;
  const nn = Math.hypot(nx, ny);
  let om = nn > 0 ? Math.atan2(ny, nx) : 0;
  if (om < 0) om += TWO_PI;
  let w = nn > 0 ? Math.acos(Math.max(-1, Math.min(1, (nx * ex + ny * ey) / (nn * e)))) : Math.atan2(ey, ex);
  if (ez < 0) w = TWO_PI - w;
  let nu = Math.acos(Math.max(-1, Math.min(1, (ex * x[0] + ey * x[1] + ez * x[2]) / (e * r))));
  if (rv < 0) nu = TWO_PI - nu;
  return { a, e, i: inc / DEG, om: om / DEG, w: w / DEG, nu: nu / DEG, energy, h };
}
