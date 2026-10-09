import { DEG, J2000, solveKepler } from "../kepler";
import { getEphemeris } from "./engine";

// Planet positions for the solar-system view. Analytic JPL "approximate
// positions" mean elements (Standish) are the always-available provider; the
// n-body library's tabulated ephemerides can replace `planetHelio` without
// changing callers.

type Elements = [number, number, number, number, number, number, number, number, number, number, number, number];

interface Planet {
  name: string;
  color: string;
  radius: number; // display radius in the wide view (scene units, not to scale)
  /** a, da, e, de, I, dI, L, dL, long.peri, dlong.peri, long.node, dlong.node (per Julian century). */
  el: Elements;
}

export const PLANETS: Planet[] = [
  { name: "Mercury", color: "#b4a99a", radius: 0.42, el: [0.38709927, 0.00000037, 0.20563593, 0.00001906, 7.00497902, -0.00594749, 252.2503235, 149472.67411175, 77.45779628, 0.16047689, 48.33076593, -0.12534081] },
  { name: "Venus", color: "#e8c98a", radius: 0.65, el: [0.72333566, 0.0000039, 0.00677672, -0.00004107, 3.39467605, -0.0007889, 181.9790995, 58517.81538729, 131.60246718, 0.00268329, 76.67984255, -0.27769418] },
  { name: "Earth", color: "#4da6ff", radius: 0.7, el: [1.00000261, 0.00000562, 0.01671123, -0.00004392, -0.00001531, -0.01294668, 100.46457166, 35999.37306329, 102.93768193, 0.32327364, 0, 0] },
  { name: "Mars", color: "#d4684a", radius: 0.5, el: [1.52371034, 0.00001847, 0.0933941, 0.00007882, 1.84969142, -0.00813131, -4.55343205, 19140.30268499, -23.94362959, 0.44441088, 49.55953891, -0.29257343] },
  { name: "Jupiter", color: "#d9b88f", radius: 1.5, el: [5.202887, -0.00011607, 0.04838624, -0.00013253, 1.30439695, -0.00183714, 34.39644051, 3034.74612775, 14.72847983, 0.21252668, 100.47390909, 0.20469106] },
  { name: "Saturn", color: "#e3d3a1", radius: 1.25, el: [9.53667594, -0.0012506, 0.05386179, -0.00050991, 2.48599187, 0.00193609, 49.95424423, 1222.49362201, 92.59887831, -0.41897216, 113.66242448, -0.28867794] },
  { name: "Uranus", color: "#9fe3e8", radius: 0.9, el: [19.18916464, -0.00196176, 0.04725744, -0.00004397, 0.77263783, -0.00242939, 313.23810451, 428.48202785, 170.9542763, 0.40805281, 74.01692503, 0.04240589] },
  { name: "Neptune", color: "#5b7bf2", radius: 0.9, el: [30.06992276, 0.00026291, 0.00859048, 0.00005105, 1.77004347, 0.00035372, -55.12002969, 218.45945325, 44.96476227, -0.32241464, 131.78422574, -0.00508664] },
];

export const EARTH_PLANET_INDEX = 2;

function orbitElements(p: Planet, T: number) {
  const e = p.el;
  return {
    a: e[0] + e[1] * T,
    e: e[2] + e[3] * T,
    i: (e[4] + e[5] * T) * DEG,
    L: e[6] + e[7] * T,
    varpi: e[8] + e[9] * T,
    om: (e[10] + e[11] * T) * DEG,
  };
}

function rotate(a: number, e: number, inc: number, om: number, w: number, E: number, out: Float64Array | number[], o: number) {
  const xp = a * (Math.cos(E) - e);
  const yp = a * Math.sqrt(1 - e * e) * Math.sin(E);
  const cw = Math.cos(w), sw = Math.sin(w);
  const cO = Math.cos(om), sO = Math.sin(om);
  const ci = Math.cos(inc), si = Math.sin(inc);
  out[o] = xp * (cw * cO - sw * sO * ci) - yp * (sw * cO + cw * sO * ci);
  out[o + 1] = xp * (cw * sO + sw * cO * ci) - yp * (sw * sO - cw * cO * ci);
  out[o + 2] = xp * sw * si + yp * cw * si;
}

/** Index of each PLANETS entry in the N-body ephemeris BODY_ORDER (sun, mercury, venus, earth, moon, mars, ...). */
const EPH_INDEX = [1, 2, 3, 5, 6, 7, 8, 9];
const SUN_INDEX = 0;
const sA = new Float64Array(6);
const sB = new Float64Array(6);

function analyticHelio(i: number, jd: number, out: Float64Array | number[], o: number): void {
  const T = (jd - J2000) / 36525;
  const el = orbitElements(PLANETS[i], T);
  const w = el.varpi * DEG - el.om;
  const M = (el.L - el.varpi) * DEG;
  rotate(el.a, el.e, el.i, el.om, w, solveKepler(M, el.e), out, o);
}

/**
 * Heliocentric ecliptic position (AU) of planet i at Julian Date jd (TDB). Taken from the ephemeris tables
 * (barycentre of the planet minus the Sun) when they are loaded and cover jd; analytic mean elements otherwise.
 * Mars-Neptune are system barycentres, as in the tables.
 */
export function planetHelio(i: number, jd: number, out: Float64Array | number[], o = 0): void {
  const eph = getEphemeris();
  if (eph) {
    const t = eph.toT(jd);
    if (eph.covers(t)) {
      eph.state(EPH_INDEX[i], t, sA);
      eph.state(SUN_INDEX, t, sB);
      out[o] = sA[0] - sB[0];
      out[o + 1] = sA[1] - sB[1];
      out[o + 2] = sA[2] - sB[2];
      return;
    }
  }
  analyticHelio(i, jd, out, o);
}

/** Closed orbit polyline (heliocentric AU) for planet i. */
export function planetOrbit(i: number, jd: number, samples: number): Float64Array {
  const T = (jd - J2000) / 36525;
  const el = orbitElements(PLANETS[i], T);
  const w = el.varpi * DEG - el.om;
  const pts = new Float64Array((samples + 1) * 3);
  for (let s = 0; s <= samples; s++) rotate(el.a, el.e, el.i, el.om, w, (s / samples) * Math.PI * 2, pts, s * 3);
  return pts;
}

// Compressed radial scale for the wide view only (Mercury ~10, Earth ~19,
// Jupiter ~39, Neptune ~62 scene units).
export function wideRadius(rAu: number): number {
  return 14 * Math.log(1 + rAu / 0.35);
}

/** Heliocentric ecliptic AU -> wide-view scene coordinates (ecliptic = XZ plane). */
export function toWide(x: number, y: number, z: number, out: Float32Array | number[], o = 0): void {
  const r = Math.sqrt(x * x + y * y + z * z);
  const f = r < 1e-12 ? 0 : wideRadius(r) / r;
  out[o] = x * f;
  out[o + 1] = z * f;
  out[o + 2] = -y * f;
}
