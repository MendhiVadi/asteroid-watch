// Force model: massless test particle in the barycentric frame under the point-mass
// gravity of Sun, Mercury, Venus, Earth, Moon, Mars, Jupiter, Saturn, Uranus, Neptune
// (positions from the Hermite ephemeris) plus an optional Sun-only 1PN term.

import { Ephemeris, IDX_SUN, N_BODIES } from './ephemeris.ts';
import type { AccFn } from './ias15.ts';

export interface ForceOptions {
  /** Include the Sun's first post-Newtonian (Schwarzschild, beta = gamma = 1) term. Default true. */
  relativity?: boolean;
  /** Include these bodies (names from BODY_ORDER). Default: all. Used for tests. */
  mask?: boolean[];
}

export function makeForce(eph: Ephemeris, opts: ForceOptions = {}): AccFn {
  const rel = opts.relativity ?? true;
  const gm = eph.gm;
  const c2 = eph.cAuDay * eph.cAuDay;
  const base = new Float64Array(3 * N_BODIES);
  const corr = new Float64Array(3 * N_BODIES);
  const sunS = new Float64Array(6);
  const use = opts.mask;
  const gmSun = gm[IDX_SUN];
  return (t, x, v, a, xRef, dx, tRef, dtOff) => {
    eph.positionsSplit(tRef, dtOff, base, corr);
    let ax = 0;
    let ay = 0;
    let az = 0;
    for (let k = 0; k < N_BODIES; k++) {
      if (use && !use[k]) continue;
      // body - particle, formed as (base - xRef) + (corr - dx) to avoid cancellation noise
      const ex = base[3 * k] - xRef[0] + (corr[3 * k] - dx[0]);
      const ey = base[3 * k + 1] - xRef[1] + (corr[3 * k + 1] - dx[1]);
      const ez = base[3 * k + 2] - xRef[2] + (corr[3 * k + 2] - dx[2]);
      const r2 = ex * ex + ey * ey + ez * ez;
      const f = gm[k] / (r2 * Math.sqrt(r2));
      ax += f * ex;
      ay += f * ey;
      az += f * ez;
    }
    if (rel) {
      eph.sunState(t, sunS);
      const rx = x[0] - base[0] - corr[0];
      const ry = x[1] - base[1] - corr[1];
      const rz = x[2] - base[2] - corr[2];
      const ux = v[0] - sunS[3];
      const uy = v[1] - sunS[4];
      const uz = v[2] - sunS[5];
      const r2 = rx * rx + ry * ry + rz * rz;
      const r = Math.sqrt(r2);
      const u2 = ux * ux + uy * uy + uz * uz;
      const ru = rx * ux + ry * uy + rz * uz;
      const k = gmSun / (c2 * r2 * r);
      const A = 4 * (gmSun / r) - u2;
      ax += k * (A * rx + 4 * ru * ux);
      ay += k * (A * ry + 4 * ru * uy);
      az += k * (A * rz + 4 * ru * uz);
    }
    a[0] = ax;
    a[1] = ay;
    a[2] = az;
  };
}
