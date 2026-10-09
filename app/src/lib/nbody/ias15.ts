// IAS15 (Rein & Spiegel 2015) - 15th-order adaptive Gauss-Radau integrator for a
// single massless particle, y'' = a(t, y, y').  Implemented from the published method:
// the Radau nodes, the Newton-form coefficients and the monomial expansions are all
// generated numerically at module load (nothing is hard-coded from memory).
//
// Features: forward and backward in time (signed dt), predictor-corrector with
// Gauss-Seidel node sweeps, step-size control from the highest-order coefficient,
// compensated (Kahan) summation of the state, and a degree-9 dense-output polynomial
// for every accepted step (StepPoly).

/**
 * Acceleration callback.  x = xRef + dx exactly (as doubles); forces that subtract large nearly
 * equal numbers (particle close to a body 1 AU from the origin) should use (xRef - bodyBase) + dx
 * to avoid the 1e-16 AU rounding noise of x itself.  t is time in days since the caller's origin.
 */
export type AccFn = (
  t: number,
  x: Float64Array,
  v: Float64Array,
  a: Float64Array,
  xRef: Float64Array,
  dx: Float64Array,
  tRef: number,
  dtOff: number,
) => void;
const ZERO3 = new Float64Array(3);

const NODES = 8; // h[0]=0 plus 7 Radau nodes
const NP = 7; // polynomial terms beta_1..beta_7 (a = a0 + sum beta_p tau^p)

/** Gauss-Radau nodes on [0,1) (h[0] = 0). Roots of P7(x)+P8(x) (excluding x=-1) mapped to [0,1]. */
function computeNodes(): number[] {
  const f = (x: number): number => legendre(7, x) + legendre(8, x);
  const roots: number[] = [];
  const N = 4000;
  let xa = -1 + 1e-6;
  let fa = f(xa);
  for (let k = 1; k <= N; k++) {
    const xb = -1 + 1e-6 + (2 - 1e-6) * (k / N);
    const fb = f(xb);
    if (fa === 0) roots.push(xa);
    else if (fa * fb < 0) {
      let lo = xa;
      let hi = xb;
      for (let it = 0; it < 200; it++) {
        const mid = 0.5 * (lo + hi);
        if (f(lo) * f(mid) <= 0) hi = mid;
        else lo = mid;
      }
      roots.push(0.5 * (lo + hi));
    }
    xa = xb;
    fa = fb;
  }
  if (roots.length !== 7) throw new Error(`ias15: expected 7 Radau roots, found ${roots.length}`);
  return [0, ...roots.map((x) => (x + 1) / 2)];
}

function legendre(n: number, x: number): number {
  let p0 = 1;
  let p1 = x;
  if (n === 0) return p0;
  for (let k = 2; k <= n; k++) {
    const p2 = ((2 * k - 1) * x * p1 - (k - 1) * p0) / k;
    p0 = p1;
    p1 = p2;
  }
  return p1;
}

export const RADAU_H: readonly number[] = computeNodes();

// n[k][p]: coefficient of tau^p in N_k(tau) = prod_{i=0}^{k-1} (tau - h_i),  k=1..7, p=1..k
const NK: number[][] = (() => {
  const out: number[][] = [[]];
  let poly: number[] = [1]; // polynomial coefficients c[p] of tau^p
  for (let k = 1; k <= NP; k++) {
    const hi = RADAU_H[k - 1];
    const next = new Array<number>(poly.length + 1).fill(0);
    for (let p = 0; p < poly.length; p++) {
      next[p + 1] += poly[p];
      next[p] -= hi * poly[p];
    }
    poly = next;
    out.push(poly.slice());
  }
  return out;
})();

// Position / velocity weights at the nodes and at the step end.
// x_j = x0 + dt h v0 + dt^2 (a0 h^2/2 + sum_p beta_p h^(p+2) / ((p+1)(p+2)))
// v_j = v0 + dt (a0 h + sum_p beta_p h^(p+1) / (p+1))
const WX: number[][] = [];
const WV: number[][] = [];
for (let j = 0; j < NODES + 1; j++) {
  const h = j < NODES ? RADAU_H[j] : 1;
  const wx = [0];
  const wv = [0];
  for (let p = 1; p <= NP; p++) {
    wx.push(Math.pow(h, p + 2) / ((p + 1) * (p + 2)));
    wv.push(Math.pow(h, p + 1) / (p + 1));
  }
  WX.push(wx);
  WV.push(wv);
}
const END = NODES; // index of the tau=1 weights

const BINOM: number[][] = (() => {
  const b: number[][] = [];
  for (let n = 0; n <= NP; n++) {
    b.push([]);
    for (let k = 0; k <= n; k++) b[n].push(k === 0 || k === n ? 1 : b[n - 1][k - 1] + b[n - 1][k]);
  }
  return b;
})();

/** Dense output of one accepted step. tau in [0,1]. */
export class StepPoly {
  t0 = 0;
  dt = 0;
  readonly x0 = new Float64Array(3);
  readonly v0 = new Float64Array(3);
  readonly a0 = new Float64Array(3);
  /** beta[3*(p-1)+c], p=1..7 */
  readonly beta = new Float64Array(3 * NP);

  get t1(): number {
    return this.t0 + this.dt;
  }

  /** Position (and velocity if vOut given) at tau. */
  eval(tau: number, xOut: Float64Array, vOut?: Float64Array): void {
    const { x0, v0, a0, beta, dt } = this;
    const t2 = tau * tau;
    for (let c = 0; c < 3; c++) {
      let sx = 0;
      let sv = 0;
      let tp = tau; // tau^p
      // Horner-free direct sum; p = 1..7
      for (let p = 1; p <= NP; p++) {
        const b = beta[3 * (p - 1) + c];
        tp *= tau; // tau^(p+1)
        sv += (b * tp) / (p + 1);
        sx += (b * tp * tau) / ((p + 1) * (p + 2));
      }
      xOut[c] = x0[c] + dt * tau * v0[c] + dt * dt * (0.5 * a0[c] * t2 + sx);
      if (vOut) vOut[c] = v0[c] + dt * (a0[c] * tau + sv);
    }
  }
}

export interface IAS15Options {
  /** Dimensionless step-size tolerance (REBOUND default 1e-9). */
  epsilon?: number;
  /** Initial step magnitude (days). */
  dtInit?: number;
  /** Largest |dt| (days). */
  maxDt?: number;
  /** Smallest |dt| (days). */
  minDt?: number;
}

export class IAS15 {
  t: number;
  readonly x: Float64Array;
  readonly v: Float64Array;
  /** acceleration at (t, x, v) */
  readonly a: Float64Array;
  /** suggested next step (signed) */
  dt: number;
  steps = 0;
  rejected = 0;
  evals = 0;

  private readonly acc: AccFn;
  private readonly epsilon: number;
  private readonly maxDt: number;
  private readonly minDt: number;
  private readonly cx = new Float64Array(3);
  private readonly cv = new Float64Array(3);
  private readonly g = new Float64Array(3 * (NP + 1)); // g[3k+c], k=1..7 (k=0 unused)
  private readonly beta = new Float64Array(3 * NP);
  private readonly prevBeta = new Float64Array(3 * NP);
  private prevDt = 0;
  /** previous accepted |dt| and its error estimate, and the running estimate of the force-noise floor */
  private prevAbsDt = 0;
  private prevErr = 0;
  private floorErr = 0;
  private hasPrev = false;
  private readonly aj = new Float64Array(3 * (NODES)); // accelerations at nodes 0..7
  private readonly xt = new Float64Array(3);
  private readonly dxt = new Float64Array(3);
  private readonly vt = new Float64Array(3);
  private readonly at = new Float64Array(3);

  constructor(acc: AccFn, t0: number, x0: ArrayLike<number>, v0: ArrayLike<number>, opts: IAS15Options = {}) {
    this.acc = acc;
    this.t = t0;
    this.x = Float64Array.from(x0);
    this.v = Float64Array.from(v0);
    this.a = new Float64Array(3);
    this.epsilon = opts.epsilon ?? 1e-9;
    this.maxDt = opts.maxDt ?? 1e9;
    this.minDt = opts.minDt ?? 1e-9;
    this.dt = opts.dtInit ?? 0.01;
    this.acc(this.t, this.x, this.v, this.a, this.x, ZERO3, this.t, 0);
    this.evals++;
  }

  /** Reset the state (e.g. to change integration direction) and forget the predictor history. */
  reset(t0: number, x0: ArrayLike<number>, v0: ArrayLike<number>, dt: number): void {
    this.t = t0;
    this.x.set(x0);
    this.v.set(v0);
    this.cx.fill(0);
    this.cv.fill(0);
    this.dt = dt;
    this.hasPrev = false;
    this.prevErr = 0;
    this.floorErr = 0;
    this.acc(this.t, this.x, this.v, this.a, this.x, ZERO3, this.t, 0);
    this.evals++;
  }

  /**
   * Take one accepted step.  |dt| is limited by `dtCap` (> 0, e.g. for close-approach
   * resolution) and the step never crosses `tLimit` (the step is clipped to land on it).
   * The returned StepPoly is reused across calls only if `poly` is passed; otherwise new.
   */
  step(tLimit: number, dtCap = Infinity, poly: StepPoly = new StepPoly()): StepPoly {
    const dir = this.dt >= 0 ? 1 : -1;
    const remaining = (tLimit - this.t) * dir;
    if (!(remaining > 0)) throw new Error('ias15: tLimit not ahead of current time');
    let dtTry = Math.abs(this.dt);
    dtTry = Math.min(dtTry, this.maxDt, dtCap);
    dtTry = Math.max(dtTry, this.minDt);
    const a0 = this.a;
    let a0max = Math.max(Math.abs(a0[0]), Math.abs(a0[1]), Math.abs(a0[2]));
    if (a0max === 0) a0max = 1e-300;

    for (let attempt = 0; attempt < 30; attempt++) {
      let clipped = false;
      let dtAbs = dtTry;
      if (dtAbs >= remaining * (1 - 1e-12)) {
        dtAbs = remaining;
        clipped = true;
      }
      const dt = dir * dtAbs;
      this.predict(dt);
      this.corrector(dt);
      // error estimate: max|b6| / max|a|
      let b6 = 0;
      for (let c = 0; c < 3; c++) b6 = Math.max(b6, Math.abs(this.g[3 * NP + c]));
      const err = b6 / a0max;
      let dtNew: number;
      if (!(err > 0)) dtNew = dtAbs * 4;
      else {
        // Force noise (rounding of the ephemeris / state, 1e-16 AU relative to a close encounter
        // distance) puts a floor under b6/a that does not shrink with dt.  When err stops
        // improving as dt is reduced we raise the effective tolerance instead of crawling at minDt.
        const epsEff = Math.min(Math.max(this.epsilon, 2 * this.floorErr), 1e-6);
        dtNew = dtAbs * Math.pow(epsEff / err, 1 / 7);
      }
      if (!Number.isFinite(dtNew)) dtNew = dtAbs * 0.25;
      if (dtNew / dtAbs < 0.25 && dtAbs > this.minDt * 1.0001) {
        // reject, retry with the smaller step
        this.rejected++;
        dtTry = Math.max(dtNew, this.minDt);
        continue;
      }
      // accept
      if (this.prevErr > 0 && dtAbs < 0.9 * this.prevAbsDt && err > this.epsilon && err > 0.4 * this.prevErr) {
        this.floorErr = Math.max(this.floorErr, err);
      }
      this.floorErr *= 0.98;
      this.prevAbsDt = dtAbs;
      this.prevErr = err;
      this.accept(dt, clipped ? tLimit : null, poly);
      const nextNatural = clipped ? Math.max(dtNew, dtTry * 0.5) : dtNew;
      this.dt = dir * Math.min(Math.max(Math.min(nextNatural, dtAbs * 4), this.minDt), this.maxDt);
      return poly;
    }
    throw new Error('ias15: step rejected too many times');
  }

  /** Initial guess for the node accelerations / g coefficients from the previous step's polynomial. */
  private predict(dt: number): void {
    const { g, beta, aj } = this;
    const a0 = this.a;
    aj[0] = a0[0];
    aj[1] = a0[1];
    aj[2] = a0[2];
    const qq = this.hasPrev ? Math.abs(dt / this.prevDt) : 1;
    if (!this.hasPrev || !(qq < 10 && qq > 0.1)) {
      for (let j = 1; j < NODES; j++) {
        aj[3 * j] = a0[0];
        aj[3 * j + 1] = a0[1];
        aj[3 * j + 2] = a0[2];
      }
    } else {
      const q = dt / this.prevDt;
      // new polynomial coefficients about the old step end: tau_old = 1 + q tau'
      const newBeta = new Float64Array(3 * NP);
      for (let c = 0; c < 3; c++) {
        for (let m = 1; m <= NP; m++) {
          let s = 0;
          for (let p = m; p <= NP; p++) s += this.prevBeta[3 * (p - 1) + c] * BINOM[p][m];
          newBeta[3 * (m - 1) + c] = s * Math.pow(q, m);
        }
      }
      for (let j = 1; j < NODES; j++) {
        const h = RADAU_H[j];
        for (let c = 0; c < 3; c++) {
          let s = a0[c];
          let hp = 1;
          for (let p = 1; p <= NP; p++) {
            hp *= h;
            s += newBeta[3 * (p - 1) + c] * hp;
          }
          aj[3 * j + c] = s;
        }
      }
    }
    // g from node accelerations (nested divided differences), then beta from g
    for (let c = 0; c < 3; c++) {
      for (let k = 1; k <= NP; k++) g[3 * k + c] = this.gFromAccels(k, c);
    }
    beta.fill(0);
    for (let k = 1; k <= NP; k++) {
      for (let c = 0; c < 3; c++) {
        const gk = g[3 * k + c];
        for (let p = 1; p <= k; p++) beta[3 * (p - 1) + c] += gk * NK[k][p];
      }
    }
  }

  /** g_k for component c from aj[0..k] and already-set g_1..g_{k-1}. */
  private gFromAccels(k: number, c: number): number {
    const { aj, g } = this;
    const h = RADAU_H;
    let val = (aj[3 * k + c] - aj[c]) / h[k];
    for (let m = 1; m < k; m++) val = (val - g[3 * m + c]) / (h[k] - h[m]);
    return val;
  }

  /** Predictor-corrector sweeps. Returns true when converged. */
  private corrector(dt: number): boolean {
    const { g, beta, aj, xt, dxt, vt, at } = this;
    const t0 = this.t;
    const x0 = this.x;
    const v0 = this.v;
    const a0 = this.a;
    let lastErr = Infinity;
    for (let iter = 0; iter < 12; iter++) {
      let dg7 = 0;
      let amax = 0;
      for (let j = 1; j < NODES; j++) {
        const wx = WX[j];
        const wv = WV[j];
        const hj = RADAU_H[j];
        for (let c = 0; c < 3; c++) {
          let sx = 0;
          let sv = 0;
          for (let p = 1; p <= NP; p++) {
            const b = beta[3 * (p - 1) + c];
            sx += b * wx[p];
            sv += b * wv[p];
          }
          dxt[c] = dt * hj * v0[c] + dt * dt * (0.5 * a0[c] * hj * hj + sx);
          xt[c] = x0[c] + dxt[c];
          vt[c] = v0[c] + dt * (a0[c] * hj + sv);
        }
        this.acc(t0 + dt * hj, xt, vt, at, x0, dxt, t0, dt * hj);
        this.evals++;
        for (let c = 0; c < 3; c++) {
          aj[3 * j + c] = at[c];
          amax = Math.max(amax, Math.abs(at[c]));
          const gNew = this.gFromAccels(j, c);
          const d = gNew - g[3 * j + c];
          if (j === NP) dg7 = Math.max(dg7, Math.abs(d));
          g[3 * j + c] = gNew;
          if (d !== 0) for (let p = 1; p <= j; p++) beta[3 * (p - 1) + c] += d * NK[j][p];
        }
      }
      const err = amax > 0 ? dg7 / amax : 0;
      if (err < 1e-16) return true;
      if (iter >= 2 && err >= lastErr) return err < 1e-10;
      lastErr = err;
    }
    return false;
  }

  private accept(dt: number, snapT: number | null, poly: StepPoly): void {
    const { beta, x, v, cx, cv } = this;
    const a0 = this.a;
    poly.t0 = this.t;
    poly.dt = dt;
    poly.x0.set(x);
    poly.v0.set(v);
    poly.a0.set(a0);
    poly.beta.set(beta);
    const wx = WX[END];
    const wv = WV[END];
    for (let c = 0; c < 3; c++) {
      let sx = 0;
      let sv = 0;
      for (let p = 1; p <= NP; p++) {
        const b = beta[3 * (p - 1) + c];
        sx += b * wx[p];
        sv += b * wv[p];
      }
      const dx = dt * v[c] + dt * dt * (0.5 * a0[c] + sx);
      const dv = dt * (a0[c] + sv);
      // Kahan-compensated accumulation
      let y = dx - cx[c];
      let tt = x[c] + y;
      cx[c] = tt - x[c] - y;
      x[c] = tt;
      y = dv - cv[c];
      tt = v[c] + y;
      cv[c] = tt - v[c] - y;
      v[c] = tt;
    }
    this.t = snapT !== null ? snapT : this.t + dt;
    this.prevBeta.set(beta);
    this.prevDt = dt;
    this.hasPrev = true;
    this.steps++;
    this.acc(this.t, this.x, this.v, this.a, this.x, ZERO3, this.t, 0);
    this.evals++;
  }
}
