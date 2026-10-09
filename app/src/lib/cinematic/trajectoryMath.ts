import { AU_KM } from "../kepler";
import type { Trajectory } from "./state";

/** Binary search for the sample interval containing t. */
export function findIndex(tr: Trajectory, t: number, hint = 0): number {
  const a = tr.t;
  let lo = 0;
  let hi = tr.n - 1;
  if (t <= a[0]) return 0;
  if (t >= a[hi]) return Math.max(0, hi - 1);
  if (hint >= 0 && hint < hi && a[hint] <= t && t < a[hint + 1]) return hint;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (a[mid] <= t) lo = mid;
    else hi = mid;
  }
  return lo;
}

/** Interpolates the asteroid (and, if sampled, Moon) Earth-relative position; returns speed in km/s. */
export function sampleAt(
  tr: Trajectory,
  t: number,
  hint: number,
  rel: Float64Array,
  moon: Float64Array | null,
): { idx: number; speedKms: number } {
  const i = findIndex(tr, t, hint);
  const t0 = tr.t[i];
  const t1 = tr.t[i + 1] ?? t0 + 1e-6;
  const u = Math.min(1, Math.max(0, (t - t0) / (t1 - t0)));
  const j = i * 3;
  const k = Math.min(i + 1, tr.n - 1) * 3;
  for (let c = 0; c < 3; c++) {
    rel[c] = tr.rel[j + c] + (tr.rel[k + c] - tr.rel[j + c]) * u;
    if (moon && tr.moon) moon[c] = tr.moon[j + c] + (tr.moon[k + c] - tr.moon[j + c]) * u;
  }
  const dt = Math.max(t1 - t0, 1e-9);
  const vx = tr.rel[k] - tr.rel[j];
  const vy = tr.rel[k + 1] - tr.rel[j + 1];
  const vz = tr.rel[k + 2] - tr.rel[j + 2];
  return { idx: i, speedKms: (Math.hypot(vx, vy, vz) / dt) * (AU_KM / 86400) };
}
