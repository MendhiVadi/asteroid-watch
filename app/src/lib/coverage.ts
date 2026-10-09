import { useStore } from "./store";

// Simulation range = ephemeris coverage, read from the manifest (and refined from the worker's `ready`).
// Nothing else in the app hardcodes an end year. All values are JD; UTC vs TDB differs by ~69 s, which is
// irrelevant for range checks.

const MANIFEST_URL = "/data/ephemeris/manifest.json";

interface ManifestLike {
  jdStart?: number;
  jdEnd?: number;
  main: { jd0: number; stepDays: number; count: number };
  moon: { jd0: number; stepDays: number; count: number };
}

export const coverage = {
  /** First JD the tables cover. */
  start: 2458849.5,
  /** Last JD the tables cover. */
  end: 2488433.5,
  loaded: false,
};

/** UTC midnight of the current day, as a JD: the earliest date the dashboard clock can show. */
export const TODAY_JD = Math.floor(Date.now() / 86400000) + 2440587.5;

export function setCoverage(start: number, end: number): void {
  if (coverage.loaded && coverage.start === start && coverage.end === end) return;
  coverage.start = start;
  coverage.end = end;
  coverage.loaded = true;
  useStore.setState({ coverageEnd: end });
}

export async function loadCoverage(): Promise<void> {
  try {
    const res = await fetch(MANIFEST_URL, { cache: "no-cache" });
    if (!res.ok) return;
    const m = (await res.json()) as ManifestLike;
    const start = m.jdStart ?? Math.max(m.main.jd0, m.moon.jd0);
    const end = m.jdEnd ?? Math.min(m.main.jd0 + (m.main.count - 1) * m.main.stepDays, m.moon.jd0 + (m.moon.count - 1) * m.moon.stepDays);
    setCoverage(start, end);
  } catch {
    /* keep defaults; the worker's `ready` will correct them */
  }
}

/** Lowest / highest JD the dashboard clock and date picker accept: [today .. end of coverage]. */
export function clockRange(): [number, number] {
  return [Math.max(TODAY_JD, coverage.start), coverage.end];
}

export function clampClock(jd: number): number {
  const [lo, hi] = clockRange();
  return Math.min(Math.max(jd, lo), hi);
}

/** YYYY-MM-DD (UTC) for a JD. */
export function jdToDateString(jd: number): string {
  return new Date((jd - 2440587.5) * 86400000).toISOString().slice(0, 10);
}
