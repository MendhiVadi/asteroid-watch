import { buildElementTable } from "./kepler";

export type Category = "impact" | "approaching" | "receding";
export const CATEGORIES: Category[] = ["impact", "approaching", "receding"];

export interface Asteroid {
  id: string;
  name: string;
  diameter_km: number;
  pha: boolean;
  category: Category;
  min_dist_au: number;
  min_date: string;
  a: number;
  e: number;
  i: number;
  om: number;
  w: number;
  ma: number;
  epoch: number;
  /** Set by the data pipeline for objects on a strict collision course. */
  strict_collision?: boolean;
  /** Change of Earth distance over the analysis window, AU (negative = closing). */
  trend_au?: number;
}

export interface AsteroidData {
  list: Asteroid[];
  count: number;
  /** Packed Keplerian table, see kepler.ts. */
  table: Float64Array;
  /** Category index per object (0 impact, 1 approaching, 2 receding). */
  cat: Uint8Array;
  /** Marker scale per object in scene units (log of diameter). */
  size: Float32Array;
  /** Lower-case "name id" strings for search. */
  haystack: string[];
  counts: Record<Category, number>;
  stats: {
    pha: number;
    closest: Asteroid | null;
    largest: Asteroid | null;
  };
  source: "real" | "mock";
}

const CAT_INDEX: Record<Category, number> = { impact: 0, approaching: 1, receding: 2 };

async function fetchJson(url: string): Promise<unknown> {
  const res = await fetch(url, { cache: "no-cache" });
  if (!res.ok) throw new Error(`${url}: HTTP ${res.status}`);
  const text = await res.text();
  // Vite's dev server answers unknown paths with index.html (HTTP 200).
  if (text.trimStart().startsWith("<")) throw new Error(`${url}: not JSON`);
  return JSON.parse(text);
}

function markerSize(diameterKm: number): number {
  const d = Number.isFinite(diameterKm) && diameterKm > 0 ? diameterKm : 0.05;
  return 0.034 + 0.058 * Math.log10(1 + d * 100);
}

function normalise(raw: unknown): Asteroid[] {
  const rows = Array.isArray(raw) ? raw : (raw as { asteroids?: unknown })?.asteroids;
  if (!Array.isArray(rows)) throw new Error("asteroid data is not an array");
  const out: Asteroid[] = [];
  for (const r of rows as Record<string, unknown>[]) {
    const category = r.category as Category;
    const a = Number(r.a), e = Number(r.e), i = Number(r.i);
    const om = Number(r.om), w = Number(r.w), ma = Number(r.ma), epoch = Number(r.epoch);
    if (!(category in CAT_INDEX)) continue;
    if (![a, e, i, om, w, ma, epoch].every(Number.isFinite) || a <= 0 || e < 0 || e >= 1) continue;
    out.push({
      id: String(r.id ?? ""),
      name: String(r.name ?? r.id ?? "Unnamed"),
      diameter_km: Number(r.diameter_km) || 0,
      pha: Boolean(r.pha),
      category,
      min_dist_au: Number(r.min_dist_au) || 0,
      min_date: String(r.min_date ?? ""),
      a, e, i, om, w, ma, epoch,
      strict_collision: r.strict_collision === undefined ? undefined : Boolean(r.strict_collision),
      trend_au: Number.isFinite(Number(r.trend_au)) ? Number(r.trend_au) : undefined,
    });
  }
  if (out.length === 0) throw new Error("no valid asteroid rows");
  return out;
}

function build(list: Asteroid[], source: "real" | "mock"): AsteroidData {
  const count = list.length;
  const cat = new Uint8Array(count);
  const size = new Float32Array(count);
  const haystack = new Array<string>(count);
  const counts: Record<Category, number> = { impact: 0, approaching: 0, receding: 0 };
  let pha = 0;
  let closest: Asteroid | null = null;
  let largest: Asteroid | null = null;
  for (let k = 0; k < count; k++) {
    const a = list[k];
    cat[k] = CAT_INDEX[a.category];
    size[k] = markerSize(a.diameter_km);
    haystack[k] = `${a.name} ${a.id}`.toLowerCase();
    counts[a.category]++;
    if (a.pha) pha++;
    if (a.min_dist_au > 0 && (!closest || a.min_dist_au < closest.min_dist_au)) closest = a;
    if (!largest || a.diameter_km > largest.diameter_km) largest = a;
  }
  return { list, count, table: buildElementTable(list), cat, size, haystack, counts, stats: { pha, closest, largest }, source };
}

/** Prefers /data/asteroids.json, falls back to the generated mock set. */
let pending: Promise<AsteroidData> | null = null;
export function loadAsteroids(): Promise<AsteroidData> {
  // Memoised so React StrictMode's double effect does not parse ~12 MB twice.
  return (pending ??= load());
}

async function load(): Promise<AsteroidData> {
  try {
    return build(normalise(await fetchJson("/data/asteroids.json")), "real");
  } catch (realError) {
    console.info("[asteroids] real dataset unavailable, using mock:", (realError as Error).message);
    return build(normalise(await fetchJson("/data/asteroids.mock.json")), "mock");
  }
}
