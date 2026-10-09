// Generates public/data/asteroids.mock.json: ~40k random but plausible NEO
// entries following the asteroids.json contract, used until the real dataset
// exists. Run: npm run gen:mock [count]
import { writeFileSync, mkdirSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";

const count = Number(process.argv[2]) || 40000;
let seed = 0x9e3779b9;
function rand() {
  seed |= 0; seed = (seed + 0x6d2b79f5) | 0;
  let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
  t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
  return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
}
const uni = (a, b) => a + (b - a) * rand();
const gauss = () => Math.sqrt(-2 * Math.log(1 - rand())) * Math.cos(2 * Math.PI * rand());
const r = (x, d) => Number(x.toFixed(d));

const EPOCH = 2461000.5; // 2025-Nov-21
const letters = "ABCDEFGHJKLMNOPQRSTUVWXY";

function elements(category) {
  let q, Q;
  if (category === "impact") {
    // orbit crosses 1 AU: perihelion inside, aphelion outside
    q = uni(0.35, 0.995);
    Q = uni(1.005, 4.2);
  } else if (category === "approaching") {
    if (rand() < 0.55) { q = uni(0.5, 1.05); Q = uni(1.0, 3.5); } else { q = uni(1.0, 1.25); Q = uni(1.3, 3.2); }
  } else {
    q = uni(1.0, 1.3);
    Q = q + Math.abs(gauss()) * 1.6 + 0.2;
  }
  const a = (q + Q) / 2;
  const e = Math.min(0.97, (Q - q) / (Q + q));
  const i = Math.min(60, -Math.log(1 - rand()) * 9);
  return { a, e, i };
}

const rows = [];
for (let n = 0; n < count; n++) {
  const u = rand();
  const category = u < 0.012 ? "impact" : u < 0.34 ? "approaching" : "receding";
  const { a, e, i } = elements(category);
  const diameter = Math.min(30, Math.max(0.001, Math.exp(Math.log(0.12) + 1.25 * gauss())));
  let min_dist, pha;
  if (category === "impact") { min_dist = uni(0.00005, 0.004); pha = rand() < 0.75; }
  else if (category === "approaching") { min_dist = uni(0.004, 0.08); pha = diameter > 0.14 && min_dist < 0.05 && rand() < 0.5; }
  else { min_dist = uni(0.05, 0.9); pha = false; }
  const year = 1990 + Math.floor(rand() * 36);
  const des = `${year} ${letters[Math.floor(rand() * 24)]}${letters[Math.floor(rand() * 24)]}${rand() < 0.5 ? Math.floor(rand() * 300) : ""}`;
  const date = new Date(Date.UTC(2026, 0, 1) + Math.floor(rand() * 3650) * 86400000).toISOString().slice(0, 10);
  rows.push({
    id: String(20000000 + n * 7 + Math.floor(rand() * 7)),
    name: `${des} (mock)`,
    diameter_km: r(diameter, 4),
    pha,
    category,
    min_dist_au: r(min_dist, 6),
    min_date: date,
    a: r(a, 5),
    e: r(e, 5),
    i: r(i, 4),
    om: r(uni(0, 360), 4),
    w: r(uni(0, 360), 4),
    ma: r(uni(0, 360), 4),
    epoch: EPOCH,
  });
}

const out = join(dirname(fileURLToPath(import.meta.url)), "..", "public", "data", "asteroids.mock.json");
mkdirSync(dirname(out), { recursive: true });
writeFileSync(out, JSON.stringify(rows));
const tally = rows.reduce((m, x) => ((m[x.category] = (m[x.category] || 0) + 1), m), {});
console.log(`wrote ${rows.length} mock asteroids to ${out}`, tally);
