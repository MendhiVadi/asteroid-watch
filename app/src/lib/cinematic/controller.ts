import { cinemaOffsetToScene } from "../kepler";
import { clock, useStore } from "../store";
import { jdFromMs } from "../kepler";
import { jdTdbFromUnixMs, unixMsFromJdTdb } from "../nbody/constants";
import type { NBodyError } from "../nbody/client";
import type { TrajectoryArrays } from "../nbody/trajectory";
import type { SimEvent, SimSummary } from "../nbody/sim";
import { coverage, getClient, loadMainEphemeris, ready } from "./engine";
import { toWide } from "./ephemeris";
import { cine, frame, type Trajectory, type ViewMode } from "./state";

let currentJob: { id: number } | null = null;
let runToken = 0;

/** JD UTC (app clock) -> JD TDB (engine time). */
export const toTdb = (jdUtc: number) => jdTdbFromUnixMs((jdUtc - 2440587.5) * 86400000);
/** JD TDB -> JD UTC (app clock). */
export const fromTdb = (jdTdb: number) => jdFromMs(unixMsFromJdTdb(jdTdb));

function emptyTrajectory(tStart: number): Trajectory {
  return {
    n: 0,
    t: new Float64Array(4096),
    rel: new Float32Array(4096 * 3),
    helio: new Float32Array(4096 * 3),
    moon: new Float32Array(4096 * 3),
    status: "window_end",
    tStart,
    reached: tStart,
    complete: false,
    tStop: tStart,
    truncated: false,
    closestAu: Infinity,
    closestJd: tStart,
    events: [],
    engine: "IAS15 N-body (Sun, planets, Earth, Moon)",
  };
}

function grow<T extends Float32Array | Float64Array>(arr: T, need: number): T {
  if (arr.length >= need) return arr;
  let cap = arr.length;
  while (cap < need) cap *= 2;
  const next = new (arr.constructor as new (n: number) => T)(cap);
  next.set(arr);
  return next;
}

function appendChunk(tr: Trajectory, a: TrajectoryArrays, events: SimEvent[], jdReached: number) {
  const m = a.jd.length;
  const need = tr.n + m;
  tr.t = grow(tr.t, need);
  tr.rel = grow(tr.rel, need * 3);
  tr.helio = grow(tr.helio, need * 3);
  if (tr.moon) tr.moon = grow(tr.moon, need * 3);
  tr.t.set(a.jd, tr.n);
  tr.rel.set(a.earth, tr.n * 3);
  tr.helio.set(a.helio, tr.n * 3);
  if (tr.moon && a.moon) tr.moon.set(a.moon, tr.n * 3);
  for (let i = 0; i < m; i++) {
    if (a.dEarth[i] < tr.closestAu) {
      tr.closestAu = a.dEarth[i];
      tr.closestJd = a.jd[i];
    }
  }
  tr.n = need;
  tr.reached = jdReached;
  tr.events.push(...events);
  extendPaths(tr, tr.n - m);
}

let localPath = new Float32Array(0);
let widePath = new Float32Array(0);

/** Maintains the scene-space polylines for samples [from, n). */
function extendPaths(tr: Trajectory, from: number) {
  localPath = grow(localPath.length ? localPath : new Float32Array(4096 * 3), tr.n * 3);
  widePath = grow(widePath.length ? widePath : new Float32Array(4096 * 3), tr.n * 3);
  for (let i = from; i < tr.n; i++) {
    cinemaOffsetToScene(tr.rel[i * 3], tr.rel[i * 3 + 1], tr.rel[i * 3 + 2], localPath, i * 3);
    toWide(tr.helio[i * 3], tr.helio[i * 3 + 1], tr.helio[i * 3 + 2], widePath, i * 3);
  }
  cine.localPath = localPath;
  cine.widePath = widePath;
  cine.pathVersion++;
}

function describe(err: unknown): { message: string; code?: string } {
  const e = err as Partial<NBodyError> & { coverageJd?: [number, number] };
  if (e?.code === "OUT_OF_COVERAGE") {
    return {
      code: e.code,
      message: "That date is outside the ephemeris coverage (2020-01-01 to 2036-12-31). Pick a date inside it.",
    };
  }
  return { code: e?.code, message: err instanceof Error ? err.message : String(err) };
}

function cancelJob() {
  if (currentJob) {
    try {
      getClient().cancel(currentJob.id);
    } catch {
      /* worker already gone */
    }
    currentJob = null;
  }
}

/**
 * Starts the N-body run for the asteroid in a worker and switches to cinematic mode as soon as the first
 * chunk arrives (playback starts while the rest is still being integrated).
 */
export async function startCinema(index: number): Promise<void> {
  const store = useStore.getState();
  const data = store.data;
  if (!data) return;
  cancelJob();
  const token = ++runToken;
  if (cine.active) {
    cine.active = false;
    cine.traj = null;
  }
  useStore.setState({ cinema: "loading", cinemaProgress: 0, cinemaError: null, cinemaEnd: null, cinemaNote: null });

  let cov: [number, number];
  try {
    [cov] = await Promise.all([ready(), loadMainEphemeris()]);
  } catch (err) {
    if (token === runToken) useStore.setState({ cinema: "error", cinemaError: describe(err).message });
    return;
  }
  if (token !== runToken) return;

  // Start at the app's current simulation date, clamped into the ephemeris coverage.
  const wanted = toTdb(clock.jd);
  const tStart = Math.min(Math.max(wanted, cov[0] + 0.5), cov[1] - 30);
  const note =
    Math.abs(tStart - wanted) > 1e-6
      ? "Start date was moved into the ephemeris coverage (2020-01-01 to 2036-12-31)."
      : null;

  const a = data.list[index];
  const tr = emptyTrajectory(tStart);
  localPath = new Float32Array(0);
  widePath = new Float32Array(0);
  let started = false;
  const begin = () => {
    started = true;
    cine.traj = tr;
    cine.asteroid = index;
    cine.t = tr.tStart;
    cine.playing = true;
    cine.ended = false;
    cine.needsCameraInit = true;
    cine.active = true;
    useStore.setState({ cinema: "ready", cinemaNote: note });
  };

  const job = getClient().run(
    { a: a.a, e: a.e, i: a.i, om: a.om, w: a.w, ma: a.ma, epoch: a.epoch, id: a.id, name: a.name },
    tStart,
    { direction: 1, chunkSamples: 1500 },
    {
      includeMoon: true,
      onChunk: (arrays, info) => {
        if (token !== runToken) return;
        appendChunk(tr, arrays, info.events, info.jdReached);
        useStore.setState({ cinemaProgress: info.progress });
        if (!started && tr.n >= 2) begin();
      },
    },
  );
  currentJob = job;

  try {
    const done = await job;
    if (token !== runToken) return;
    finish(tr, done.summary);
    if (!started) begin();
  } catch (err) {
    if (token !== runToken) return;
    const d = describe(err);
    if (d.code === "CANCELLED" || d.code === "DISPOSED") return;
    if (started) {
      // Keep what was integrated, but make clear the run did not finish.
      tr.complete = true;
      tr.truncated = true;
      tr.tStop = tr.reached;
      useStore.setState({ cinemaNote: `Simulation stopped early: ${d.message}` });
    } else {
      useStore.setState({ cinema: "error", cinemaError: d.message });
    }
  } finally {
    if (currentJob === job) currentJob = null;
  }
}

function finish(tr: Trajectory, summary: SimSummary) {
  tr.status = summary.status;
  tr.tStop = summary.stopJd;
  tr.reached = summary.stopJd;
  tr.truncated = summary.truncated;
  tr.escapeJd = summary.escapeJd;
  tr.closestAu = summary.minEarthDistAu;
  tr.closestJd = summary.minEarthJd;
  tr.events = summary.events;
  tr.complete = true;
  useStore.setState({ cinemaProgress: 1 });
}

export function exitCinema(): void {
  runToken++;
  cancelJob();
  if (cine.active) clock.jd = Math.min(Math.max(fromTdb(cine.t), fromTdb(coverage[0])), fromTdb(coverage[1]));
  cine.active = false;
  cine.traj = null;
  useStore.setState({ cinema: "off", cinemaEnd: null, cinemaError: null, cinemaNote: null });
  useStore.getState().view("home");
}

export function setPlaying(playing: boolean): void {
  if (!cine.traj) return;
  if (playing && cine.ended) replay();
  cine.playing = playing;
}

export function replay(): void {
  const tr = cine.traj;
  if (!tr) return;
  cine.t = tr.tStart;
  cine.ended = false;
  cine.playing = true;
  useStore.setState({ cinemaEnd: null });
}

export function scrub(t: number): void {
  const tr = cine.traj;
  if (!tr) return;
  const end = tr.complete ? tr.tStop : tr.reached;
  cine.t = Math.min(Math.max(t, tr.tStart), end);
  if (cine.ended && cine.t < tr.tStop) {
    cine.ended = false;
    useStore.setState({ cinemaEnd: null });
  }
}

export function setView(view: ViewMode): void {
  cine.view = view;
}

// Dev-only handle for poking at the running app from the console.
if (import.meta.env.DEV) {
  (window as unknown as { __cine: unknown }).__cine = { cine, frame, startCinema, exitCinema, toTdb, fromTdb };
}
