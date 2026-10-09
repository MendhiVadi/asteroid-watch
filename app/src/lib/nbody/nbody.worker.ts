// Web Worker wrapper.  Usage (Vite):
//   const w = new Worker(new URL('./nbody.worker.ts', import.meta.url), { type: 'module' });
// or use NBodyClient from client.ts.
//
// Messages in :
//   { type:'init', baseUrl?: string }                          load ephemeris (default '/data/ephemeris/')
//   { type:'run', id, asteroid, jdStart, opts?, includeMoon? }  stream a simulation
//   { type:'cancel', id }
// Messages out:
//   { type:'ready', coverage:[jdStart, jdEnd] }
//   { type:'chunk', id, jd, earth, helio, moon?, dEarth, dMoon, events, jdReached, progress }   (typed arrays transferred)
//   { type:'done', id, summary }
//   { type:'error', id?, code?, message }

import { Ephemeris, loadEphemeris } from './ephemeris.ts';
import { NBodyRangeError, simulateStream } from './sim.ts';
import type { SimEvent, SimOptions, SimSummary } from './sim.ts';
import type { AsteroidElements, AsteroidState } from './elements.ts';
import { samplesToArrays } from './trajectory.ts';

export type WorkerIn =
  | { type: 'init'; baseUrl?: string }
  | {
      type: 'run';
      id: number;
      asteroid: AsteroidElements | AsteroidState;
      jdStart: number;
      opts?: Omit<SimOptions, 'ephemeris'>;
      includeMoon?: boolean;
    }
  | { type: 'cancel'; id: number };

export type WorkerOut =
  | { type: 'ready'; coverage: [number, number] }
  | {
      type: 'chunk';
      id: number;
      jd: Float64Array;
      earth: Float32Array;
      helio: Float32Array;
      moon?: Float32Array;
      dEarth: Float32Array;
      dMoon: Float32Array;
      events: SimEvent[];
      jdReached: number;
      progress: number;
    }
  | { type: 'done'; id: number; summary: SimSummary }
  | { type: 'error'; id?: number; code?: string; message: string };

interface WorkerScope {
  onmessage: ((ev: MessageEvent<WorkerIn>) => void) | null;
  postMessage(msg: WorkerOut, transfer?: Transferable[]): void;
}
const ctx = self as unknown as WorkerScope;

let eph: Ephemeris | null = null;
let ephPromise: Promise<Ephemeris> | null = null;
const cancelled = new Set<number>();

const yieldToEventLoop = (): Promise<void> => new Promise((r) => setTimeout(r, 0));

async function run(m: Extract<WorkerIn, { type: 'run' }>): Promise<void> {
  try {
    const e = eph ?? (await (ephPromise ?? (ephPromise = loadEphemeris())));
    eph = e;
    const gen = simulateStream(m.asteroid, m.jdStart, { ...(m.opts ?? {}), ephemeris: e });
    for (;;) {
      if (cancelled.has(m.id)) {
        cancelled.delete(m.id);
        return;
      }
      const r = gen.next();
      if (r.done) {
        ctx.postMessage({ type: 'done', id: m.id, summary: r.value });
        return;
      }
      const a = samplesToArrays(r.value.samples, e, !!m.includeMoon);
      const transfer: Transferable[] = [a.jd.buffer, a.earth.buffer, a.helio.buffer, a.dEarth.buffer, a.dMoon.buffer];
      if (a.moon) transfer.push(a.moon.buffer);
      ctx.postMessage(
        {
          type: 'chunk',
          id: m.id,
          ...a,
          events: r.value.events,
          jdReached: r.value.jdReached,
          progress: r.value.progress,
        },
        transfer,
      );
      await yieldToEventLoop(); // lets 'cancel' messages in
    }
  } catch (err) {
    const msg = err instanceof Error ? err.message : String(err);
    ctx.postMessage({ type: 'error', id: m.id, code: err instanceof NBodyRangeError ? err.code : undefined, message: msg });
  }
}

ctx.onmessage = (ev) => {
  const m = ev.data;
  if (m.type === 'init') {
    ephPromise = loadEphemeris(m.baseUrl);
    ephPromise
      .then((e) => {
        eph = e;
        ctx.postMessage({ type: 'ready', coverage: [e.jdStart, e.jdEnd] });
      })
      .catch((err: unknown) => ctx.postMessage({ type: 'error', message: `ephemeris load failed: ${String(err)}` }));
  } else if (m.type === 'run') {
    cancelled.delete(m.id);
    void run(m);
  } else if (m.type === 'cancel') {
    cancelled.add(m.id);
  }
};
