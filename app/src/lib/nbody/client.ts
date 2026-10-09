// Main-thread helper around nbody.worker.ts.
//
//   const nb = new NBodyClient({ baseUrl: '/data/ephemeris/' });
//   await nb.ready;                                   // [coverageStartJd, coverageEndJd]
//   const job = nb.run(elements, jdStart, { direction: 1 }, { includeMoon: true,
//                  onChunk: (arrays, info) => ... }); // progressive playback while it streams
//   const traj = await job;  nb.cancel(job.id);  nb.dispose();

import type { AsteroidElements, AsteroidState } from './elements.ts';
import type { SimEvent, SimOptions } from './sim.ts';
import { concatArrays, makeTrajectory } from './trajectory.ts';
import type { Trajectory, TrajectoryArrays } from './trajectory.ts';
import type { WorkerIn, WorkerOut } from './nbody.worker.ts';

export interface RunOptions {
  includeMoon?: boolean;
  onChunk?: (arrays: TrajectoryArrays, info: { events: SimEvent[]; jdReached: number; progress: number }) => void;
}

export class NBodyError extends Error {
  readonly code?: string;
  constructor(message: string, code?: string) {
    super(message);
    this.name = 'NBodyError';
    this.code = code;
  }
}

interface Pending {
  parts: TrajectoryArrays[];
  opts: RunOptions;
  resolve: (t: Trajectory) => void;
  reject: (e: Error) => void;
}

export class NBodyClient {
  private worker: Worker;
  private nextId = 1;
  private pending = new Map<number, Pending>();
  /** resolves with the ephemeris coverage [jdStart, jdEnd] once the tables are loaded */
  readonly ready: Promise<[number, number]>;

  constructor(opts: { baseUrl?: string; worker?: Worker } = {}) {
    this.worker = opts.worker ?? new Worker(new URL('./nbody.worker.ts', import.meta.url), { type: 'module' });
    this.ready = new Promise((resolve, reject) => {
      this.worker.onmessage = (ev: MessageEvent<WorkerOut>): void => {
        const m = ev.data;
        if (m.type === 'ready') resolve(m.coverage);
        else if (m.type === 'error' && m.id === undefined) reject(new NBodyError(m.message, m.code));
        this.handle(m);
      };
    });
    this.post({ type: 'init', baseUrl: opts.baseUrl });
  }

  private post(m: WorkerIn): void {
    this.worker.postMessage(m);
  }

  private handle(m: WorkerOut): void {
    if (m.type === 'chunk') {
      const p = this.pending.get(m.id);
      if (!p) return;
      const arrays: TrajectoryArrays = {
        jd: m.jd,
        earth: m.earth,
        helio: m.helio,
        moon: m.moon,
        dEarth: m.dEarth,
        dMoon: m.dMoon,
      };
      p.parts.push(arrays);
      p.opts.onChunk?.(arrays, { events: m.events, jdReached: m.jdReached, progress: m.progress });
    } else if (m.type === 'done') {
      const p = this.pending.get(m.id);
      if (!p) return;
      this.pending.delete(m.id);
      p.resolve(makeTrajectory(concatArrays(p.parts), m.summary));
    } else if (m.type === 'error' && m.id !== undefined) {
      const p = this.pending.get(m.id);
      if (!p) return;
      this.pending.delete(m.id);
      p.reject(new NBodyError(m.message, m.code));
    }
  }

  /** Start a simulation; the returned promise resolves with the full Trajectory (also has `.id`). */
  run(
    asteroid: AsteroidElements | AsteroidState,
    jdStart: number,
    simOpts: Omit<SimOptions, 'ephemeris'> = {},
    opts: RunOptions = {},
  ): Promise<Trajectory> & { id: number } {
    const id = this.nextId++;
    const p = new Promise<Trajectory>((resolve, reject) => {
      this.pending.set(id, { parts: [], opts, resolve, reject });
    });
    this.post({ type: 'run', id, asteroid, jdStart, opts: simOpts, includeMoon: opts.includeMoon });
    return Object.assign(p, { id });
  }

  cancel(id: number): void {
    const p = this.pending.get(id);
    this.post({ type: 'cancel', id });
    if (p) {
      this.pending.delete(id);
      p.reject(new NBodyError('cancelled', 'CANCELLED'));
    }
  }

  dispose(): void {
    this.worker.terminate();
    for (const p of this.pending.values()) p.reject(new NBodyError('disposed', 'DISPOSED'));
    this.pending.clear();
  }
}
