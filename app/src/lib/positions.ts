import { computeScenePositions } from "./kepler";

interface WorkerReply {
  type: "ready" | "positions";
  jd?: number;
  buf?: Float32Array;
}

// Double-buffered position propagation. The Kepler solve for ~40k objects runs
// in a worker; the render loop keeps reading the last completed buffer, so the
// main thread never blocks on propagation. Falls back to inline compute if the
// worker cannot be created.
export class PositionEngine {
  readonly count: number;
  private table: Float64Array;
  private worker: Worker | null = null;
  private current: Float32Array;
  private spare: Float32Array | null;
  private pending = false;
  private workerReady = false;
  private lastJd = Number.NaN;
  /** True once a worker has permanently failed; computation then stays inline. */
  private inline = false;
  ready = false;

  // The worker is created lazily by start() (called from an effect and, as a safety
  // net, from update()). dispose() only stops the worker, so a disposed engine that is
  // kept alive by a memo (React StrictMode mount/unmount/mount) restarts cleanly.
  constructor(table: Float64Array, count: number) {
    this.table = table;
    this.count = count;
    this.current = new Float32Array(count * 3);
    this.spare = new Float32Array(count * 3);
  }

  get positions(): Float32Array {
    return this.current;
  }

  /** Idempotent: spin up the worker unless one is running or inline mode is forced. */
  start(): void {
    if (this.worker || this.inline) return;
    this.workerReady = false;
    this.pending = false;
    this.lastJd = Number.NaN;
    this.spare ??= new Float32Array(this.count * 3);
    try {
      const w = new Worker(new URL("./keplerWorker.ts", import.meta.url), { type: "module" });
      this.worker = w;
      w.onmessage = (ev: MessageEvent<WorkerReply>) => {
        if (this.worker === w) this.onMessage(ev.data);
      };
      w.onerror = () => {
        if (this.worker === w) this.fallbackInline();
      };
      w.postMessage({ type: "init", table: this.table.slice(), count: this.count });
    } catch {
      this.worker = null;
      this.inline = true;
    }
  }

  private onMessage(msg: WorkerReply) {
    if (msg.type === "ready") {
      this.workerReady = true;
    } else if (msg.buf) {
      this.spare = this.current;
      this.current = msg.buf;
      this.pending = false;
      this.ready = true;
    }
  }

  private fallbackInline() {
    this.worker?.terminate();
    this.worker = null;
    this.inline = true;
    this.pending = false;
    this.lastJd = Number.NaN;
    this.spare ??= new Float32Array(this.count * 3);
  }

  /** Request positions for jd (cheap no-op if unchanged or a request is in flight). */
  update(jd: number): void {
    if (!this.worker && !this.inline) this.start();
    if (this.worker) {
      if (!this.workerReady || this.pending || !this.spare) return;
      if (jd === this.lastJd && this.ready) return;
      this.lastJd = jd;
      this.pending = true;
      const buf = this.spare;
      this.spare = null;
      this.worker.postMessage({ type: "step", jd, buf }, [buf.buffer]);
    } else {
      if (jd === this.lastJd && this.ready) return;
      this.lastJd = jd;
      computeScenePositions(this.table, this.count, jd, this.current);
      this.ready = true;
    }
  }

  /** Stops the worker. The engine can be started again; an in-flight buffer is replaced. */
  dispose() {
    this.worker?.terminate();
    this.worker = null;
    this.workerReady = false;
    this.pending = false;
    this.lastJd = Number.NaN;
    this.spare ??= new Float32Array(this.count * 3);
  }
}
