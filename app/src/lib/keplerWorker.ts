import { computeScenePositions } from "./kepler";

interface InitMessage {
  type: "init";
  table: Float64Array;
  count: number;
}
interface StepMessage {
  type: "step";
  jd: number;
  buf: Float32Array;
}

const ctx = self as unknown as {
  onmessage: ((ev: MessageEvent<InitMessage | StepMessage>) => void) | null;
  postMessage: (message: unknown, transfer?: Transferable[]) => void;
};

let table: Float64Array = new Float64Array(0);
let count = 0;

ctx.onmessage = (ev) => {
  const msg = ev.data;
  if (msg.type === "init") {
    table = msg.table;
    count = msg.count;
    ctx.postMessage({ type: "ready" });
  } else {
    computeScenePositions(table, count, msg.jd, msg.buf);
    ctx.postMessage({ type: "positions", jd: msg.jd, buf: msg.buf }, [msg.buf.buffer]);
  }
};
