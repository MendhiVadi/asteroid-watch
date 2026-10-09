import { create } from "zustand";
import { jdFromMs } from "./kepler";
import type { AsteroidData, Category } from "./data";

/** Simulation speeds in days of simulated time per real second. 0 = paused. */
export const SPEEDS = [-100, -30, -10, -3, -1, 0, 0.25, 1, 3, 10, 30, 100];
export const PAUSED_INDEX = 5;
export const DEFAULT_SPEED_INDEX = 7;

/** Mutable simulation clock - read every frame, so it lives outside React state. */
export const clock = { jd: jdFromMs(Date.now()) };

/** Per-frame shared values written by the scene and read by overlays. */
export const live = {
  positions: null as Float32Array | null,
  visible: new Uint32Array(0),
  visibleCount: 0,
};

export type ViewCommand = { kind: "earth" | "wide" | "home" | "none"; n: number };
export type CinemaPhase = "off" | "loading" | "ready" | "error";

interface State {
  data: AsteroidData | null;
  loadState: "loading" | "ready" | "error";
  error: string | null;
  show: Record<Category, boolean>;
  selected: number;
  hover: number;
  follow: boolean;
  /** Increments whenever a selection is made, so the camera can fly to it. */
  flyTick: number;
  speedIdx: number;
  lastSpeedIdx: number;
  reducedMotion: boolean;
  viewCommand: ViewCommand;
  cinema: CinemaPhase;
  cinemaProgress: number;
  cinemaError: string | null;
  /** Terminal status once playback has stopped (not window_end). */
  cinemaEnd: string | null;
  /** Informational note shown in the cinematic HUD (clamped start date, early stop...). */
  cinemaNote: string | null;
  setData: (data: AsteroidData) => void;
  setError: (message: string) => void;
  toggle: (cat: Category) => void;
  select: (index: number, enableCategory?: boolean) => void;
  setHover: (index: number) => void;
  setFollow: (follow: boolean) => void;
  setSpeedIdx: (idx: number) => void;
  togglePause: () => void;
  view: (kind: "earth" | "wide" | "home") => void;
}

export const useStore = create<State>((set, get) => ({
  data: null,
  loadState: "loading",
  error: null,
  show: { impact: true, approaching: true, receding: true },
  selected: -1,
  hover: -1,
  follow: false,
  flyTick: 0,
  speedIdx: DEFAULT_SPEED_INDEX,
  lastSpeedIdx: DEFAULT_SPEED_INDEX,
  reducedMotion: false,
  viewCommand: { kind: "none", n: 0 },
  cinema: "off",
  cinemaProgress: 0,
  cinemaError: null,
  cinemaEnd: null,
  cinemaNote: null,
  setData: (data) => set({ data, loadState: "ready", error: null }),
  setError: (message) => set({ loadState: "error", error: message }),
  toggle: (cat) =>
    set((s) => {
      const show = { ...s.show, [cat]: !s.show[cat] };
      const sel = s.selected;
      const dropSelection = sel >= 0 && s.data && s.data.list[sel].category === cat && !show[cat];
      return dropSelection ? { show, selected: -1, follow: false } : { show };
    }),
  select: (index, enableCategory = false) =>
    set((s) => {
      if (index < 0 || !s.data) return { selected: -1, follow: false };
      const cat = s.data.list[index].category;
      const show = enableCategory && !s.show[cat] ? { ...s.show, [cat]: true } : s.show;
      return { selected: index, follow: true, flyTick: s.flyTick + 1, show };
    }),
  setHover: (index) => set({ hover: index }),
  setFollow: (follow) => set({ follow }),
  setSpeedIdx: (idx) => set({ speedIdx: idx, lastSpeedIdx: idx === PAUSED_INDEX ? get().lastSpeedIdx : idx }),
  togglePause: () => {
    const { speedIdx, lastSpeedIdx } = get();
    set({ speedIdx: speedIdx === PAUSED_INDEX ? (lastSpeedIdx === PAUSED_INDEX ? DEFAULT_SPEED_INDEX : lastSpeedIdx) : PAUSED_INDEX });
  },
  view: (kind) => set((s) => ({ viewCommand: { kind, n: s.viewCommand.n + 1 }, follow: false })),
}));

// Dev-only handle for poking at the running app from the console.
if (import.meta.env.DEV) (window as unknown as { __aw: unknown }).__aw = { useStore, live, clock };
