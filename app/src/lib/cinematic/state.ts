import { Vector3 } from "three";
import type { SimEvent } from "../nbody/sim";

export type TerminalStatus = "impact_earth" | "impact_moon" | "escaped" | "window_end";
export type ViewMode = "auto" | "earth" | "solar";

/**
 * Earth-centred trajectory of one asteroid (ecliptic J2000, AU, times are JD TDB), streamed in from the
 * N-body worker. Sampling is non-uniform (1 min inside 0.01 AU, hourly to 0.05 AU, daily beyond), so
 * everything interpolates by time. Arrays are over-allocated; only the first `n` samples are valid.
 */
export interface Trajectory {
  n: number;
  t: Float64Array;
  /** Asteroid position relative to Earth, xyz interleaved. */
  rel: Float32Array;
  /** Asteroid heliocentric position (relative to the Sun), xyz interleaved. */
  helio: Float32Array;
  /** Moon position relative to Earth, xyz interleaved (null if not available). */
  moon: Float32Array | null;
  /** Terminal status; only final once `complete` is true. */
  status: TerminalStatus;
  tStart: number;
  /** Latest simulated time delivered so far. */
  reached: number;
  /** True when the run finished (tStop/status/closest approach are final). */
  complete: boolean;
  /** Last simulated time; playback stops here unless status is window_end. */
  tStop: number;
  truncated: boolean;
  closestAu: number;
  closestJd: number;
  escapeJd?: number;
  events: SimEvent[];
  engine: string;
}

/** Mutable cinematic playback state, read every frame by the scene. */
export const cine = {
  active: false,
  traj: null as Trajectory | null,
  asteroid: -1,
  t: 0,
  playing: true,
  /** Base playback rate in simulated days per real second. */
  baseRate: 60,
  autoWarp: true,
  view: "auto" as ViewMode,
  ended: false,
  endedAt: 0,
  /** Effective rate this frame (days/s) after time-warp. */
  rate: 0,
  needsCameraInit: false,
  /** Bumped whenever the path arrays below grow (the scene re-uploads them). */
  pathVersion: 0,
  /** Precomputed polylines in scene units, xyz interleaved. */
  localPath: new Float32Array(0),
  widePath: new Float32Array(0),
};

/** Derived world-space values for this frame (written by the driver). */
export const frame = {
  /** 0 = Earth-centred scale, 1 = heliocentric compressed. */
  w: 0,
  earth: new Vector3(),
  earthScale: 1,
  ast: new Vector3(),
  astLocal: new Vector3(),
  moon: new Vector3(),
  moonScale: 1,
  sun: new Vector3(),
  sunScale: 30,
  distAu: 0,
  speedKms: 0,
  idx: 0,
};
