import { NBodyClient } from "../nbody/client";
import type { Ephemeris } from "../nbody/ephemeris";
import { loadEphemeris } from "../nbody/ephemeris";

// Thin owner of the N-body worker (one NBodyClient for the app) and of the
// main-thread copy of the ephemeris tables used to place the Sun and planets.

const BASE_URL = "/data/ephemeris/";

/** Ephemeris coverage, JD TDB (2020-01-01 .. 2036-12-31); refined from the worker's `ready`. */
export const coverage: [number, number] = [2458849.5, 2465058.5];

let client: NBodyClient | null = null;
let ephemeris: Ephemeris | null = null;
let ephemerisPromise: Promise<Ephemeris> | null = null;

export function getClient(): NBodyClient {
  client ??= new NBodyClient({ baseUrl: BASE_URL });
  return client;
}

/** Waits for the worker to load the tables and returns the coverage. */
export async function ready(): Promise<[number, number]> {
  const cov = await getClient().ready;
  coverage[0] = cov[0];
  coverage[1] = cov[1];
  return cov;
}

/** Main-thread ephemeris (Hermite interpolation of planets.bin / moon.bin); null until loaded. */
export function getEphemeris(): Ephemeris | null {
  return ephemeris;
}

export function loadMainEphemeris(): Promise<Ephemeris> {
  ephemerisPromise ??= loadEphemeris(BASE_URL).then((e) => {
    ephemeris = e;
    return e;
  });
  return ephemerisPromise;
}

export function disposeEngine(): void {
  client?.dispose();
  client = null;
}
