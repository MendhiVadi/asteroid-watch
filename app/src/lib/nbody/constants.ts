// Shared constants and time helpers for the N-body engine.
// Units everywhere: AU, day, radian unless a name says otherwise (Km, KmS, Deg).
// Time scale: JD in TDB (what JPL Horizons and SBDB use).

export const DEG = Math.PI / 180;
export const TWO_PI = 2 * Math.PI;
export const AU_KM = 149597870.7;
export const DAY_S = 86400;
export const C_KM_S = 299792.458;
/** Speed of light in AU/day. */
export const C_AU_D = (C_KM_S * DAY_S) / AU_KM;
/** 1 AU/day in km/s. */
export const AUD_TO_KMS = AU_KM / DAY_S;

export const EARTH_RADIUS_KM = 6371.0;
export const MOON_RADIUS_KM = 1737.4;
export const EARTH_RADIUS_AU = EARTH_RADIUS_KM / AU_KM;
export const MOON_RADIUS_AU = MOON_RADIUS_KM / AU_KM;

/** Earth Hill-sphere radius used for the escape test (AU). */
export const HILL_RADIUS_AU = 0.01;
/** Look-ahead "return" radius for the escape test (AU). */
export const RETURN_RADIUS_AU = 0.05;

/** TT - UTC in seconds for 2017-2026+ (37 leap seconds + 32.184). TDB - TT is < 2 ms. */
export const TT_MINUS_UTC_S = 69.184;

/** Unix time (ms, UTC) -> JD (TDB). */
export function jdTdbFromUnixMs(ms: number): number {
  return (ms + TT_MINUS_UTC_S * 1000) / 86400000 + 2440587.5;
}

/** JD (TDB) -> Unix time (ms, UTC). */
export function unixMsFromJdTdb(jd: number): number {
  return (jd - 2440587.5) * 86400000 - TT_MINUS_UTC_S * 1000;
}

export type Vec3 = [number, number, number];
