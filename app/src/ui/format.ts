import { useEffect, useState } from "react";
import { AU_KM, LUNAR_DISTANCE_AU } from "../lib/kepler";

export function useTick(ms: number): number {
  const [n, setN] = useState(0);
  useEffect(() => {
    const id = setInterval(() => setN((v) => v + 1), ms);
    return () => clearInterval(id);
  }, [ms]);
  return n;
}

export const fmtInt = (n: number) => n.toLocaleString("en-US");

export function fmtDiameter(km: number): string {
  if (!km) return "unknown";
  if (km < 1) return `${Math.round(km * 1000).toLocaleString("en-US")} m`;
  return `${km.toFixed(km < 10 ? 2 : 1)} km`;
}

export function fmtDistance(au: number): string {
  if (!au) return "n/a";
  const ld = au / LUNAR_DISTANCE_AU;
  if (ld < 100) return `${ld.toFixed(ld < 10 ? 2 : 1)} LD`;
  return `${au.toFixed(3)} AU`;
}

export function fmtKm(au: number): string {
  const km = au * AU_KM;
  if (km < 1e6) return `${Math.round(km).toLocaleString("en-US")} km`;
  return `${(km / 1e6).toFixed(2)} million km`;
}

export function fmtSpeed(days: number): string {
  if (days === 0) return "Paused";
  const mag = Math.abs(days);
  const label = `${mag < 1 ? mag : mag.toString()} day${mag === 1 ? "" : "s"}/s`;
  return days < 0 ? `Reverse ${label}` : label;
}

export function fmtPeriod(aAu: number): string {
  const days = 365.256 * Math.pow(aAu, 1.5);
  return days > 700 ? `${(days / 365.256).toFixed(2)} yr` : `${Math.round(days)} d`;
}
