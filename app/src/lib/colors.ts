import type { Category } from "./data";

export const CATEGORY_COLOR: Record<Category, string> = {
  impact: "#ff4d5a",
  approaching: "#ffb020",
  receding: "#3cd6f0",
};

/** Same colours as linear-ish floats, indexed like AsteroidData.cat. */
export const CATEGORY_RGB: [number, number, number][] = [
  [1.0, 0.3, 0.35],
  [1.0, 0.69, 0.13],
  [0.24, 0.84, 0.94],
];

export const CATEGORY_LABEL: Record<Category, string> = {
  impact: "Impact-capable (diminishing radius)",
  approaching: "Approaching (closing, no impact path)",
  receding: "Receding (never crash)",
};
