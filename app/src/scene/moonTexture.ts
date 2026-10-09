import { useEffect, useState } from "react";
import { CanvasTexture, SRGBColorSpace, Texture } from "three";

declare const __HAS_MOON_TEXTURE__: boolean;

// Deterministic PRNG so the Moon looks the same every load.
function mulberry32(seed: number) {
  return () => {
    seed = (seed + 0x6d2b79f5) | 0;
    let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// Procedural grey cratered surface (colour + height) painted on a canvas.
export function makeProceduralMoon(): { map: CanvasTexture; bump: CanvasTexture } {
  const W = 1024, H = 512;
  const rand = mulberry32(2024);

  // Multi-octave value noise on small grids, bilinearly sampled.
  const octave = (gw: number, gh: number) => {
    const g = new Float32Array((gw + 1) * (gh + 1));
    for (let y = 0; y <= gh; y++) {
      for (let x = 0; x <= gw; x++) g[y * (gw + 1) + x] = x === gw ? g[y * (gw + 1)] : rand();
    }
    return (u: number, v: number) => {
      const fx = u * gw, fy = v * gh;
      const x0 = Math.floor(fx), y0 = Math.floor(fy);
      const tx = fx - x0, ty = fy - y0;
      const sx = tx * tx * (3 - 2 * tx), sy = ty * ty * (3 - 2 * ty);
      const i = y0 * (gw + 1) + x0;
      const a = g[i] + (g[i + 1] - g[i]) * sx;
      const b = g[i + gw + 1] + (g[i + gw + 2] - g[i + gw + 1]) * sx;
      return a + (b - a) * sy;
    };
  };
  const o1 = octave(6, 3), o2 = octave(16, 8), o3 = octave(48, 24), o4 = octave(128, 64);

  const colour = document.createElement("canvas");
  colour.width = W;
  colour.height = H;
  const cctx = colour.getContext("2d")!;
  const img = cctx.createImageData(W, H);
  for (let y = 0; y < H; y++) {
    for (let x = 0; x < W; x++) {
      const u = x / W, v = y / H;
      const n = 0.4 * o1(u, v) + 0.3 * o2(u, v) + 0.2 * o3(u, v) + 0.1 * o4(u, v);
      // Darker "maria" patches where the broad octave is low.
      const mare = Math.max(0, 0.42 - o1(u, v)) * 1.6;
      const g = Math.round(Math.min(235, Math.max(40, 95 + n * 90 - mare * 70)));
      const k = (y * W + x) * 4;
      img.data[k] = g;
      img.data[k + 1] = g - 1;
      img.data[k + 2] = g - 4;
      img.data[k + 3] = 255;
    }
  }
  cctx.putImageData(img, 0, 0);

  const bumpCanvas = document.createElement("canvas");
  bumpCanvas.width = W;
  bumpCanvas.height = H;
  const bctx = bumpCanvas.getContext("2d")!;
  bctx.drawImage(colour, 0, 0);

  // Craters: dark bowl, bright rim; drawn with wrap-around at the seam.
  for (let c = 0; c < 650; c++) {
    const radius = 2 + Math.pow(rand(), 4.2) * 40;
    const cx = rand() * W;
    const cy = (0.1 + 0.8 * rand()) * H;
    for (const ox of [-W, 0, W]) {
      const x = cx + ox;
      if (x + radius < 0 || x - radius > W) continue;
      const bowl = cctx.createRadialGradient(x, cy, radius * 0.15, x, cy, radius);
      bowl.addColorStop(0, "rgba(40,40,42,0.28)");
      bowl.addColorStop(0.75, "rgba(60,60,62,0.14)");
      bowl.addColorStop(0.9, "rgba(235,235,230,0.16)");
      bowl.addColorStop(1, "rgba(235,235,230,0)");
      cctx.fillStyle = bowl;
      cctx.beginPath();
      cctx.arc(x, cy, radius, 0, Math.PI * 2);
      cctx.fill();

      const relief = bctx.createRadialGradient(x, cy, radius * 0.1, x, cy, radius);
      relief.addColorStop(0, "rgba(0,0,0,0.75)");
      relief.addColorStop(0.8, "rgba(40,40,40,0.4)");
      relief.addColorStop(0.92, "rgba(255,255,255,0.6)");
      relief.addColorStop(1, "rgba(255,255,255,0)");
      bctx.fillStyle = relief;
      bctx.beginPath();
      bctx.arc(x, cy, radius, 0, Math.PI * 2);
      bctx.fill();
    }
  }

  const map = new CanvasTexture(colour);
  map.colorSpace = SRGBColorSpace;
  const bump = new CanvasTexture(bumpCanvas);
  return { map, bump };
}

// If a real public/moon.jpg exists at build time it is used (detected in
// vite.config.ts, so there is no 404 probing); otherwise the procedural one.
export function useMoonTextures(): { map: Texture; bump: Texture | null } {
  const [procedural] = useState(makeProceduralMoon);
  const [photo, setPhoto] = useState<Texture | null>(null);

  useEffect(() => {
    if (!__HAS_MOON_TEXTURE__) return;
    const img = new Image();
    img.onload = () => {
      const tex = new Texture(img);
      tex.colorSpace = SRGBColorSpace;
      tex.needsUpdate = true;
      setPhoto(tex);
    };
    img.src = "/moon.jpg";
  }, []);

  return photo ? { map: photo, bump: null } : procedural;
}
