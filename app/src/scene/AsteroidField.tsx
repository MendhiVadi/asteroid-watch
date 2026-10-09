import { useEffect, useLayoutEffect, useMemo, useRef } from "react";
import { useFrame, useThree } from "@react-three/fiber";
import { Color, DoubleSide, Matrix4, MeshStandardMaterial, type InstancedMesh, type Mesh } from "three";
import { CATEGORIES, type AsteroidData } from "../lib/data";
import { CATEGORY_RGB } from "../lib/colors";
import { PositionEngine } from "../lib/positions";
import { clock, live, useStore } from "../lib/store";
import { cine } from "../lib/cinematic/state";
import { makeRockGeometry } from "./rockGeometry";
import { EARTH_RADIUS } from "./Globe";

function rockMaterial(): MeshStandardMaterial {
  const mat = new MeshStandardMaterial({ color: "#ffffff", flatShading: true, roughness: 0.85, metalness: 0.05 });
  // Self-illuminate a little so the category colour stays readable on the night side.
  mat.onBeforeCompile = (shader) => {
    shader.fragmentShader = shader.fragmentShader.replace(
      "#include <emissivemap_fragment>",
      "#include <emissivemap_fragment>\n  totalEmissiveRadiance += diffuseColor.rgb * 0.38;",
    );
  };
  return mat;
}

/** Renders every asteroid through one InstancedMesh and handles picking. */
export function AsteroidField() {
  const data = useStore((s) => s.data);
  return data ? <Field data={data} /> : null;
}

function Field({ data }: { data: AsteroidData }) {
  const meshRef = useRef<InstancedMesh>(null);
  const show = useStore((s) => s.show);
  const camera = useThree((s) => s.camera);
  const gl = useThree((s) => s.gl);

  const engine = useMemo(() => new PositionEngine(data.table, data.count), [data]);
  // start() is idempotent and dispose() is restartable, so StrictMode's mount/unmount/mount keeps one worker.
  useEffect(() => {
    engine.start();
    return () => engine.dispose();
  }, [engine]);

  const geometry = useMemo(() => makeRockGeometry(data.count > 12000 ? 0 : 1), [data]);
  const material = useMemo(rockMaterial, []);
  useEffect(
    () => () => {
      geometry.dispose();
      material.dispose();
    },
    [geometry, material],
  );

  // Per-instance spin axis / phase / rate, fixed for the session.
  const spin = useMemo(() => {
    const n = data.count;
    const axis = new Float32Array(n * 3);
    const phase = new Float32Array(n);
    const rate = new Float32Array(n);
    let s = 1234567;
    const rnd = () => {
      s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
      return s / 4294967296;
    };
    for (let k = 0; k < n; k++) {
      const z = rnd() * 2 - 1;
      const t = rnd() * Math.PI * 2;
      const r = Math.sqrt(1 - z * z);
      axis[k * 3] = r * Math.cos(t);
      axis[k * 3 + 1] = z;
      axis[k * 3 + 2] = r * Math.sin(t);
      phase[k] = rnd() * Math.PI * 2;
      rate[k] = 0.25 + rnd() * 1.4;
    }
    return { axis, phase, rate };
  }, [data]);

  // Indices of objects whose category is currently switched on.
  const visible = useMemo(() => {
    const idx = new Uint32Array(data.count);
    let n = 0;
    for (let k = 0; k < data.count; k++) {
      if (show[CATEGORIES[data.cat[k]]]) idx[n++] = k;
    }
    return idx.subarray(0, n);
  }, [data, show]);

  const visibleRef = useRef(visible);
  visibleRef.current = visible;

  useLayoutEffect(() => {
    const mesh = meshRef.current;
    if (!mesh) return;
    if (!mesh.instanceColor) mesh.setColorAt(0, new Color());
    const colors = mesh.instanceColor!.array as Float32Array;
    for (let j = 0; j < visible.length; j++) {
      const c = CATEGORY_RGB[data.cat[visible[j]]];
      colors[j * 3] = c[0];
      colors[j * 3 + 1] = c[1];
      colors[j * 3 + 2] = c[2];
    }
    mesh.instanceColor!.needsUpdate = true;
    mesh.count = engine.ready ? visible.length : 0;
    live.visible = visible;
    live.visibleCount = visible.length;
  }, [visible, data, engine]);

  useFrame((state) => {
    const mesh = meshRef.current;
    if (!mesh) return;
    mesh.visible = !cine.active; // cinematic mode shows only the selected asteroid
    if (cine.active) return;
    engine.update(clock.jd);
    if (!engine.ready) return;
    const pos = engine.positions;
    live.positions = pos;

    const cam = state.camera.position;
    const camDist = Math.hypot(cam.x, cam.y, cam.z);
    // Keep distant markers legible when zoomed out.
    const boost = Math.min(2.5, Math.max(1, camDist / 30));
    const t = state.clock.elapsedTime;
    const sel = useStore.getState().selected;
    const vis = visibleRef.current;
    mesh.count = vis.length;
    const arr = mesh.instanceMatrix.array as Float32Array;
    const { axis, phase, rate } = spin;
    const size = data.size;

    for (let j = 0; j < vis.length; j++) {
      const k = vis[j];
      const s = k === sel ? 0 : size[k] * boost; // the selected rock is drawn separately
      const ang = phase[k] + rate[k] * t;
      const c = Math.cos(ang), sn = Math.sin(ang), om = 1 - c;
      const x = axis[k * 3], y = axis[k * 3 + 1], z = axis[k * 3 + 2];
      const o = j * 16;
      arr[o] = (c + x * x * om) * s;
      arr[o + 1] = (x * y * om + z * sn) * s;
      arr[o + 2] = (x * z * om - y * sn) * s;
      arr[o + 3] = 0;
      arr[o + 4] = (x * y * om - z * sn) * s;
      arr[o + 5] = (c + y * y * om) * s;
      arr[o + 6] = (y * z * om + x * sn) * s;
      arr[o + 7] = 0;
      arr[o + 8] = (x * z * om + y * sn) * s;
      arr[o + 9] = (y * z * om - x * sn) * s;
      arr[o + 10] = (c + z * z * om) * s;
      arr[o + 11] = 0;
      arr[o + 12] = pos[k * 3];
      arr[o + 13] = pos[k * 3 + 1];
      arr[o + 14] = pos[k * 3 + 2];
      arr[o + 15] = 1;
    }
    mesh.instanceMatrix.needsUpdate = true;
  });

  // Picking: screen-space nearest instance. Raycasting 40k instances against
  // their meshes would take hundreds of ms per click, this is < 1 ms.
  useEffect(() => {
    const el = gl.domElement;
    const vp = new Matrix4();
    const THRESH_MIN = 11;

    const pick = (clientX: number, clientY: number): number => {
      const pos = live.positions;
      if (!pos) return -1;
      const rect = el.getBoundingClientRect();
      const px = clientX - rect.left, py = clientY - rect.top;
      camera.updateMatrixWorld();
      vp.multiplyMatrices(camera.projectionMatrix, camera.matrixWorldInverse);
      const m = vp.elements;
      const w2 = rect.width / 2, h2 = rect.height / 2;
      const cam = camera.position;
      const camDist = Math.hypot(cam.x, cam.y, cam.z);
      const boost = Math.min(2.5, Math.max(1, camDist / 30));
      const focal = h2 / Math.tan(((camera as unknown as { fov: number }).fov * Math.PI) / 360);
      const vis = live.visible;
      const n = live.visibleCount;
      const size = data.size;
      let best = -1;
      let bestScore = Infinity;
      for (let j = 0; j < n; j++) {
        const k = vis[j];
        const x = pos[k * 3], y = pos[k * 3 + 1], z = pos[k * 3 + 2];
        const cw = m[3] * x + m[7] * y + m[11] * z + m[15];
        if (cw <= 0.01) continue;
        const sx = w2 + ((m[0] * x + m[4] * y + m[8] * z + m[12]) / cw) * w2;
        const sy = h2 - ((m[1] * x + m[5] * y + m[9] * z + m[13]) / cw) * h2;
        const dx = sx - px, dy = sy - py;
        const thresh = Math.max(THRESH_MIN, (size[k] * boost * focal) / cw + 5);
        const d2 = dx * dx + dy * dy;
        if (d2 > thresh * thresh) continue;
        const score = d2 / (thresh * thresh) + cw * 1e-4;
        if (score >= bestScore) continue;
        // Skip anything hidden behind the Earth.
        const lx = x - cam.x, ly = y - cam.y, lz = z - cam.z;
        const len = Math.hypot(lx, ly, lz);
        const dxn = lx / len, dyn = ly / len, dzn = lz / len;
        const b = -(cam.x * dxn + cam.y * dyn + cam.z * dzn);
        const c = cam.x * cam.x + cam.y * cam.y + cam.z * cam.z - EARTH_RADIUS * EARTH_RADIUS;
        const disc = b * b - c;
        if (disc > 0) {
          const tHit = -b - Math.sqrt(disc);
          if (tHit > 0 && tHit < len) continue;
        }
        best = k;
        bestScore = score;
      }
      return best;
    };

    let downX = 0, downY = 0, downT = 0, down = false;
    const onDown = (e: PointerEvent) => {
      if (cine.active) return;
      down = true;
      downX = e.clientX;
      downY = e.clientY;
      downT = performance.now();
    };
    const onUp = (e: PointerEvent) => {
      if (!down) return;
      down = false;
      if (Math.hypot(e.clientX - downX, e.clientY - downY) > 5 || performance.now() - downT > 700) return;
      const k = pick(e.clientX, e.clientY);
      const store = useStore.getState();
      if (k >= 0) store.select(k);
      else if (store.selected >= 0) store.select(-1);
    };

    let raf = 0;
    let lastX = 0, lastY = 0;
    const onMove = (e: PointerEvent) => {
      if (e.buttons !== 0 || cine.active) return;
      lastX = e.clientX;
      lastY = e.clientY;
      if (raf) return;
      raf = requestAnimationFrame(() => {
        raf = 0;
        const k = pick(lastX, lastY);
        const store = useStore.getState();
        if (k !== store.hover) store.setHover(k);
        el.style.cursor = k >= 0 ? "pointer" : "";
      });
    };
    const onLeave = () => useStore.getState().setHover(-1);

    el.addEventListener("pointerdown", onDown);
    el.addEventListener("pointerup", onUp);
    el.addEventListener("pointermove", onMove);
    el.addEventListener("pointerleave", onLeave);
    return () => {
      el.removeEventListener("pointerdown", onDown);
      el.removeEventListener("pointerup", onUp);
      el.removeEventListener("pointermove", onMove);
      el.removeEventListener("pointerleave", onLeave);
      if (raf) cancelAnimationFrame(raf);
    };
  }, [camera, gl, data]);

  return (
    <>
      <instancedMesh
        ref={meshRef}
        args={[geometry, material, data.count]}
        frustumCulled={false}
        raycast={() => null}
      />
      <SelectedRock data={data} />
    </>
  );
}

// High-detail version of the selected asteroid + a screen-space ring marker.
function SelectedRock({ data }: { data: AsteroidData }) {
  const rockRef = useRef<Mesh>(null);
  const ringRef = useRef<Mesh>(null);
  const geometry = useMemo(() => makeRockGeometry(3), []);
  const material = useMemo(() => {
    const m = rockMaterial();
    m.side = DoubleSide;
    return m;
  }, []);
  const color = useMemo(() => new Color(), []);

  useEffect(
    () => () => {
      geometry.dispose();
      material.dispose();
    },
    [geometry, material],
  );

  useFrame((state) => {
    const rock = rockRef.current;
    const ring = ringRef.current;
    const sel = useStore.getState().selected;
    const pos = live.positions;
    if (!rock || !ring) return;
    if (sel < 0 || !pos || cine.active) {
      rock.visible = false;
      ring.visible = false;
      return;
    }
    rock.visible = ring.visible = true;
    const cam = state.camera.position;
    const camDist = Math.hypot(cam.x, cam.y, cam.z);
    const boost = Math.min(2.5, Math.max(1, camDist / 30));
    const x = pos[sel * 3], y = pos[sel * 3 + 1], z = pos[sel * 3 + 2];
    const s = data.size[sel] * boost * 1.15;
    rock.position.set(x, y, z);
    rock.scale.setScalar(s);
    const t = state.clock.elapsedTime;
    rock.rotation.set(t * 0.35, t * 0.5, 0);
    const c = CATEGORY_RGB[data.cat[sel]];
    color.setRGB(c[0], c[1], c[2]);
    material.color.copy(color);
    ring.position.set(x, y, z);
    ring.quaternion.copy(state.camera.quaternion);
    const d = cam.distanceTo(ring.position);
    ring.scale.setScalar(Math.max(s * 1.9, d * 0.028));
  });

  return (
    <>
      <mesh ref={rockRef} geometry={geometry} material={material} visible={false} raycast={() => null} />
      <mesh ref={ringRef} visible={false} renderOrder={10} raycast={() => null}>
        <ringGeometry args={[0.86, 1, 56]} />
        <meshBasicMaterial color="#ffffff" transparent opacity={0.9} depthTest={false} toneMapped={false} />
      </mesh>
    </>
  );
}
