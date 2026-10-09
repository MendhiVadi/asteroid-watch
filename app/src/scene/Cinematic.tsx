import { useEffect, useMemo, useRef } from "react";
import { useFrame, useThree } from "@react-three/fiber";
import {
  AdditiveBlending,
  BufferAttribute,
  BufferGeometry,
  CanvasTexture,
  Color,
  Line,
  LineBasicMaterial,
  PerspectiveCamera,
  PointLight,
  Vector3,
  type DirectionalLight,
  type Group,
  type Mesh,
  type MeshBasicMaterial,
  type Sprite,
} from "three";
import { AU_KM, cinemaOffsetToScene } from "../lib/kepler";
import { CATEGORY_RGB } from "../lib/colors";
import { cine, frame } from "../lib/cinematic/state";
import { findIndex, sampleAt } from "../lib/cinematic/trajectoryMath";
import { EARTH_PLANET_INDEX, PLANETS, planetHelio, planetOrbit, toWide } from "../lib/cinematic/ephemeris";
import { clock, useStore } from "../lib/store";
import { makeRockGeometry } from "./rockGeometry";

const LOCAL_TRAIL_DAYS = 45;

const smooth = (a: number, b: number, x: number) => {
  const t = Math.min(1, Math.max(0, (x - a) / (b - a)));
  return t * t * (3 - 2 * t);
};

const rel = new Float64Array(3);
const moonRel = new Float64Array(3);
const e3 = new Float64Array(3);
const tmpA = new Float32Array(3);
const tmpB = new Float32Array(3);
const tmpM = new Float32Array(3);
const tmpW = new Float32Array(3);

/** Advances cinematic time (with auto time-warp) and computes this frame's world-space layout. */
export function CinematicDriver() {
  useFrame((_, delta) => {
    const tr = cine.traj;
    if (!cine.active || !tr) return;
    const dt = Math.min(delta, 0.05);

    let s = sampleAt(tr, cine.t, frame.idx, rel, moonRel);
    let dist = Math.hypot(rel[0], rel[1], rel[2]);

    // Rate: base speed, slowed so the remaining approach/departure time is
    // never covered in less than ~4 real seconds (slow-motion near the planet).
    let rate = cine.baseRate;
    if (cine.autoWarp) {
      const vAuDay = (s.speedKms * 86400) / AU_KM;
      const tau = dist / Math.max(vAuDay, 1e-9);
      rate = Math.min(rate, Math.max(tau / 4, 0.002));
    }
    cine.rate = rate;

    if (cine.playing && !cine.ended) {
      cine.t += rate * dt;
      const end = tr.complete ? tr.tStop : tr.reached;
      if (cine.t >= end) {
        if (!tr.complete) {
          cine.t = end; // wait for the integrator to deliver more samples
        } else if (tr.status === "window_end") {
          cine.t = tr.tStart; // loop to the start of the window
        } else {
          cine.t = tr.tStop;
          cine.playing = false;
          cine.ended = true;
          cine.endedAt = performance.now();
          useStore.setState({ cinemaEnd: tr.status });
        }
      }
      s = sampleAt(tr, cine.t, s.idx, rel, moonRel);
      dist = Math.hypot(rel[0], rel[1], rel[2]);
    }
    frame.idx = s.idx;
    frame.distAu = dist;
    frame.speedKms = s.speedKms;

    const wTarget = cine.view === "solar" ? 1 : cine.view === "earth" ? 0 : smooth(0.02, 0.25, dist);
    if (cine.needsCameraInit) frame.w = wTarget;
    frame.w += (wTarget - frame.w) * (1 - Math.exp(-dt * 2.2));
    const w = frame.w;

    // Earth in the heliocentric (wide) frame.
    planetHelio(EARTH_PLANET_INDEX, cine.t, e3);
    toWide(e3[0], e3[1], e3[2], tmpW);
    frame.earth.set(tmpW[0], tmpW[1], tmpW[2]);
    frame.earthScale = 1 + (0.3 - 1) * w;

    // Asteroid: Earth-centred true-scale offset blended with its compressed heliocentric position.
    cinemaOffsetToScene(rel[0], rel[1], rel[2], tmpA);
    frame.astLocal.set(tmpA[0], tmpA[1], tmpA[2]);
    toWide(rel[0] + e3[0], rel[1] + e3[1], rel[2] + e3[2], tmpB);
    frame.ast.set(
      frame.earth.x + tmpA[0] + (tmpB[0] - frame.earth.x - tmpA[0]) * w,
      frame.earth.y + tmpA[1] + (tmpB[1] - frame.earth.y - tmpA[1]) * w,
      frame.earth.z + tmpA[2] + (tmpB[2] - frame.earth.z - tmpA[2]) * w,
    );

    // Moon: true-scale offset, or a fixed small orbit radius in the wide view.
    cinemaOffsetToScene(moonRel[0], moonRel[1], moonRel[2], tmpM);
    const mr = Math.hypot(moonRel[0], moonRel[1], moonRel[2]) || 1;
    const wideR = 2.1 / mr;
    frame.moon.set(
      frame.earth.x + tmpM[0] + (moonRel[0] * wideR - tmpM[0]) * w,
      frame.earth.y + tmpM[1] + (moonRel[2] * wideR - tmpM[1]) * w,
      frame.earth.z + tmpM[2] + (-moonRel[1] * wideR - tmpM[2]) * w,
    );
    frame.moonScale = 1 + (0.45 - 1) * w;

    // Sun: far away along the real sun direction in Earth view, at the origin in the wide view.
    const len = frame.earth.length() || 1;
    const k = 1 - 330 / len;
    frame.sun.set(frame.earth.x * k * (1 - w), frame.earth.y * k * (1 - w), frame.earth.z * k * (1 - w));
    frame.sunScale = 28 + (5 - 28) * w;
  });
  return null;
}

/** Sun-direction light in cinematic mode; the fixed original light otherwise. */
export function SceneLights() {
  const dir = useRef<DirectionalLight>(null);
  const point = useMemo(() => {
    const p = new PointLight("#fff3df", 3.4, 0, 0);
    return p;
  }, []);
  useFrame(() => {
    if (dir.current) dir.current.visible = !cine.active;
    point.visible = cine.active;
    if (cine.active) point.position.copy(frame.sun);
  });
  return (
    <>
      <ambientLight intensity={0.6} />
      <directionalLight ref={dir} position={[5, 5, 5]} intensity={1} />
      <primitive object={point} />
    </>
  );
}

function glowTexture() {
  const c = document.createElement("canvas");
  c.width = c.height = 128;
  const g = c.getContext("2d")!;
  const grad = g.createRadialGradient(64, 64, 0, 64, 64, 64);
  grad.addColorStop(0, "rgba(255,240,200,1)");
  grad.addColorStop(0.25, "rgba(255,200,110,0.55)");
  grad.addColorStop(1, "rgba(255,160,60,0)");
  g.fillStyle = grad;
  g.fillRect(0, 0, 128, 128);
  return new CanvasTexture(c);
}

/** Sun, Mercury-Neptune (Earth is the main Globe) and their orbit lines. */
export function SolarSystem() {
  const sunGroup = useRef<Group>(null);
  const glow = useRef<Sprite>(null);
  const planetsGroup = useRef<Group>(null);
  const planetMeshes = useRef<(Mesh | null)[]>([]);
  const texture = useMemo(glowTexture, []);

  const orbits = useMemo(() => {
    const jd = clock.jd;
    return PLANETS.map((p, i) => {
      const pts = planetOrbit(i, jd, 256);
      const out = new Float32Array(pts.length);
      for (let s = 0; s < pts.length; s += 3) toWide(pts[s], pts[s + 1], pts[s + 2], out, s);
      const g = new BufferGeometry();
      g.setAttribute("position", new BufferAttribute(out, 3));
      const m = new LineBasicMaterial({ color: new Color(p.color), transparent: true, opacity: 0, depthWrite: false });
      const line = new Line(g, m);
      line.frustumCulled = false;
      line.raycast = () => {};
      return line;
    });
  }, []);

  useEffect(
    () => () => {
      orbits.forEach((l) => {
        l.geometry.dispose();
        (l.material as LineBasicMaterial).dispose();
      });
      texture.dispose();
    },
    [orbits, texture],
  );

  useFrame(() => {
    const on = cine.active;
    if (sunGroup.current) {
      sunGroup.current.visible = on;
      sunGroup.current.position.copy(frame.sun);
      sunGroup.current.scale.setScalar(frame.sunScale);
    }
    if (glow.current) glow.current.scale.setScalar(5.5);
    const w = frame.w;
    if (planetsGroup.current) planetsGroup.current.visible = on && w > 0.02;
    if (!on) return;
    orbits.forEach((l, i) => {
      const mat = l.material as LineBasicMaterial;
      mat.opacity = w * (i === EARTH_PLANET_INDEX ? 0.7 : 0.35);
      l.visible = w > 0.02;
    });
    for (let i = 0; i < PLANETS.length; i++) {
      const mesh = planetMeshes.current[i];
      if (!mesh || i === EARTH_PLANET_INDEX) continue;
      planetHelio(i, cine.t, e3);
      toWide(e3[0], e3[1], e3[2], tmpW);
      mesh.position.set(tmpW[0], tmpW[1], tmpW[2]);
    }
  });

  return (
    <>
      <group ref={sunGroup} visible={false}>
        <mesh raycast={() => null}>
          <sphereGeometry args={[1, 40, 40]} />
          <meshBasicMaterial color="#ffd27a" toneMapped={false} />
        </mesh>
        <sprite ref={glow}>
          <spriteMaterial map={texture} blending={AdditiveBlending} depthWrite={false} transparent toneMapped={false} />
        </sprite>
      </group>
      {orbits.map((l, i) => (
        <primitive key={PLANETS[i].name} object={l} />
      ))}
      <group ref={planetsGroup} visible={false}>
        {PLANETS.map((p, i) =>
          i === EARTH_PLANET_INDEX ? null : (
            <mesh
              key={p.name}
              ref={(m) => {
                planetMeshes.current[i] = m;
              }}
              raycast={() => null}
            >
              <sphereGeometry args={[p.radius, 24, 24]} />
              <meshStandardMaterial color={p.color} roughness={0.9} />
            </mesh>
          ),
        )}
      </group>
    </>
  );
}

/** The selected asteroid: tumbling rock, marker ring, Earth-relative and heliocentric trails, impact flash. */
export function AsteroidActor() {
  const data = useStore((s) => s.data);
  const rock = useRef<Mesh>(null);
  const ring = useRef<Mesh>(null);
  const flash = useRef<Mesh>(null);
  const localGroup = useRef<Group>(null);
  const seen = useRef<unknown>(null);
  const seenVersion = useRef(-1);
  const geometry = useMemo(() => makeRockGeometry(3), []);

  const lines = useMemo(() => {
    const mk = (color: string, opacity: number) => {
      const m = new LineBasicMaterial({ color, transparent: true, opacity, depthWrite: false });
      const l = new Line(new BufferGeometry(), m);
      l.frustumCulled = false;
      l.raycast = () => {};
      return l;
    };
    return {
      localPast: mk("#ffffff", 0.9),
      localFuture: mk("#ffffff", 0.22),
      widePast: mk("#ffffff", 0.9),
      wideFuture: mk("#ffffff", 0.22),
    };
  }, []);

  useEffect(
    () => () => {
      geometry.dispose();
      Object.values(lines).forEach((l) => {
        l.geometry.dispose();
        (l.material as LineBasicMaterial).dispose();
      });
    },
    [geometry, lines],
  );

  useFrame((state) => {
    const tr = cine.traj;
    const show = cine.active && tr && data;
    for (const m of [rock.current, ring.current, flash.current]) if (m) m.visible = !!show;
    if (localGroup.current) localGroup.current.visible = !!show;
    lines.widePast.visible = lines.wideFuture.visible = !!show;
    if (!show || !tr || !data) return;

    if (seen.current !== tr) {
      seen.current = tr;
      const rgb = CATEGORY_RGB[data.cat[cine.asteroid]];
      const col = new Color(rgb[0], rgb[1], rgb[2]);
      (lines.localPast.material as LineBasicMaterial).color.copy(col);
      (lines.widePast.material as LineBasicMaterial).color.copy(col);
      (rock.current!.material as MeshBasicMaterial).color.copy(col);
      seenVersion.current = -1;
    }
    if (seenVersion.current !== cine.pathVersion) {
      // The streamed path grew (possibly into new, larger arrays): re-upload it.
      seenVersion.current = cine.pathVersion;
      for (const l of Object.values(lines)) l.geometry.dispose();
      lines.localPast.geometry.setAttribute("position", new BufferAttribute(cine.localPath, 3));
      lines.localFuture.geometry.setAttribute("position", new BufferAttribute(cine.localPath, 3));
      lines.widePast.geometry.setAttribute("position", new BufferAttribute(cine.widePath, 3));
      lines.wideFuture.geometry.setAttribute("position", new BufferAttribute(cine.widePath, 3));
    }

    const w = frame.w;
    const i = frame.idx;
    // The Earth-relative trail only covers +-LOCAL_TRAIL_DAYS: over longer spans it wraps into clutter.
    const i0 = findIndex(tr, cine.t - LOCAL_TRAIL_DAYS, 0);
    const i1 = findIndex(tr, cine.t + LOCAL_TRAIL_DAYS, 0);
    lines.localPast.geometry.setDrawRange(i0, Math.max(0, i - i0 + 1));
    lines.localFuture.geometry.setDrawRange(i, Math.max(0, i1 - i + 1));
    lines.widePast.geometry.setDrawRange(0, i + 1);
    lines.wideFuture.geometry.setDrawRange(i, tr.n - i);
    (lines.localPast.material as LineBasicMaterial).opacity = 0.9 * (1 - w);
    (lines.localFuture.material as LineBasicMaterial).opacity = 0.22 * (1 - w);
    (lines.widePast.material as LineBasicMaterial).opacity = 0.9 * w;
    (lines.wideFuture.material as LineBasicMaterial).opacity = 0.22 * w;
    lines.localPast.visible = lines.localFuture.visible = w < 0.98;
    lines.widePast.visible = lines.wideFuture.visible = w > 0.02;
    // Earth-relative trail is drawn around the Earth's current position.
    if (localGroup.current) localGroup.current.position.copy(frame.earth);

    const camDist = state.camera.position.distanceTo(frame.ast);
    const base = 0.12 + data.size[cine.asteroid] * 1.2;
    const s = Math.max(base, camDist * 0.007);
    const t = state.clock.elapsedTime;
    const r = rock.current!;
    r.position.copy(frame.ast);
    r.scale.setScalar(s);
    r.rotation.set(t * 0.4, t * 0.55, 0);
    const rg = ring.current!;
    rg.position.copy(frame.ast);
    rg.quaternion.copy(state.camera.quaternion);
    rg.scale.setScalar(Math.max(s * 2, camDist * 0.03));

    const f = flash.current!;
    const impact = cine.ended && (tr.status === "impact_earth" || tr.status === "impact_moon");
    f.visible = impact;
    if (impact) {
      const age = (performance.now() - cine.endedAt) / 1000;
      f.position.copy(frame.ast);
      f.scale.setScalar(0.6 + age * 7);
      (f.material as MeshBasicMaterial).opacity = Math.max(0, 1 - age / 2.4);
    }
  });

  return (
    <>
      <group ref={localGroup} visible={false}>
        <primitive object={lines.localFuture} />
        <primitive object={lines.localPast} />
      </group>
      <primitive object={lines.wideFuture} />
      <primitive object={lines.widePast} />
      <mesh ref={rock} geometry={geometry} visible={false} raycast={() => null}>
        <meshBasicMaterial color="#ffffff" toneMapped={false} />
      </mesh>
      <mesh ref={ring} visible={false} renderOrder={10} raycast={() => null}>
        <ringGeometry args={[0.86, 1, 56]} />
        <meshBasicMaterial color="#ffffff" transparent opacity={0.85} depthTest={false} toneMapped={false} />
      </mesh>
      <mesh ref={flash} visible={false} renderOrder={11} raycast={() => null}>
        <sphereGeometry args={[1, 24, 24]} />
        <meshBasicMaterial color="#ffd9a0" transparent opacity={1} blending={AdditiveBlending} depthWrite={false} toneMapped={false} />
      </mesh>
    </>
  );
}

const camDir = new Vector3(0.55, 0.42, 0.72).normalize();
const tmpV = new Vector3();
const cur = new Vector3();
const targetE = new Vector3();
const targetW = new Vector3();
const goalTarget = new Vector3();
const newTarget = new Vector3();

type ControlsLike = { target: Vector3; minDistance: number; maxDistance: number; update: () => void };

/** Auto-framing camera: Earth-centred close-up, solar-system overview, or auto between them. */
export function CinematicCamera() {
  const controls = useThree((s) => s.controls) as unknown as ControlsLike | null;
  const camera = useThree((s) => s.camera) as PerspectiveCamera;
  const size = useThree((s) => s.size);
  const gl = useThree((s) => s.gl);
  const st = useRef({ dSmooth: 0, mul: 1, lastInput: 0, dir: camDir.clone() });

  useEffect(() => {
    const el = gl.domElement;
    const mark = () => {
      st.current.lastInput = performance.now();
    };
    el.addEventListener("pointerdown", mark);
    el.addEventListener("wheel", mark, { passive: true });
    return () => {
      el.removeEventListener("pointerdown", mark);
      el.removeEventListener("wheel", mark);
    };
  }, [gl]);

  useFrame((_, delta) => {
    const c = controls;
    if (!c || !cine.active) return;
    const dt = Math.min(delta, 0.05);
    const s = st.current;
    c.minDistance = 0.5;
    c.maxDistance = 900;

    const aspect = size.width / Math.max(1, size.height);
    const tanHalf = Math.tan((camera.fov * Math.PI) / 360);
    // Frame Earth plus the asteroid, but never pull back past ~40 units: a distant asteroid just leaves the frame.
    const sepFull = frame.astLocal.length();
    const sep = Math.min(sepFull, 40);
    const dEarth = Math.max(7, (sep / 2 + 2.5) / (tanHalf * Math.min(1, aspect)) * 1.2 + 2);
    targetE.copy(frame.earth).addScaledVector(frame.astLocal, sepFull > 1e-6 ? (0.5 * sep) / sepFull : 0);
    targetW.copy(frame.earth).add(frame.ast).multiplyScalar(0.3);
    const dWide = Math.max(120, frame.ast.length() * 2.6, 70 / Math.min(1, aspect));

    const wc = frame.w * frame.w * (3 - 2 * frame.w);
    goalTarget.lerpVectors(targetE, targetW, wc);
    const goalD = Math.exp(Math.log(dEarth) * (1 - wc) + Math.log(dWide) * wc);

    if (cine.needsCameraInit) {
      s.dSmooth = goalD;
      s.mul = 1;
      s.dir.copy(camDir);
      c.target.copy(goalTarget);
      camera.position.copy(goalTarget).addScaledVector(s.dir, goalD);
      c.update();
      cine.needsCameraInit = false;
      return;
    }

    // Preserve the user's orbit direction and zoom ratio.
    cur.copy(camera.position).sub(c.target);
    const curD = cur.length();
    if (curD > 1e-6) s.dir.copy(cur).multiplyScalar(1 / curD);
    if (s.dSmooth > 0 && Math.abs(curD - s.dSmooth) > 1e-3 * s.dSmooth) {
      s.mul *= curD / s.dSmooth; // the user zoomed: keep their ratio to the auto framing
      s.dSmooth = curD;
    }
    if (performance.now() - s.lastInput > 3500) {
      const a = 0.05 * dt;
      const x = s.dir.x * Math.cos(a) + s.dir.z * Math.sin(a);
      const z = -s.dir.x * Math.sin(a) + s.dir.z * Math.cos(a);
      s.dir.x = x;
      s.dir.z = z;
    }
    const k = 1 - Math.exp(-dt * 3.5);
    newTarget.copy(c.target).lerp(goalTarget, k);
    s.dSmooth += (goalD * s.mul - s.dSmooth) * k;
    c.target.copy(newTarget);
    tmpV.copy(s.dir).multiplyScalar(s.dSmooth);
    camera.position.copy(newTarget).add(tmpV);
    c.update();
  });

  return null;
}
