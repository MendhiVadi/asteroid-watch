import { useEffect, useRef, type ComponentRef } from "react";
import { useFrame, useThree } from "@react-three/fiber";
import { OrbitControls } from "@react-three/drei";
import gsap from "gsap";
import { Vector3 } from "three";
import { clock, live, SPEEDS, useStore } from "../lib/store";
import { cine } from "../lib/cinematic/state";
import { clampClock } from "../lib/coverage";

const HOME = new Vector3(0, 4.5, 16);
const WIDE = new Vector3(0, 34, 62);
const EARTH_CLOSE = new Vector3(0, 1.4, 4.4);
const FOLLOW_DISTANCE = 5.5;

// Advances the shared simulation clock (days per real second).
function ClockDriver() {
  useFrame((_, delta) => {
    if (cine.active) return;
    const speed = SPEEDS[useStore.getState().speedIdx];
    // Stay inside [today .. end of ephemeris coverage] so cinematic mode can always start.
    if (speed !== 0) clock.jd = clampClock(clock.jd + speed * Math.min(delta, 0.1));
  });
  return null;
}

// OrbitControls with a wide zoom range (just above the Earth's surface out to
// the whole scene), an intro fly-in, view presets, and optional follow of the
// selected asteroid.
export function CameraRig() {
  const controls = useRef<ComponentRef<typeof OrbitControls>>(null);
  const camera = useThree((s) => s.camera);
  const fly = useRef({ until: 0 });
  const tween = useRef<gsap.core.Tween | null>(null);

  const goTo = (position: Vector3, duration: number) => {
    const c = controls.current;
    if (!c) return;
    tween.current?.kill();
    const tl = gsap.timeline({ onUpdate: () => c.update() });
    tl.to(c.target, { x: 0, y: 0, z: 0, duration, ease: "power2.inOut" }, 0);
    tl.to(camera.position, { x: position.x, y: position.y, z: position.z, duration, ease: "power2.inOut" }, 0);
    tween.current = tl as unknown as gsap.core.Tween;
  };

  // Intro: start far out and fly in to the home view.
  useEffect(() => {
    camera.position.set(0, 20, 82);
    const reduced = useStore.getState().reducedMotion;
    const id = requestAnimationFrame(() => goTo(HOME, reduced ? 0.01 : 3));
    return () => {
      cancelAnimationFrame(id);
      tween.current?.kill();
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  useEffect(
    () =>
      useStore.subscribe((s, prev) => {
        if (s.viewCommand.n !== prev.viewCommand.n) {
          if (s.viewCommand.kind === "home") goTo(HOME, 1.4);
          else if (s.viewCommand.kind === "earth") goTo(EARTH_CLOSE, 1.6);
          else if (s.viewCommand.kind === "wide") goTo(WIDE, 1.8);
        }
        if (s.flyTick !== prev.flyTick) {
          tween.current?.kill();
          fly.current.until = performance.now() + 1400;
        }
      }),
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [],
  );

  const delta = useRef(new Vector3());
  const dir = useRef(new Vector3());
  useFrame((_, dt) => {
    const c = controls.current;
    if (!c || cine.active) return;
    c.maxDistance = 95;
    const { selected, follow } = useStore.getState();
    const following = follow && selected >= 0 && live.positions;
    // Allow getting close to a small rock, or to the Earth's surface otherwise.
    c.minDistance = following ? 0.35 : 2.6;
    if (!following) return;
    const p = live.positions!;
    const target = c.target;
    const k = 1 - Math.exp(-dt * 7);
    delta.current.set(p[selected * 3] - target.x, p[selected * 3 + 1] - target.y, p[selected * 3 + 2] - target.z).multiplyScalar(k);
    target.add(delta.current);
    camera.position.add(delta.current);
    if (performance.now() < fly.current.until) {
      dir.current.copy(camera.position).sub(target);
      const dist = dir.current.length();
      const nd = dist + (FOLLOW_DISTANCE - dist) * k;
      camera.position.copy(target).addScaledVector(dir.current.normalize(), nd);
    }
    c.update();
  });

  return (
    <>
      <ClockDriver />
      <OrbitControls
        ref={controls}
        makeDefault
        enableDamping
        dampingFactor={0.08}
        minDistance={2.6}
        maxDistance={95}
        zoomSpeed={0.9}
        rotateSpeed={0.6}
        panSpeed={0.7}
        onStart={() => {
          // The user took over: stop automated camera motion.
          tween.current?.kill();
          fly.current.until = 0;
        }}
      />
    </>
  );
}
