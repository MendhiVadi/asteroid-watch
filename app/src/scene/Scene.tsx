import { Suspense, useEffect, useRef, useState } from "react";
import { Canvas, useFrame } from "@react-three/fiber";
import { Stars } from "@react-three/drei";
import { Globe } from "./Globe";
import { Moon } from "./Moon";
import { AsteroidField } from "./AsteroidField";
import { OrbitTrail } from "./OrbitTrail";
import { CameraRig } from "./CameraRig";
import { AsteroidActor, CinematicCamera, CinematicDriver, SceneLights, SolarSystem } from "./Cinematic";
import { cine, frame } from "../lib/cinematic/state";
import type { Group } from "three";

// Places the Earth (and Moon) at its heliocentric position in cinematic mode.
function EarthRig({ children }: { children: React.ReactNode }) {
  const ref = useRef<Group>(null);
  useFrame(() => {
    if (!ref.current) return;
    if (cine.active) ref.current.position.copy(frame.earth);
    else ref.current.position.set(0, 0, 0);
  });
  return <group ref={ref}>{children}</group>;
}

// Full-viewport canvas. Lights and camera match the original globe scene
// (ambient 0.6, directional [5,5,5], fov 50); the frame loop pauses whenever
// the scene is scrolled/hidden off-screen.
export default function Scene() {
  const hostRef = useRef<HTMLDivElement>(null);
  const [onScreen, setOnScreen] = useState(true);

  useEffect(() => {
    const node = hostRef.current;
    if (!node || typeof IntersectionObserver === "undefined") return;
    const observer = new IntersectionObserver(([entry]) => setOnScreen(entry.isIntersecting), { threshold: 0 });
    observer.observe(node);
    return () => observer.disconnect();
  }, []);

  return (
    <div ref={hostRef} className="scene" aria-hidden="true">
      <Canvas
        fallback={null}
        dpr={[1, 1.5]}
        frameloop={onScreen ? "always" : "never"}
        camera={{ position: [0, 20, 82], fov: 50, near: 0.05, far: 2000 }}
      >
        <color attach="background" args={["#040a16"]} />
        <SceneLights />
        <Stars radius={400} depth={80} count={3500} factor={7} saturation={0} fade speed={0.4} />
        <CameraRig />
        <CinematicDriver />
        <CinematicCamera />
        <EarthRig>
          <Suspense fallback={null}>
            <Globe />
          </Suspense>
          <Moon />
        </EarthRig>
        <SolarSystem />
        <AsteroidActor />
        <AsteroidField />
        <OrbitTrail />
      </Canvas>
    </div>
  );
}
