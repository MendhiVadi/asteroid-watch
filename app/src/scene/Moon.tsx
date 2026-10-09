import { useMemo, useRef } from "react";
import { useFrame } from "@react-three/fiber";
import { BufferGeometry, Float32BufferAttribute, type Group, type Mesh } from "three";
import { clock } from "../lib/store";
import { cine, frame } from "../lib/cinematic/state";
import { DEG, MOON_ORBIT_SCENE, moonLongitude } from "../lib/kepler";
import { EARTH_RADIUS } from "./Globe";
import { useMoonTextures } from "./moonTexture";

const MOON_RADIUS = EARTH_RADIUS * 0.2727; // ~0.27 Earth radii
const MOON_INCLINATION = 5.145 * DEG;

// The Moon rides the same compressed distance scale as the asteroids
// (its true 0.00257 AU maps to MOON_ORBIT_SCENE), advances with the simulation
// clock (27.3 d period), and is tidally locked because the mesh is a fixed
// child of the rotating pivot.
export function Moon() {
  const pivot = useRef<Group>(null);
  const tilt = useRef<Group>(null);
  const moonMesh = useRef<Mesh>(null);
  const ringRef = useRef<Group>(null);
  const { map, bump } = useMoonTextures();

  const ring = useMemo(() => {
    const pts: number[] = [];
    for (let i = 0; i <= 128; i++) {
      const a = (i / 128) * Math.PI * 2;
      pts.push(Math.cos(a) * MOON_ORBIT_SCENE, 0, -Math.sin(a) * MOON_ORBIT_SCENE);
    }
    const g = new BufferGeometry();
    g.setAttribute("position", new Float32BufferAttribute(pts, 3));
    return g;
  }, []);

  useFrame(() => {
    const m = moonMesh.current;
    if (!m || !pivot.current || !tilt.current) return;
    if (cine.active) {
      // Cinematic: position comes from the trajectory (relative to the Earth rig).
      tilt.current.rotation.x = 0;
      pivot.current.rotation.y = 0;
      if (ringRef.current) ringRef.current.visible = false;
      m.position.set(frame.moon.x - frame.earth.x, frame.moon.y - frame.earth.y, frame.moon.z - frame.earth.z);
      m.scale.setScalar(frame.moonScale);
      m.rotation.y = Math.atan2(-m.position.z, m.position.x) + Math.PI;
    } else {
      tilt.current.rotation.x = MOON_INCLINATION;
      pivot.current.rotation.y = moonLongitude(clock.jd);
      if (ringRef.current) ringRef.current.visible = true;
      m.position.set(MOON_ORBIT_SCENE, 0, 0);
      m.scale.setScalar(1);
      m.rotation.y = Math.PI;
    }
  });

  return (
    <group ref={tilt} rotation={[MOON_INCLINATION, 0, 0]}>
      <group ref={ringRef}>
        <lineLoop geometry={ring}>
          <lineBasicMaterial color="#5b7aa6" transparent opacity={0.35} />
        </lineLoop>
      </group>
      <group ref={pivot}>
        <mesh ref={moonMesh} position={[MOON_ORBIT_SCENE, 0, 0]} rotation={[0, Math.PI, 0]}>
          <sphereGeometry args={[MOON_RADIUS, 48, 48]} />
          <meshStandardMaterial map={map} bumpMap={bump} bumpScale={bump ? 1.2 : 0} roughness={0.95} metalness={0} />
        </mesh>
      </group>
    </group>
  );
}
