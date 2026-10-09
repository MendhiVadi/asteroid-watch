import { useRef } from "react";
import { useFrame } from "@react-three/fiber";
import { useTexture } from "@react-three/drei";
import { BackSide, SRGBColorSpace, type Group } from "three";
import { useStore } from "../lib/store";
import { cine, frame } from "../lib/cinematic/state";

const AXIAL_TILT = (23.5 * Math.PI) / 180;
// Scene axes are the ecliptic (x, y, z) -> (x, z, -y), so ecliptic north is +Y. Earth's pole sits at
// ecliptic longitude 90 deg, i.e. ecliptic (0, sin e, cos e) = scene (0, cos e, -sin e): the pole
// leans toward scene -Z, which is a rotation of -e about X (not about Z).
const ROTATION_SPEED = 0.18; // rad/s, continuous free spin (no region lock)
export const EARTH_RADIUS = 2;

// Ported from the Land Dispute globe: radius-2 textured sphere on a 23.5 degree
// axial tilt, spinning at 0.18 rad/s, wrapped in a faint back-side blue
// atmosphere. The India lock/scroll hand-off is removed - Earth just spins.
export function Globe() {
  const orientRef = useRef<Group>(null);
  const scaleRef = useRef<Group>(null);
  const spin = useRef(0);
  // colorSpace is applied in the loader callback, before three uploads the
  // texture to the GPU; setting it later renders washed-out sRGB data.
  const earthTexture = useTexture("/earth.jpg", (texture) => {
    const loaded = Array.isArray(texture) ? texture : [texture];
    loaded.forEach((t) => {
      t.colorSpace = SRGBColorSpace;
      t.needsUpdate = true;
    });
  });

  useFrame((_, delta) => {
    scaleRef.current?.scale.setScalar(cine.active ? frame.earthScale : 1);
    if (useStore.getState().reducedMotion) return;
    // Clamp delta so a resumed tab/offscreen renderer does not jump.
    spin.current += Math.min(delta, 0.05) * ROTATION_SPEED; // prograde (west to east)
    if (orientRef.current) orientRef.current.rotation.y = spin.current;
  });

  return (
    <group ref={scaleRef}>
    <group rotation={[-AXIAL_TILT, 0, 0]}>
      <group ref={orientRef}>
        <mesh>
          <sphereGeometry args={[EARTH_RADIUS, 64, 64]} />
          <meshStandardMaterial map={earthTexture} roughness={0.7} metalness={0.05} />
        </mesh>
      </group>
      <mesh scale={1.08}>
        <sphereGeometry args={[EARTH_RADIUS, 48, 48]} />
        <meshBasicMaterial color="#4da6ff" transparent opacity={0.1} side={BackSide} depthWrite={false} />
      </mesh>
    </group>
    </group>
  );
}
