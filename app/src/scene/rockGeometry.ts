import { BufferAttribute, IcosahedronGeometry, type BufferGeometry } from "three";

// Low-poly asteroid: an icosahedron whose vertices are pushed in/out by a
// smooth function of position (so shared corners move together - no cracks)
// and squashed into an irregular lump. Unit-ish radius; per-instance scale
// sets the on-screen size.
export function makeRockGeometry(detail: number): BufferGeometry {
  const geo = new IcosahedronGeometry(1, detail);
  const pos = geo.getAttribute("position") as BufferAttribute;
  for (let i = 0; i < pos.count; i++) {
    const x = pos.getX(i), y = pos.getY(i), z = pos.getZ(i);
    const n =
      0.16 * Math.sin(3.1 * x + 1.7 * y) * Math.cos(2.3 * z - 0.9 * x) +
      0.09 * Math.sin(6.3 * y + 4.1 * z + 1.3) +
      0.05 * Math.cos(11 * x - 9 * z + 2.1) * Math.sin(8 * y);
    const r = 0.92 + n;
    pos.setXYZ(i, x * r * 1.12, y * r * 0.78, z * r * 0.95);
  }
  geo.computeVertexNormals();
  return geo;
}
