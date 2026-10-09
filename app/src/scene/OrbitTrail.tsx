import { useEffect, useMemo, useRef } from "react";
import { useFrame } from "@react-three/fiber";
import { BufferAttribute, BufferGeometry, Line, LineBasicMaterial } from "three";
import { samplePath } from "../lib/kepler";
import { CATEGORY_RGB } from "../lib/colors";
import { clock, useStore } from "../lib/store";
import { cine } from "../lib/cinematic/state";

const SAMPLES = 480;
const PAST_DAYS = 200;
const FUTURE_DAYS = 400;
const REGEN_DAYS = 110;

// Earth-relative path of the selected asteroid (what you actually see it do
// around Earth), through the same compressed mapping as the markers. Dim to
// bright along the timeline; regenerated as the simulation clock drifts.
export function OrbitTrail() {
  const state = useRef({ sel: -2, centre: 0 });
  const { line, positions, colors } = useMemo(() => {
    const positions = new Float32Array(SAMPLES * 3);
    const colors = new Float32Array(SAMPLES * 3);
    const geometry = new BufferGeometry();
    geometry.setAttribute("position", new BufferAttribute(positions, 3));
    geometry.setAttribute("color", new BufferAttribute(colors, 3));
    const material = new LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.95, depthWrite: false });
    const line = new Line(geometry, material);
    line.frustumCulled = false;
    line.visible = false;
    line.raycast = () => {};
    return { line, positions, colors };
  }, []);

  useEffect(
    () => () => {
      line.geometry.dispose();
      (line.material as LineBasicMaterial).dispose();
    },
    [line],
  );

  useFrame(() => {
    const { selected, data } = useStore.getState();
    if (selected < 0 || !data || cine.active) {
      line.visible = false;
      state.current.sel = -2;
      return;
    }
    const st = state.current;
    if (st.sel === selected && Math.abs(clock.jd - st.centre) < REGEN_DAYS) return;
    st.sel = selected;
    st.centre = clock.jd;
    const t0 = clock.jd - PAST_DAYS;
    const t1 = clock.jd + FUTURE_DAYS;
    samplePath(data.table, selected, t0, t1, SAMPLES, positions);
    const rgb = CATEGORY_RGB[data.cat[selected]];
    for (let s = 0; s < SAMPLES; s++) {
      const day = t0 + ((t1 - t0) * s) / (SAMPLES - 1) - clock.jd;
      // Past is dimmer; both ends fade out.
      const fade = day < 0 ? 0.25 + 0.5 * (1 + day / PAST_DAYS) : 1 - 0.65 * (day / FUTURE_DAYS);
      colors[s * 3] = rgb[0] * fade;
      colors[s * 3 + 1] = rgb[1] * fade;
      colors[s * 3 + 2] = rgb[2] * fade;
    }
    line.geometry.attributes.position.needsUpdate = true;
    line.geometry.attributes.color.needsUpdate = true;
    line.visible = true;
  });

  return <primitive object={line} />;
}
