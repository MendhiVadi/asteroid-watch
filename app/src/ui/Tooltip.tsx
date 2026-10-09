import { useEffect, useRef } from "react";
import { CATEGORY_COLOR } from "../lib/colors";
import { useStore } from "../lib/store";
import { fmtDiameter } from "./format";

export function Tooltip() {
  const hover = useStore((s) => s.hover);
  const data = useStore((s) => s.data);
  const ref = useRef<HTMLDivElement>(null);

  useEffect(() => {
    const move = (e: PointerEvent) => {
      const el = ref.current;
      if (el) el.style.transform = `translate(${e.clientX + 14}px, ${e.clientY + 14}px)`;
    };
    window.addEventListener("pointermove", move);
    return () => window.removeEventListener("pointermove", move);
  }, []);

  if (hover < 0 || !data) return null;
  const a = data.list[hover];
  return (
    <div ref={ref} className="tooltip" role="presentation">
      <i style={{ background: CATEGORY_COLOR[a.category] }} />
      {a.name} <small>{fmtDiameter(a.diameter_km)}</small>
    </div>
  );
}
