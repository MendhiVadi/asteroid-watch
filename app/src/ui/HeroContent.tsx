import { useStore } from "../lib/store";
import { fmtDiameter, fmtDistance, fmtInt } from "./format";

// Title + headline stats, always rendered (also in the static / no-WebGL view).
export function HeroContent() {
  const data = useStore((s) => s.data);
  const loadState = useStore((s) => s.loadState);
  const error = useStore((s) => s.error);

  return (
    <header className="hero">
      <p className="eyebrow">Near-Earth object monitor</p>
      <h1>Asteroid Watch</h1>
      <p className="hero-description">
        Every known near-Earth asteroid, propagated along its orbit around a spinning Earth and Moon.
        Drag to orbit, scroll to zoom, click a rock to follow its path.
      </p>
      {loadState === "loading" && <p className="status">Loading asteroid catalogue...</p>}
      {loadState === "error" && <p className="status status--error">Could not load asteroid data: {error}</p>}
      {data && (
        <dl className="stats">
          <div>
            <dt>Tracked</dt>
            <dd>{fmtInt(data.count)}</dd>
          </div>
          <div>
            <dt>Hazardous</dt>
            <dd>{fmtInt(data.stats.pha)}</dd>
          </div>
          <div>
            <dt>Closest pass</dt>
            <dd>{data.stats.closest ? fmtDistance(data.stats.closest.min_dist_au) : "n/a"}</dd>
          </div>
          <div>
            <dt>Largest</dt>
            <dd>{data.stats.largest ? fmtDiameter(data.stats.largest.diameter_km) : "n/a"}</dd>
          </div>
        </dl>
      )}
      {data?.source === "mock" && (
        <p className="badge-mock">Demo data: randomly generated placeholder orbits, not real asteroids.</p>
      )}
    </header>
  );
}
