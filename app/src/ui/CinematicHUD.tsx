import { CATEGORY_COLOR, CATEGORY_LABEL } from "../lib/colors";
import { AU_KM, LUNAR_DISTANCE_AU } from "../lib/kepler";
import { unixMsFromJdTdb } from "../lib/nbody/constants";
import { cine, frame, type ViewMode } from "../lib/cinematic/state";
import { exitCinema, replay, scrub, setPlaying, setView } from "../lib/cinematic/controller";
import { useStore } from "../lib/store";
import { fmtDiameter, useTick } from "./format";
import { useState } from "react";

const RATES = [5, 20, 60, 180, 365];
const fmtUtc = (jdTdb: number) => new Date(unixMsFromJdTdb(jdTdb)).toISOString();
const VIEWS: { id: ViewMode; label: string; hint: string }[] = [
  { id: "auto", label: "Auto", hint: "Pulls back to the solar system, pushes in for the approach" },
  { id: "earth", label: "Earth-centred", hint: "Close-up of the approach" },
  { id: "solar", label: "Solar system", hint: "Heliocentric overview" },
];

function fmtDist(au: number) {
  const km = au * AU_KM;
  const ld = au / LUNAR_DISTANCE_AU;
  const main = ld < 100 ? `${ld.toFixed(ld < 10 ? 3 : 2)} LD` : `${au.toFixed(3)} AU`;
  const sub = km < 1e7 ? `${Math.round(km).toLocaleString("en-US")} km` : `${(km / 1e6).toFixed(1)} million km`;
  return { main, sub };
}

function fmtRate(d: number) {
  if (d >= 1) return `${d.toFixed(d < 10 ? 1 : 0)} d/s`;
  if (d >= 1 / 24) return `${(d * 24).toFixed(1)} h/s`;
  return `${(d * 1440).toFixed(0)} min/s`;
}

function Overlays() {
  const phase = useStore((s) => s.cinema);
  const progress = useStore((s) => s.cinemaProgress);
  const error = useStore((s) => s.cinemaError);
  const end = useStore((s) => s.cinemaEnd);

  if (phase === "loading")
    return (
      <div className="cine-center" role="status">
        <p className="cine-title">Simulating trajectory</p>
        <div className="cine-bar">
          <i style={{ width: `${Math.round(progress * 100)}%` }} />
        </div>
        <p className="cine-sub">Integrating the orbit with Sun, planets, Earth and Moon...</p>
      </div>
    );
  if (phase === "error")
    return (
      <div className="cine-center cine-center--error" role="alert">
        <p className="cine-title">Simulation failed</p>
        <p className="cine-sub">{error}</p>
        <button type="button" className="chip" onClick={exitCinema}>
          Back
        </button>
      </div>
    );
  if (end === "impact_earth" || end === "impact_moon")
    return (
      <>
        <div className="cine-flash" aria-hidden="true" />
        <div className="cine-banner cine-banner--impact" role="alert">
          <p className="cine-title">{end === "impact_earth" ? "Impact" : "Lunar impact"}</p>
          <p className="cine-sub">
            {end === "impact_earth" ? "The asteroid strikes Earth." : "The asteroid strikes the Moon."} Playback stopped.
          </p>
          <div className="cine-actions">
            <button type="button" className="chip" onClick={replay}>
              Replay
            </button>
            <button type="button" className="chip" onClick={exitCinema}>
              Exit
            </button>
          </div>
        </div>
      </>
    );
  if (end === "escaped")
    return (
      <div className="cine-banner" role="status">
        <p className="cine-title">Escaped Earth&apos;s influence</p>
        <p className="cine-sub">The asteroid leaves Earth&apos;s pull and does not return. Playback stopped.</p>
        <div className="cine-actions">
          <button type="button" className="chip" onClick={replay}>
            Replay
          </button>
          <button type="button" className="chip" onClick={exitCinema}>
            Exit
          </button>
        </div>
      </div>
    );
  return null;
}

export function CinematicHUD() {
  const data = useStore((s) => s.data);
  const phase = useStore((s) => s.cinema);
  const progress = useStore((s) => s.cinemaProgress);
  const note = useStore((s) => s.cinemaNote);
  const [view, setViewState] = useState<ViewMode>(cine.view);
  useTick(100);

  const tr = cine.traj;
  const a = data && cine.asteroid >= 0 ? data.list[cine.asteroid] : null;
  const ready = phase === "ready" && tr && a;
  const d = fmtDist(frame.distAu);
  const closest = tr ? fmtDist(tr.closestAu) : null;

  return (
    <div className="cine">
      <div className="cine-top">
        <button type="button" className="chip cine-exit" onClick={exitCinema}>
          &larr; Exit cinematic
        </button>
        {a && (
          <div className="cine-name" style={{ "--c": CATEGORY_COLOR[a.category] } as React.CSSProperties}>
            <span className="pill">{CATEGORY_LABEL[a.category]}</span>
            <strong>{a.name}</strong>
            <small>{fmtDiameter(a.diameter_km)}</small>
          </div>
        )}
      </div>

      {ready && (
        <>
          <div className="cine-hud" aria-live="off">
            <div>
              <span className="field-label">Distance to Earth</span>
              <strong className="cine-big">{d.main}</strong>
              <small>{d.sub}</small>
            </div>
            <div>
              <span className="field-label">Relative speed</span>
              <strong>{frame.speedKms.toFixed(1)} km/s</strong>
            </div>
            <div>
              <span className="field-label">Date (UTC)</span>
              <strong>{fmtUtc(cine.t).slice(0, 16).replace("T", " ")}</strong>
            </div>
            <div>
              <span className="field-label">Closest approach</span>
              <strong>{closest?.main}</strong>
              <small>{fmtUtc(tr.closestJd).slice(0, 10)}</small>
            </div>
            <div>
              <span className="field-label">Playback</span>
              <strong>{fmtRate(cine.rate)}</strong>
              <small>{cine.autoWarp && cine.rate < cine.baseRate - 1e-9 ? "slow-motion" : "base speed"}</small>
            </div>
          </div>

          <div className="cine-controls">
            <div className="cine-row">
              <button
                type="button"
                className="icon-btn"
                onClick={() => setPlaying(!cine.playing)}
                aria-label={cine.playing ? "Pause" : "Play"}
              >
                {cine.playing ? "❚❚" : "▶"}
              </button>
              <input
                type="range"
                aria-label="Scrub time"
                min={tr.tStart}
                max={tr.complete ? tr.tStop : tr.reached}
                step={(tr.tStop - tr.tStart) / 2000}
                value={cine.t}
                onChange={(e) => scrub(Number(e.target.value))}
              />
              <span className="cine-time">
                {Math.round(cine.t - tr.tStart)} / {Math.round((tr.complete ? tr.tStop : tr.reached) - tr.tStart)}
                {tr.complete ? "" : "+"} d
              </span>
            </div>
            <div className="cine-row cine-row--opts">
              <label className="cine-opt">
                Speed
                <select
                  value={RATES.includes(cine.baseRate) ? cine.baseRate : 60}
                  onChange={(e) => {
                    cine.baseRate = Number(e.target.value);
                  }}
                >
                  {RATES.map((r) => (
                    <option key={r} value={r}>
                      {r} d/s
                    </option>
                  ))}
                </select>
              </label>
              <label className="cine-opt">
                <input
                  type="checkbox"
                  defaultChecked={cine.autoWarp}
                  onChange={(e) => {
                    cine.autoWarp = e.target.checked;
                  }}
                />
                Slow near approach
              </label>
              <div className="cine-views" role="group" aria-label="Camera view">
                {VIEWS.map((v) => (
                  <button
                    key={v.id}
                    type="button"
                    className="chip"
                    title={v.hint}
                    aria-pressed={view === v.id}
                    onClick={() => {
                      setView(v.id);
                      setViewState(v.id);
                    }}
                  >
                    {v.label}
                  </button>
                ))}
              </div>
            </div>
            {!tr.complete && (
              <p className="cine-note cine-note--live">
                Integrating the orbit... {Math.round(progress * 100)}% (playback continues as samples arrive)
              </p>
            )}
            {tr.truncated && (
              <p className="cine-note cine-note--warn">Run truncated: the simulation hit its sample/step limit and stops early.</p>
            )}
            {note && <p className="cine-note cine-note--warn">{note}</p>}
            <p className="cine-note">
              Solar-system view: radial distances are log-compressed and planet sizes are not to scale. Earth-centred view
              is close to true scale. Engine: {tr.engine}.
            </p>
          </div>
        </>
      )}
      <Overlays />
    </div>
  );
}
