import { useDeferredValue, useMemo, useState } from "react";
import { CATEGORIES, type Category } from "../lib/data";
import { CATEGORY_COLOR } from "../lib/colors";
import { clock, PAUSED_INDEX, SPEEDS, useStore } from "../lib/store";
import { jdFromMs, msFromJd } from "../lib/kepler";
import { clampClock, clockRange, jdToDateString } from "../lib/coverage";

import { fmtDiameter, fmtInt, fmtSpeed, useTick } from "./format";

const CATEGORY_INFO: Record<Category, { label: string; sub: string }> = {
  impact: { label: "Impact-capable (diminishing radius)", sub: "Distance to Earth shrinks toward an impact path" },
  approaching: { label: "Approaching (closing, no impact path)", sub: "Closing in, passes at a safe distance" },
  receding: { label: "Receding (never crash)", sub: "Moving away from Earth, never crashes" },
};

function Toggle({ cat }: { cat: Category }) {
  const on = useStore((s) => s.show[cat]);
  const toggle = useStore((s) => s.toggle);
  const count = useStore((s) => s.data?.counts[cat] ?? 0);
  const info = CATEGORY_INFO[cat];
  return (
    <button
      type="button"
      role="switch"
      aria-checked={on}
      className={`toggle${on ? " toggle--on" : ""}`}
      style={{ "--c": CATEGORY_COLOR[cat] } as React.CSSProperties}
      onClick={() => toggle(cat)}
    >
      <span className="toggle-track" aria-hidden="true">
        <span className="toggle-thumb" />
      </span>
      <span className="toggle-text">
        <span className="toggle-label">{info.label}</span>
        <span className="toggle-sub">{info.sub}</span>
      </span>
      <span className="toggle-count">{fmtInt(count)}</span>
    </button>
  );
}

function Search() {
  const data = useStore((s) => s.data);
  const select = useStore((s) => s.select);
  const [query, setQuery] = useState("");
  const [open, setOpen] = useState(false);
  const deferred = useDeferredValue(query.trim().toLowerCase());

  const results = useMemo(() => {
    if (!data || deferred.length < 2) return [];
    const out: number[] = [];
    for (let k = 0; k < data.count && out.length < 8; k++) {
      if (data.haystack[k].includes(deferred)) out.push(k);
    }
    return out;
  }, [data, deferred]);

  const choose = (k: number) => {
    select(k, true);
    setQuery("");
    setOpen(false);
  };

  return (
    <div className="search">
      <label className="field-label" htmlFor="asteroid-search">
        Search
      </label>
      <input
        id="asteroid-search"
        type="search"
        placeholder="Name or designation, e.g. 2019 AB"
        autoComplete="off"
        value={query}
        disabled={!data}
        onChange={(e) => {
          setQuery(e.target.value);
          setOpen(true);
        }}
        onFocus={() => setOpen(true)}
        onBlur={() => setTimeout(() => setOpen(false), 120)}
        onKeyDown={(e) => {
          if (e.key === "Enter" && results.length) choose(results[0]);
          if (e.key === "Escape") setOpen(false);
        }}
      />
      {open && deferred.length >= 2 && (
        <ul className="search-results" role="listbox">
          {results.length === 0 && <li className="search-empty">No matches</li>}
          {results.map((k) => {
            const a = data!.list[k];
            return (
              <li key={k} role="option" aria-selected="false">
                <button type="button" onMouseDown={(e) => e.preventDefault()} onClick={() => choose(k)}>
                  <i style={{ background: CATEGORY_COLOR[a.category] }} />
                  <span>{a.name}</span>
                  <small>{fmtDiameter(a.diameter_km)}</small>
                </button>
              </li>
            );
          })}
        </ul>
      )}
    </div>
  );
}

function TimeControls() {
  useStore((s) => s.coverageEnd);
  const speedIdx = useStore((s) => s.speedIdx);
  const setSpeedIdx = useStore((s) => s.setSpeedIdx);
  const togglePause = useStore((s) => s.togglePause);
  useTick(250);
  const date = new Date(msFromJd(clock.jd)).toISOString().slice(0, 10);

  return (
    <div className="time">
      <div className="time-head">
        <span className="field-label">Time</span>
        <span className="time-date">{date}</span>
      </div>
      <div className="time-row">
        <button type="button" className="icon-btn" onClick={togglePause} aria-label={speedIdx === PAUSED_INDEX ? "Play" : "Pause"}>
          {speedIdx === PAUSED_INDEX ? "▶" : "❚❚"}
        </button>
        <input
          type="range"
          min={0}
          max={SPEEDS.length - 1}
          step={1}
          value={speedIdx}
          aria-label="Time speed"
          aria-valuetext={fmtSpeed(SPEEDS[speedIdx])}
          onChange={(e) => setSpeedIdx(Number(e.target.value))}
        />
      </div>
      <div className="time-foot">
        <span>{fmtSpeed(SPEEDS[speedIdx])}</span>
        <span className="time-actions">
          <input
            type="date"
            aria-label="Jump to date"
            min={jdToDateString(clockRange()[0])}
            max={jdToDateString(clockRange()[1])}
            value={date}
            onChange={(e) => {
              const ms = Date.parse(`${e.target.value}T00:00:00Z`);
              if (Number.isFinite(ms)) clock.jd = clampClock(jdFromMs(ms));
            }}
          />
          <button type="button" className="link-btn" onClick={() => (clock.jd = clampClock(jdFromMs(Date.now())))}>
            Now
          </button>
        </span>
      </div>
    </div>
  );
}

export function ControlsPanel() {
  useStore((s) => s.coverageEnd); // re-render when the ephemeris coverage is read
  const data = useStore((s) => s.data);
  const show = useStore((s) => s.show);
  const view = useStore((s) => s.view);
  const [open, setOpen] = useState(() => !matchMedia("(max-width: 820px)").matches);
  const shown = data ? CATEGORIES.reduce((n, c) => n + (show[c] ? data.counts[c] : 0), 0) : 0;

  return (
    <aside className={`panel${open ? "" : " panel--closed"}`} aria-label="Controls">
      <button type="button" className="panel-handle" aria-expanded={open} onClick={() => setOpen((v) => !v)}>
        <span>Filters &amp; time</span>
        <span className="panel-summary">{data ? `${fmtInt(shown)} shown` : "..."}</span>
      </button>
      <div className="panel-body">
        <div className="toggles">
          {CATEGORIES.map((c) => (
            <Toggle key={c} cat={c} />
          ))}
          <p className="toggle-note">Categories are computed for the {jdToDateString(clockRange()[0])} to {jdToDateString(clockRange()[1])} window.</p>
        </div>
        <Search />
        <TimeControls />
        <div className="view-row">
          <span className="field-label">View</span>
          <button type="button" className="chip" onClick={() => view("earth")}>
            Earth
          </button>
          <button type="button" className="chip" onClick={() => view("wide")}>
            Full scene
          </button>
        </div>
      </div>
    </aside>
  );
}
