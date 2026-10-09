import { CATEGORY_COLOR, CATEGORY_LABEL } from "../lib/colors";
import { AU_KM, earthDistanceAu } from "../lib/kepler";
import { clock, useStore } from "../lib/store";
import { startCinema } from "../lib/cinematic/controller";
import { fmtDiameter, fmtDistance, fmtKm, fmtPeriod, useTick } from "./format";

/** Real designation: the parenthesised provisional designation of the name, else the name itself. */
function designation(name: string): string {
  const m = /\(([^)]+)\)\s*$/.exec(name);
  return m ? m[1] : name;
}

const DT_DAYS = 0.01;

export function InfoCard() {
  const selected = useStore((s) => s.selected);
  const data = useStore((s) => s.data);
  const follow = useStore((s) => s.follow);
  const select = useStore((s) => s.select);
  const setFollow = useStore((s) => s.setFollow);
  useTick(250);
  if (selected < 0 || !data) return null;

  const a = data.list[selected];
  const now = earthDistanceAu(data.table, selected, clock.jd);
  // Sign of d(distance)/dt at the current clock time, from positions at jd and jd + 0.01 d.
  const rateKms = ((earthDistanceAu(data.table, selected, clock.jd + DT_DAYS) - now) / DT_DAYS) * (AU_KM / 86400);
  const closing = rateKms < 0;
  const color = CATEGORY_COLOR[a.category];

  return (
    <section className="card" aria-label={`Asteroid ${a.name}`} style={{ "--c": color } as React.CSSProperties}>
      <div className="card-head">
        <div>
          <span className="pill">{CATEGORY_LABEL[a.category]}</span>
          {a.pha && <span className="pill pill--warn">Potentially hazardous</span>}
          {a.strict_collision && <span className="pill pill--warn">Collision course</span>}
          <h2>{a.name}</h2>
        </div>
        <button type="button" className="icon-btn" onClick={() => select(-1)} aria-label="Close">
          &#10005;
        </button>
      </div>
      <dl className="card-grid">
        <div>
          <dt>Diameter</dt>
          <dd>{fmtDiameter(a.diameter_km)}</dd>
        </div>
        <div>
          <dt>Distance now</dt>
          <dd>
            {fmtDistance(now)} <small>{fmtKm(now)}</small>
          </dd>
        </div>
        <div>
          <dt>Right now</dt>
          <dd className={closing ? "trend trend--closing" : "trend trend--away"}>
            {closing ? "Closing in" : "Moving away"} <small>{Math.abs(rateKms).toFixed(1)} km/s radial</small>
          </dd>
        </div>
        <div>
          <dt>Closest approach</dt>
          <dd>
            {fmtDistance(a.min_dist_au)} <small>{a.min_date.slice(0, 10)}</small>
          </dd>
        </div>
        <div>
          <dt>Orbital period</dt>
          <dd>{fmtPeriod(a.a)}</dd>
        </div>
        <div>
          <dt>a / e / i</dt>
          <dd>
            {a.a.toFixed(3)} AU / {a.e.toFixed(3)} / {a.i.toFixed(1)}&deg;
          </dd>
        </div>
        <div>
          <dt>Designation</dt>
          <dd>{designation(a.name)}</dd>
        </div>
      </dl>
      <div className="card-actions">
        <button type="button" className="chip chip--accent" onClick={() => void startCinema(selected)}>
          Cinematic &#9654;
        </button>
        <button type="button" className="chip" onClick={() => setFollow(!follow)} aria-pressed={follow}>
          {follow ? "Following" : "Follow"}
        </button>
        <a
          className="link-btn"
          href={`https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?sstr=${encodeURIComponent(a.id)}`}
          target="_blank"
          rel="noreferrer"
        >
          JPL record &#8599;
        </a>
      </div>
      <p className="card-note">Showing the Earth-relative path, 200 days back to 400 days ahead.</p>
    </section>
  );
}
