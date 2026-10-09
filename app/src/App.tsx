import { lazy, Suspense, useEffect, useState } from "react";
import { loadAsteroids } from "./lib/data";
import { PAUSED_INDEX, useStore } from "./lib/store";
import { ErrorBoundary } from "./ErrorBoundary";
import { disposeEngine } from "./lib/cinematic/engine";
import { HeroContent } from "./ui/HeroContent";
import { ControlsPanel } from "./ui/ControlsPanel";
import { InfoCard } from "./ui/InfoCard";
import { Tooltip } from "./ui/Tooltip";
import { CinematicHUD } from "./ui/CinematicHUD";

const Scene = lazy(() => import("./scene/Scene"));

function supportsWebGL() {
  try {
    const canvas = document.createElement("canvas");
    const context = canvas.getContext("webgl2") || canvas.getContext("webgl");
    context?.getExtension("WEBGL_lose_context")?.loseContext();
    return Boolean(context);
  } catch {
    return false;
  }
}

// Static fallback: a plain textured Earth disc behind the same overlay UI.
function StaticBackdrop() {
  return (
    <div className="scene scene--static" aria-hidden="true">
      <div className="static-earth" />
    </div>
  );
}

export default function App() {
  const [hasWebGL] = useState(supportsWebGL);
  const setData = useStore((s) => s.setData);
  const setError = useStore((s) => s.setError);
  const cinema = useStore((s) => s.cinema);
  const cinematic = cinema !== "off";

  useEffect(() => {
    const media = matchMedia("(prefers-reduced-motion: reduce)");
    const apply = () => {
      useStore.setState({ reducedMotion: media.matches });
      if (media.matches) useStore.getState().setSpeedIdx(PAUSED_INDEX);
    };
    apply();
    media.addEventListener("change", apply);
    return () => media.removeEventListener("change", apply);
  }, []);

  useEffect(() => () => disposeEngine(), []);

  useEffect(() => {
    let cancelled = false;
    loadAsteroids().then(
      (data) => !cancelled && setData(data),
      (err: unknown) => !cancelled && setError(err instanceof Error ? err.message : String(err)),
    );
    return () => {
      cancelled = true;
    };
  }, [setData, setError]);

  return (
    <main className="app">
      {hasWebGL ? (
        <ErrorBoundary fallback={<StaticBackdrop />}>
          <Suspense fallback={<StaticBackdrop />}>
            <Scene />
          </Suspense>
        </ErrorBoundary>
      ) : (
        <StaticBackdrop />
      )}
      <div className="overlay">
        {cinematic ? (
          <CinematicHUD />
        ) : (
          <>
            <HeroContent />
            <ControlsPanel />
            <InfoCard />
            <p className="hint">
              Drag to orbit &middot; scroll to zoom &middot; click a rock to inspect. Distances are log-compressed around
              Earth; sizes are not to scale.
            </p>
          </>
        )}
      </div>
      {!cinematic && <Tooltip />}
    </main>
  );
}
