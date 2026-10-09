import { Component, type ReactNode } from "react";

// Minimal error boundary: if the 3D scene throws, show the static fallback.
export class ErrorBoundary extends Component<{ fallback: ReactNode; children: ReactNode }, { failed: boolean }> {
  state = { failed: false };
  static getDerivedStateFromError() {
    return { failed: true };
  }
  componentDidCatch(error: unknown) {
    console.error("[scene] failed, showing static view", error);
  }
  render() {
    return this.state.failed ? this.props.fallback : this.props.children;
  }
}
