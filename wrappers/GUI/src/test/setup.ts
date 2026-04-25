import "@testing-library/jest-dom/vitest";
import { afterEach, vi } from "vitest";
import { cleanup } from "@testing-library/react";

afterEach(() => {
  cleanup();
  vi.clearAllMocks();
});

// Tauri's IPC layer isn't available in jsdom; tests opt in to specific
// behavior by re-mocking individual functions in the relevant suite.
vi.mock("@tauri-apps/api/core", () => ({
  invoke: vi.fn(async () => {
    throw new Error("invoke not mocked in this test");
  }),
}));

// react-plotly.js pulls in plotly.js-dist-min, which references browser
// APIs jsdom doesn't fully implement. Replace with a stub so component
// trees that include <Plot/> can still mount under test.
vi.mock("react-plotly.js", () => ({
  default: () => null,
}));
