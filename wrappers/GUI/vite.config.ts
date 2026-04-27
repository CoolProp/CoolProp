/// <reference types="vitest" />
import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [react()],
  clearScreen: false,
  // Tauri serves the bundled frontend via a custom protocol
  // (https://tauri.localhost on Windows WebView2). The default Vite-emitted
  // absolute asset paths (`/assets/...`) don't resolve under that scheme,
  // leaving Windows users with a blank CoolProp window even though the
  // process stays alive (gh-2825). Relative paths fix it.
  base: "./",
  server: {
    port: 1420,
    strictPort: true,
    watch: { ignored: ["**/src-tauri/**"] },
  },
  test: {
    environment: "jsdom",
    globals: true,
    setupFiles: ["./src/test/setup.ts"],
    include: ["src/**/*.{test,spec}.{ts,tsx}"],
    css: false,
  },
});
