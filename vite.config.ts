import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";
import { resolve } from "node:path";

export default defineConfig({
  plugins: [react()],
  base: "./",
  worker: { format: "es" },
  build: { target: "es2022", sourcemap: true, rollupOptions: { input: {
    app: resolve(__dirname, "index.html"), methods: resolve(__dirname, "methods.html"),
    preprocessingMethods: resolve(__dirname, "methods-preprocessing.html"),
    consensusMethods: resolve(__dirname, "methods-consensus.html"), downstreamMethods: resolve(__dirname, "methods-downstream.html"),
  } } },
  server: {
    headers: {
      "Cross-Origin-Opener-Policy": "same-origin",
      "Cross-Origin-Embedder-Policy": "require-corp",
    },
  },
});
