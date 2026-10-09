import { existsSync, rmSync } from 'node:fs'
import { fileURLToPath } from 'node:url'
import { resolve } from 'node:path'
import react from '@vitejs/plugin-react'
import { defineConfig } from 'vite'

// moon.jpg is optional: when it exists in public/ the Moon uses it, otherwise a
// procedural cratered texture is generated at runtime (nothing is downloaded).
const hasMoonTexture = existsSync(fileURLToPath(new URL('./public/moon.jpg', import.meta.url)))

// asteroids.mock.json is only fetched when the real dataset fails to load, so it is dropped from
// the production output whenever the real file ships.
const dropMockData = {
  name: 'drop-mock-data',
  apply: 'build' as const,
  closeBundle() {
    const dir = fileURLToPath(new URL('./dist/data', import.meta.url))
    if (existsSync(resolve(dir, 'asteroids.json'))) rmSync(resolve(dir, 'asteroids.mock.json'), { force: true })
  },
}

export default defineConfig({
  plugins: [react(), dropMockData],
  define: { __HAS_MOON_TEXTURE__: JSON.stringify(hasMoonTexture) },
  worker: { format: 'es' },
})
