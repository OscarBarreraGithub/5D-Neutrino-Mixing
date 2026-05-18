// @ts-check
import { defineConfig } from 'astro/config';

// Static-only build for Cloudflare Pages.
// KaTeX rendering happens client-side via the layout's auto-render script.
export default defineConfig({
  output: 'static',
  trailingSlash: 'always',
  build: {
    format: 'directory',
  },
  vite: {
    server: {
      fs: {
        // Allow Vite to read content JSON sitting under src/content/
        strict: true,
      },
    },
  },
});
