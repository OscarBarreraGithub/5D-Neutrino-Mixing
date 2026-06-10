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
      // Allow access through VS Code Remote / SSH port forwarding (otherwise Vite
      // returns "Blocked request. This host is not allowed."). Dev-only setting.
      allowedHosts: true,
      fs: {
        // Allow Vite to read content JSON sitting under src/content/
        strict: true,
      },
    },
  },
});
