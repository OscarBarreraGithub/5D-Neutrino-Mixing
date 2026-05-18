// Minimal CDP screenshot helper that FORCES the mobile viewport via
// Emulation.setDeviceMetricsOverride, so Chrome on macOS stops snapping to its
// 500-pixel minimum window width.  Uses Node's built-in WebSocket (Node 22+).
//
// Usage:
//   node scripts/screenshot_cdp.mjs <url> <outPng> <w> <h> [openAnchorSubstring] [waitMs]
import { spawn } from 'node:child_process';
import { setTimeout as wait } from 'node:timers/promises';
import { writeFileSync } from 'node:fs';

const CHROME = '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome';
const PORT = 9222 + Math.floor(Math.random() * 1000);
const [, , url, out, w, h, openAnchor, waitMs] = process.argv;
const W = parseInt(w, 10);
const H = parseInt(h, 10);

const proc = spawn(CHROME, [
  '--headless=new',
  '--disable-gpu',
  '--no-sandbox',
  '--hide-scrollbars',
  '--user-data-dir=/tmp/cdp-' + PORT,
  '--remote-debugging-port=' + PORT,
  '--window-size=' + W + ',' + H,
  'about:blank',
], { stdio: 'pipe' });

let target = null;
for (let i = 0; i < 80; i++) {
  await wait(120);
  try {
    const r = await fetch('http://127.0.0.1:' + PORT + '/json');
    if (r.ok) {
      const list = await r.json();
      target = list.find((t) => t.type === 'page') || list[0];
      if (target && target.webSocketDebuggerUrl) break;
    }
  } catch {}
}
if (!target) {
  console.error('no devtools target');
  proc.kill('SIGKILL');
  process.exit(1);
}

const ws = new WebSocket(target.webSocketDebuggerUrl);
await new Promise((res, rej) => {
  ws.addEventListener('open', () => res());
  ws.addEventListener('error', rej, { once: true });
});

let nextId = 1;
const pending = new Map();
ws.addEventListener('message', (ev) => {
  const m = JSON.parse(ev.data);
  if (m.id && pending.has(m.id)) {
    const { res, rej } = pending.get(m.id);
    pending.delete(m.id);
    if (m.error) rej(new Error(JSON.stringify(m.error)));
    else res(m.result);
  }
});
function send(method, params = {}) {
  const id = nextId++;
  return new Promise((res, rej) => {
    pending.set(id, { res, rej });
    ws.send(JSON.stringify({ id, method, params }));
  });
}

await send('Page.enable');
await send('Emulation.setDeviceMetricsOverride', {
  width: W,
  height: H,
  deviceScaleFactor: 1,
  mobile: W < 768,
});
await send('Page.navigate', { url });
for (let i = 0; i < 100; i++) {
  await wait(120);
  const r = await send('Runtime.evaluate', { expression: 'document.readyState' });
  if (r.result?.value === 'complete') break;
}
await wait(parseInt(waitMs || '1200', 10));

if (openAnchor) {
  await send('Runtime.evaluate', {
    expression: `(() => { let btn=null; document.querySelectorAll('[data-cite-open]').forEach(b=>{ if(!btn && (b.getAttribute('data-cite-open')||'').indexOf(${JSON.stringify(openAnchor)})!==-1) btn=b; }); if (btn) btn.click(); return btn ? btn.getAttribute('data-cite-open') : 'NONE'; })()`,
  });
  await wait(700);
}

const cap = await send('Page.captureScreenshot', { format: 'png' });
writeFileSync(out, Buffer.from(cap.data, 'base64'));
console.log('saved ' + out);

try { ws.close(); } catch {}
proc.kill('SIGKILL');
process.exit(0);
