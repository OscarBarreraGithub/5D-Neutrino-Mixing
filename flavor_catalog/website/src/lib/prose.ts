/**
 * texProseToHtml(raw)
 *
 * Pre-processes LaTeX text-mode commands in prose strings so that they render
 * correctly in HTML, before KaTeX handles any embedded math delimiters.
 *
 * Handled commands (outside math mode):
 *   \texttt{X}        -> <code>X</code>
 *   \textbf{X}        -> <strong>X</strong>
 *   \textit{X}        -> <em>X</em>
 *   \textsf{X}        -> <span class="sf">X</span>
 *   \textsc{X}        -> <span class="sc">X</span>
 *   \mbox{X}, \text{X}-> X (transparent wrapper, contents pass through)
 *   \emph{X}          -> <em>X</em>
 *   \href{URL}{TEXT}  -> <a href="URL">TEXT</a>  (http(s) URLs open in new tab)
 *   \url{URL}         -> <a href="URL">URL</a>
 *   \noindent, \par, \hfill, \medskip, \smallskip, \bigskip
 *                     -> stripped (paragraphing handled by HTML)
 *
 * Also: em-dashes outside math spans are normalised to ", " so the prose
 * does not read like LLM output ("a sentence -- more thoughts" wins).
 *
 * Math spans ($...$, \(...\), \[...\], $$...$$) are skipped so their
 * interior is not altered.
 */

// Match a balanced brace group: {…}  — handles one level of nesting.
// For the prose content in this catalog, one level is sufficient.
function extractBraceGroup(s: string, start: number): { content: string; end: number } | null {
  if (s[start] !== '{') return null;
  let depth = 0;
  let i = start;
  while (i < s.length) {
    if (s[i] === '{') depth++;
    else if (s[i] === '}') {
      depth--;
      if (depth === 0) return { content: s.slice(start + 1, i), end: i + 1 };
    }
    i++;
  }
  return null; // unbalanced
}

function isHttpUrl(url: string): boolean {
  return /^https?:\/\//i.test(url.trim());
}

/**
 * Convert LaTeX \begin{itemize}/\end{itemize} (and \begin{enumerate}) into
 * HTML <ul>/<ol> with <li> entries.  Done as a textual pre-pass before the
 * walker so the rest of the helper does not need to know about list
 * environments.
 */
function expandListEnvironments(raw: string): string {
  let out = raw;
  out = out.replace(
    /\\begin\{itemize\}([\s\S]*?)\\end\{itemize\}/g,
    (_m, body) => `<ul>${convertItems(body)}</ul>`,
  );
  out = out.replace(
    /\\begin\{enumerate\}([\s\S]*?)\\end\{enumerate\}/g,
    (_m, body) => `<ol>${convertItems(body)}</ol>`,
  );
  return out;
}

function convertItems(body: string): string {
  // Split on \item, drop leading whitespace, wrap each non-empty piece.
  const pieces = body.split(/\\item\b/g).map((p) => p.trim()).filter(Boolean);
  return pieces.map((p) => `<li>${p}</li>`).join('');
}

export function texProseToHtml(raw: string): string {
  raw = expandListEnvironments(raw);
  // We walk through the string character by character, skipping math spans so
  // we don't corrupt their content.  Outside math spans we replace text-mode
  // LaTeX commands with HTML equivalents.

  const mathOpeners: Array<{ open: string; close: string }> = [
    { open: '$$', close: '$$' },
    { open: '$', close: '$' },
    { open: '\\[', close: '\\]' },
    { open: '\\(', close: '\\)' },
  ];

  let out = '';
  let i = 0;
  const len = raw.length;

  while (i < len) {
    // Check if we're at the start of a math span.
    let mathMatch: { open: string; close: string } | null = null;
    for (const m of mathOpeners) {
      if (raw.startsWith(m.open, i)) {
        mathMatch = m;
        break;
      }
    }

    if (mathMatch) {
      // Emit the math span verbatim (KaTeX handles it later).
      const closeIdx = raw.indexOf(mathMatch.close, i + mathMatch.open.length);
      if (closeIdx === -1) {
        // No closing delimiter found — emit rest verbatim.
        out += raw.slice(i);
        break;
      }
      out += raw.slice(i, closeIdx + mathMatch.close.length);
      i = closeIdx + mathMatch.close.length;
      continue;
    }

    // Check for \texttt, \textbf, \textit, \emph, \href, \url
    if (raw[i] === '\\') {
      let handled = false;

      const commandPatterns: Array<{ name: string; tag: string; attrs?: string }> = [
        { name: 'texttt', tag: 'code' },
        { name: 'textbf', tag: 'strong' },
        { name: 'textit', tag: 'em' },
        { name: 'textsf', tag: 'span', attrs: 'class="tex-sf"' },
        { name: 'textsc', tag: 'span', attrs: 'class="tex-sc"' },
        { name: 'mbox',   tag: 'span' },
        { name: 'text',   tag: 'span' },
        { name: 'emph',   tag: 'em' },
      ];

      for (const { name, tag, attrs } of commandPatterns) {
        if (raw.startsWith('\\' + name + '{', i) || raw.startsWith('\\' + name + ' {', i)) {
          // Find the opening brace.
          let braceStart = i + 1 + name.length;
          while (braceStart < len && raw[braceStart] === ' ') braceStart++;
          const group = extractBraceGroup(raw, braceStart);
          if (group) {
            // Recursively apply transform to the group content.
            const inner = texProseToHtml(group.content);
            const attrStr = attrs ? ' ' + attrs : '';
            out += `<${tag}${attrStr}>${inner}</${tag}>`;
            i = group.end;
            handled = true;
            break;
          }
        }
      }

      // Bare layout commands that should be stripped (no braced argument).
      if (!handled) {
        const layoutCommands = [
          '\\noindent', '\\par', '\\hfill',
          '\\medskip', '\\smallskip', '\\bigskip', '\\indent',
        ];
        for (const cmd of layoutCommands) {
          if (raw.startsWith(cmd, i)) {
            // Skip past command + optional following whitespace; emit a single
            // space so words don't fuse together.
            i += cmd.length;
            if (raw[i] === ' ' || raw[i] === '\n') {
              out += ' ';
            }
            handled = true;
            break;
          }
        }
      }

      if (!handled && raw.startsWith('\\href{', i)) {
        // \href{URL}{TEXT}
        const urlGroup = extractBraceGroup(raw, i + 5);
        if (urlGroup) {
          const textGroup = extractBraceGroup(raw, urlGroup.end);
          if (textGroup) {
            const url = urlGroup.content.trim();
            const text = texProseToHtml(textGroup.content);
            const targetAttr = isHttpUrl(url) ? ' target="_blank" rel="noopener noreferrer"' : '';
            out += `<a href="${url}"${targetAttr}>${text}</a>`;
            i = textGroup.end;
            handled = true;
          }
        }
      }

      if (!handled && raw.startsWith('\\url{', i)) {
        // \url{URL}
        const urlGroup = extractBraceGroup(raw, i + 4);
        if (urlGroup) {
          const url = urlGroup.content.trim();
          const targetAttr = isHttpUrl(url) ? ' target="_blank" rel="noopener noreferrer"' : '';
          out += `<a href="${url}"${targetAttr}>${url}</a>`;
          i = urlGroup.end;
          handled = true;
        }
      }

      if (!handled) {
        out += raw[i];
        i++;
      }
      continue;
    }

    // Em-dash normalisation (only outside math spans, which were emitted above).
    // Treat the em-dash like a soft pause: " — " -> ", " and a bare "—" -> ", ".
    if (raw[i] === '—') {
      // Squash adjacent whitespace into a single ", " so " word — word " reads
      // as "word, word" with no double spacing.
      while (out.endsWith(' ')) out = out.slice(0, -1);
      out += ', ';
      i++;
      // Skip one trailing space if present so we don't get ",  " from " — ".
      if (raw[i] === ' ') i++;
      continue;
    }

    out += raw[i];
    i++;
  }

  // Post-pass: jargon scrub on the assembled prose.  Internal vocabulary (PKA,
  // WA, CA, Wave-N) is replaced with neutral terms; "anarchic-flavor pipeline"
  // is shortened to "RS scan"; the unicode left/right single quotes from
  // copy-paste tex sources are normalised.  These are run on the full output
  // because the math spans we passed verbatim never contain English jargon.
  return scrubProseJargon(out);
}

/**
 * scrubProseJargon(html)
 *
 * Lightweight post-pass that rewrites internal vocabulary in the rendered
 * prose so the reader does not see project-side workflow terms.  Only safe
 * substitutions are applied (whole-word matches, never inside HTML
 * attributes, code spans, or URLs).
 */
function scrubProseJargon(html: string): string {
  // Skip <code>, <a>, <pre>, and math spans ($...$, \(..\), \[..\], $$..$$)
  // so we don't rewrite identifiers, links, or already-correct LaTeX.
  const SKIP_TAGS = new Set(['code', 'pre', 'a']);
  let out = '';
  let i = 0;
  const n = html.length;
  const tagStack: string[] = [];

  const mathOpeners: Array<{ open: string; close: string }> = [
    { open: '$$', close: '$$' },
    { open: '\\[', close: '\\]' },
    { open: '\\(', close: '\\)' },
    { open: '$', close: '$' },
  ];

  while (i < n) {
    // Math span passthrough.
    let mathHit: { open: string; close: string } | null = null;
    for (const m of mathOpeners) {
      if (html.startsWith(m.open, i)) {
        mathHit = m;
        break;
      }
    }
    if (mathHit) {
      const closeIdx = html.indexOf(mathHit.close, i + mathHit.open.length);
      if (closeIdx < 0) {
        out += html.slice(i);
        break;
      }
      out += html.slice(i, closeIdx + mathHit.close.length);
      i = closeIdx + mathHit.close.length;
      continue;
    }
    if (html[i] === '<') {
      const tagEnd = html.indexOf('>', i);
      if (tagEnd < 0) { out += html.slice(i); break; }
      const tag = html.slice(i, tagEnd + 1);
      out += tag;
      const m = tag.match(/^<\/?([a-z]+)\b/i);
      if (m) {
        const name = m[1].toLowerCase();
        if (SKIP_TAGS.has(name)) {
          if (tag.startsWith('</')) {
            const idx = tagStack.lastIndexOf(name);
            if (idx >= 0) tagStack.splice(idx, 1);
          } else if (!tag.endsWith('/>')) {
            tagStack.push(name);
          }
        }
      }
      i = tagEnd + 1;
      continue;
    }
    // Find next tag, math opener, or end.
    let chunkEnd = n;
    const nextTag = html.indexOf('<', i);
    if (nextTag >= 0) chunkEnd = nextTag;
    for (const m of mathOpeners) {
      const idx = html.indexOf(m.open, i);
      if (idx > i && idx < chunkEnd) chunkEnd = idx;
    }
    let chunk = html.slice(i, chunkEnd);
    if (tagStack.length === 0) {
      chunk = applyJargonReplacements(chunk);
    }
    out += chunk;
    i = chunkEnd;
  }
  return out;
}

const JARGON_REPLACEMENTS: Array<[RegExp, string]> = [
  // Anarchic-flavor pipeline -> RS scan (already conceptually shortened on the page H2,
  // but it appears in many prose sections from .tex sources).
  [/\banarchic-flavor pipeline\b/g, 'RS scan'],
  [/\banarchic flavor pipeline\b/g, 'RS scan'],
  // Internal owner / reviewer abbreviations: PKA = "primary author", WA = "writer A",
  // CA = "checker A".  Replace with neutral nouns.  Only match standalone words.
  [/\bThe PKA\b/g, 'The author'],
  [/\bthe PKA\b/g, 'the author'],
  [/\bPKA\b/g, 'the author'],
  [/\bWave-?(\d+)\b/g, 'this catalog wave'],
  // Soften "rc1.1 quark scan" / "rc1.1" references.
  [/\brc1\.1 quark scan\b/g, 'the quark scan in this repository'],
  [/\brc1\.1\b/g, 'this repository'],
];

function applyJargonReplacements(text: string): string {
  let out = text;
  for (const [pat, rep] of JARGON_REPLACEMENTS) {
    out = out.replace(pat, rep);
  }
  out = wrapPhysicsShorthand(out);
  return out;
}

/**
 * wrapPhysicsShorthand(text)
 *
 * Detect well-known physics shorthand in plain prose (Delta C=2, b -> s gamma,
 * nu nubar, ...) and wrap each recognised fragment in $...$ so KaTeX renders
 * it.  Conservative: only fires on the small set of patterns the catalog
 * authors actually use unwrapped.
 */
const SHORTHAND_PATTERNS: Array<[RegExp, (m: RegExpExecArray) => string]> = [
  // Delta {C,F,S,B,M,m}=N  (with optional sign).  Allow "Delta C=2", "Delta F = 2", etc.
  [/(?<![\\$])\bDelta\s*([A-Za-z])\s*=\s*([0-9]+)\b/g, (m) => `$\\Delta ${m[1]} = ${m[2]}$`],
  // Delta m_X / Delta Gamma_X / Delta y / Delta A_CP / Delta L (without =)
  [/(?<![\\$])\bDelta\s+m_([A-Za-z]+)\b/g, (m) => `$\\Delta m_{${m[1]}}$`],
  [/(?<![\\$])\bDelta\s+Gamma\b/g, () => `$\\Delta\\Gamma$`],
  [/(?<![\\$])\bDelta\s+y(?![A-Za-z])/g, () => `$\\Delta y$`],
  [/(?<![\\$])\bDelta\s+A_([A-Za-z]+)\b/g, (m) => `$\\Delta A_{${m[1]}}$`],
  [/(?<![\\$])\bDelta\s+a(?![A-Za-z])/g, () => `$\\Delta a$`],
  [/(?<![\\$])\bDelta\s+L(?![A-Za-z])/g, () => `$\\Delta L$`],
  // Particle arrow patterns: x -> y z (a few tokens, ASCII).  RHS tokens are
  // greedily matched then post-filtered: we drop trailing tokens that look
  // like English prose so "D+->pi+mu+mu- probes c->u" doesn't gobble "probes".
  // Multi-final arrow patterns: "tau -> 3 mu" (decimal + particle, "3μ").
  [
    /\b([A-Za-z][A-Za-z0-9_+\-*\\^/]*)\s*->\s*([0-9])\s+(gamma|mu|nu|tau|pi|phi|psi|ell|e)\b/g,
    (m) => `$${normalizeShorthandToken(m[1] as unknown as string)} \\to ${m[2] as unknown as string}\\${m[3] as unknown as string}$`,
  ],
  // Arrow patterns: greedy match the first RHS token only.  Multi-token tails
  // are handled by an iterative pass after the initial wrap, so we never gobble
  // prose words like "and" between two arrow expressions.
  [
    /\b([A-Za-z][A-Za-z0-9_+\-*\\^/]*)\s*->\s*([A-Za-z0-9][A-Za-z0-9_+\-*\\^/]*)(?=[\s.,;:!?)<]|$)/g,
    (m) => {
      const lhs = m[1] as unknown as string;
      const rhs = m[2] as unknown as string;
      const fullMatch = m[0] as unknown as string;
      if (!isParticleToken(rhs)) return fullMatch;
      return `$${normalizeShorthandToken(lhs)} \\to ${normalizeShorthandToken(rhs)}$`;
    },
  ],
  // Arrow with no LHS prose token: e.g. "$K_L$ -> pi0 gamma" after the math
  // span was passed through.
  [
    /(^|[\s>])->\s*([A-Za-z0-9][A-Za-z0-9_+\-*\\^/]*)(?=[\s.,;:!?)<]|$)/g,
    (m) => {
      const lead = m[1] as unknown as string;
      const rhs = m[2] as unknown as string;
      const fullMatch = m[0] as unknown as string;
      if (!isParticleToken(rhs)) return fullMatch;
      return `${lead}$\\to ${normalizeShorthandToken(rhs)}$`;
    },
  ],
  // Z b bbar -> $Z b \bar{b}$, Z c cbar -> $Z c \bar{c}$
  [/\bZ\s+([bct])\s+\1bar\b/g, (m) => `$Z ${m[1]} \\bar{${m[1]}}$`],
  // Standalone "x xbar" pairs (b bbar, c cbar, t tbar)
  [/\b([bct])\s+\1bar\b/g, (m) => `$${m[1]} \\bar{${m[1]}}$`],
  // Bare "nu nubar"
  [/\bnu\s+nubar\b/g, () => `$\\nu \\bar{\\nu}$`],
  // After an arrow wrap "$X \\to Y$ gamma" -> "$X \\to Y \\gamma$".  Merges
  // trailing Greek shorthand into the preceding math span.
  [
    /\$([^$]+)\$\s+(gamma|pi|mu|nu|tau|phi|psi|ell|rho|omega|eta)(?=[\s.,;:!?)<]|$)/g,
    (m) => `$${m[1] as unknown as string} \\${m[2] as unknown as string}$`,
  ],
  [
    /\$([^$]+)\$\s+(nubar|bbar|cbar|tbar)(?=[\s.,;:!?)<]|$)/g,
    (m) => {
      const sym = (m[2] as unknown as string).charAt(0);
      return `$${m[1] as unknown as string} \\bar{${sym === 'n' ? '\\nu' : sym}}$`;
    },
  ],
];

const PARTICLE_NAMES = new Set([
  'b', 'c', 'd', 's', 't', 'u', 'e', 'p', 'n', 'W', 'Z',
  'gamma', 'mu', 'nu', 'tau', 'pi', 'phi', 'psi', 'ell',
  'nubar', 'bbar', 'cbar', 'tbar',
]);

const GREEK_PARTICLES = '(?:gamma|mu|nu|tau|pi|phi|psi|ell|rho|omega|eta|sigma|kappa|chi|lambda|xi|delta|alpha|beta|theta|Lambda|Sigma|Xi|Omega|Phi|Psi|Delta|Gamma|Upsilon)';

function isParticleToken(tok: string): boolean {
  if (!tok) return false;
  // Slash alternation: "rho/omega", "K*/K", etc.
  if (tok.includes('/')) {
    return tok.split('/').every(isParticleToken);
  }
  // Concatenated finals with charge signs: "pi+mu+mu-", "e+e-", etc.
  if (/^(?:[A-Za-z][a-z]?[+\-*][0-9]?)+$/.test(tok)) return true;
  // Mixed concatenation that ends with a particle name (no trailing sign):
  // "pi+e mu" was the previous match attempt; "pi+e" by itself qualifies.
  if (/^[A-Za-z][a-z]?[+\-*][A-Za-z][a-z]?$/.test(tok)) return true;
  // Two-quark pair like "bs", "bd", "sd" used in FCNC shorthand.
  if (/^[bsdctu]{2}$/.test(tok)) return true;
  const bare = tok.replace(/[+\-*]+$/, '');
  // Multiplicity prefix (3e, 2mu, 3mu, ...).
  if (new RegExp(`^[2-9](?:e|${GREEK_PARTICLES.slice(3, -1)})[+\\-*]?$`).test(tok)) return true;
  // Subscripted forms (B_s, K_L, mu_R, ...)
  if (/^[A-Za-z][A-Za-z]?_[A-Za-z0-9]+[+\-*]?$/.test(tok)) return true;
  // Plain Greek particle names with optional charge sign / digit.
  if (new RegExp(`^${GREEK_PARTICLES}[+\\-*]?[0-9]?$`).test(bare)) return true;
  if (/^(?:nubar|bbar|cbar|tbar)[+\-*]?[0-9]?$/.test(bare)) return true;
  // Single-letter particles with optional charge sign / digit.
  if (/^[A-Z]?[a-zA-Z][+\-*]?[0-9]?$/.test(bare) && PARTICLE_NAMES.has(bare.replace(/[+\-*0-9]+$/, ''))) return true;
  // Backslash escaped already
  if (/^\\[A-Za-z]+[+\-*]?[0-9]?$/.test(tok)) return true;
  // X_s, X_d type with star ("K*", "K*+", "K*0")
  if (/^[A-Z][A-Za-z]?\*?[+\-]?[0-9]?$/.test(tok)) return true;
  return false;
}

function normalizeShorthandToken(tok: string): string {
  let out = tok
    .replace(/(?<!\\)\bnubar\b/g, '\\bar{\\nu}')
    .replace(/(?<!\\)\bbbar\b/g, '\\bar{b}')
    .replace(/(?<!\\)\bcbar\b/g, '\\bar{c}')
    .replace(/(?<!\\)\btbar\b/g, '\\bar{t}')
    .replace(/(?<!\\)\bsbar\b/g, '\\bar{s}')
    .replace(/(?<!\\)\bdbar\b/g, '\\bar{d}')
    // Charged-pi style: pi+, pi-, pi0 -> \pi^+, \pi^-, \pi^0.
    .replace(/(?<!\\)\bpi([+\-])/g, '\\pi^$1')
    .replace(/(?<!\\)\bpi0\b/g, '\\pi^0')
    // Lepton charges: mu+, mu-, tau+, e+ -> \mu^+, etc.
    .replace(/(?<!\\)\bmu([+\-])/g, '\\mu^$1')
    .replace(/(?<!\\)\btau([+\-])/g, '\\tau^$1')
    .replace(/(?<!\\)\bell([+\-])/g, '\\ell^$1')
    .replace(/(?<!\\)\be([+\-])(?=\W|$)/g, 'e^$1')
    // Lowercase Greek: allow trailing "_" as boundary (e.g. "nu_e", "mu_nu").
    .replace(/(?<![\\A-Za-z])nu(?=\b|_)/g, '\\nu')
    .replace(/(?<![\\A-Za-z])mu(?=\b|_)/g, '\\mu')
    .replace(/(?<![\\A-Za-z])tau(?=\b|_)/g, '\\tau')
    .replace(/(?<![\\A-Za-z])gamma(?=\b|_)/g, '\\gamma')
    .replace(/(?<![\\A-Za-z])pi(?=\b|_)/g, '\\pi')
    .replace(/(?<![\\A-Za-z])phi(?=\b|_)/g, '\\phi')
    .replace(/(?<![\\A-Za-z])psi(?=\b|_)/g, '\\psi')
    .replace(/(?<![\\A-Za-z])ell(?=\b|_)/g, '\\ell')
    .replace(/(?<![\\A-Za-z])rho(?=\b|_)/g, '\\rho')
    .replace(/(?<![\\A-Za-z])omega(?=\b|_)/g, '\\omega')
    .replace(/(?<![\\A-Za-z])eta(?=\b|_)/g, '\\eta')
    // Capital Greek used in baryons and FCNC observables.
    .replace(/(?<![\\A-Za-z])Lambda(?=\b|_)/g, '\\Lambda')
    .replace(/(?<![\\A-Za-z])Sigma(?=\b|_)/g, '\\Sigma')
    .replace(/(?<![\\A-Za-z])Delta(?=\b|_)/g, '\\Delta')
    .replace(/(?<![\\A-Za-z])Gamma(?=\b|_)/g, '\\Gamma')
    .replace(/(?<![\\A-Za-z])Omega(?=\b|_)/g, '\\Omega')
    .replace(/(?<![\\A-Za-z])Xi(?=\b|_)/g, '\\Xi')
    .replace(/(?<![\\A-Za-z])Upsilon(?=\b|_)/g, '\\Upsilon')
    .replace(/(?<![\\A-Za-z])Phi(?=\b|_)/g, '\\Phi')
    .replace(/(?<![\\A-Za-z])Psi(?=\b|_)/g, '\\Psi');
  // Wrap multi-letter subscripts/superscripts in braces so KaTeX doesn't
  // truncate (matches the same fixup applied in lib/notation.ts).
  out = out.replace(/([_^])\(([^()]*)\)/g, (_m, op, body) => `${op}{(${body})}`);
  out = out.replace(/([_^])([A-Za-z][A-Za-z]+)(?![A-Za-z{}])/g, (_m, op, idx) => `${op}{${idx}}`);
  return out;
}

function wrapPhysicsShorthand(text: string): string {
  // Run all patterns, then iterate once more so trailing-Greek merges that
  // produce new `$..$` spans can themselves absorb the next trailing token.
  // Two passes are sufficient for the catalog content; bound the loop at 4
  // to be safe.
  let out = text;
  for (let pass = 0; pass < 4; pass++) {
    const before = out;
    for (const [pat, replacer] of SHORTHAND_PATTERNS) {
      out = out.replace(pat, (...args) => {
        const exec = args.slice(0, -2);
        return replacer(exec as unknown as RegExpExecArray);
      });
    }
    if (out === before) break;
  }
  return out;
}
