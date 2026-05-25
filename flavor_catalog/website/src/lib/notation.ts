const ALREADY_LATEX_PATTERN =
  /\\(?:to|bar|mathrm|mathcal|varepsilon|pi|nu|mu|ell|rightarrow|frac)|\\[A-Za-z]/;

const BAR_REPLACEMENTS: Array<[RegExp, string]> = [
  [/\btibar(?=\b|_)/g, '\\bar{D}'],
  [/\bDbar(?=\b|_)/g, '\\bar{D}'],
  [/\bBbar(?=\b|_)/g, '\\bar{B}'],
  [/\bTbar(?=\b|_)/g, '\\bar{T}'],
  [/\bCbar(?=\b|_)/g, '\\bar{C}'],
  [/\bKbar(?=\b|_)/g, '\\bar{K}'],
  [/\bnubar(?=\b|_)/g, '\\bar{\\nu}'],
  [/\bbbar(?=\b|_)/gi, '\\bar{b}'],
  [/\btbar(?=\b|_)/g, '\\bar{t}'],
  [/\bcbar(?=\b|_)/g, '\\bar{c}'],
  [/\bsbar(?=\b|_)/g, '\\bar{s}'],
  [/\bdbar(?=\b|_)/g, '\\bar{d}'],
  [/\bubar(?=\b|_)/g, '\\bar{u}'],
];

// Match a Greek-letter shorthand even when the surrounding boundary is "_" (a
// word character in JS regex) -- so "nu_e", "mu_nu", "pi_0", "tau_e", and
// "_mu" (as a subscript fragment) all get the Greek replacement.  Without the
// "_" relaxation, those tokens rendered the prefix as italic Latin letters
// (n, m, p, t) with a single-character subscript glued on.
const GREEK_REPLACEMENTS: Array<[RegExp, string]> = [
  [/(?<![\\A-Za-z])pi(?=\b|_)/g, '\\pi'],
  [/(?<![\\A-Za-z])nu(?=\b|_)/g, '\\nu'],
  [/(?<![\\A-Za-z])mu(?=\b|_)/g, '\\mu'],
  [/(?<![\\A-Za-z])ell(?=\b|_)/g, '\\ell'],
  [/(?<![\\A-Za-z])gamma(?=\b|_)/g, '\\gamma'],
  [/(?<![\\A-Za-z])tau(?=\b|_)/g, '\\tau'],
  [/(?<![\\A-Za-z])phi(?=\b|_)/g, '\\phi'],
  [/(?<![\\A-Za-z])psi(?=\b|_)/g, '\\psi'],
  [/(?<![\\A-Za-z])upsilon(?=\b|_)/g, '\\upsilon'],
  [/(?<![\\A-Za-z])eta(?=\b|_)/g, '\\eta'],
  [/(?<![\\A-Za-z])rho(?=\b|_)/g, '\\rho'],
  [/(?<![\\A-Za-z])sigma(?=\b|_)/g, '\\sigma'],
  [/(?<![\\A-Za-z])omega(?=\b|_)/g, '\\omega'],
  [/(?<![\\A-Za-z])lambda(?=\b|_)/g, '\\lambda'],
  [/(?<![\\A-Za-z])xi(?=\b|_)/g, '\\xi'],
  [/(?<![\\A-Za-z])chi(?=\b|_)/g, '\\chi'],
  [/(?<![\\A-Za-z])theta(?=\b|_)/g, '\\theta'],
  [/(?<![\\A-Za-z])epsilon(?=\b|_)/g, '\\epsilon'],
  [/(?<![\\A-Za-z])varepsilon(?=\b|_)/g, '\\varepsilon'],
  [/(?<![\\A-Za-z])beta(?=\b|_)/g, '\\beta'],
  [/(?<![\\A-Za-z])alpha(?=\b|_)/g, '\\alpha'],
  [/(?<![\\A-Za-z])delta(?=\b|_)/g, '\\delta'],
  [/(?<![\\A-Za-z])kappa(?=\b|_)/g, '\\kappa'],
  // Capital Greek letters that appear in physics shorthand without
  // backslashes (PDG-table observables like "Delta m_D", "Gamma_D",
  // "Sigma", etc).
  [/(?<![\\A-Za-z])Upsilon(?=\b|_)/g, '\\Upsilon'],
  [/(?<![\\A-Za-z])Delta(?=\b|_)/g, '\\Delta'],
  [/(?<![\\A-Za-z])Gamma(?=\b|_)/g, '\\Gamma'],
  [/(?<![\\A-Za-z])Sigma(?=\b|_)/g, '\\Sigma'],
  [/(?<![\\A-Za-z])Omega(?=\b|_)/g, '\\Omega'],
  [/(?<![\\A-Za-z])Lambda(?=\b|_)/g, '\\Lambda'],
  [/(?<![\\A-Za-z])Phi(?=\b|_)/g, '\\Phi'],
  [/(?<![\\A-Za-z])Psi(?=\b|_)/g, '\\Psi'],
  [/(?<![\\A-Za-z])Theta(?=\b|_)/g, '\\Theta'],
  [/(?<![\\A-Za-z])Xi(?=\b|_)/g, '\\Xi'],
];

interface TextPhrase {
  phrase: string;
  replacement: string;
  trimAdjacentSpaces?: boolean;
}

const TEXT_PHRASES: TextPhrase[] = [
  {
    phrase: 'pole observables:',
    replacement: '\\,\\text{pole observables: }\\,',
    trimAdjacentSpaces: true,
  },
  {
    phrase: 'observables:',
    replacement: '\\,\\text{observables: }\\,',
    trimAdjacentSpaces: true,
  },
  { phrase: 'experimentally', replacement: '\\,\\text{experimentally }\\,' },
  {
    phrase: '(high-mass tail, EFT contact-operator)',
    replacement: '\\,\\text{(high-mass tail, EFT contact-operator)}\\,',
  },
  {
    phrase: '(high-mass spin-0/spin-2 resonance)',
    replacement: '\\,\\text{(high-mass spin-0/spin-2 resonance)}\\,',
  },
  {
    phrase: 'high-mass spin-0/spin-2 resonance',
    replacement: '\\,\\text{high-mass spin-0/spin-2 resonance}\\,',
  },
  { phrase: '(singlet)', replacement: '\\,\\text{(singlet)}\\,' },
  { phrase: 'doublet', replacement: '\\,\\text{doublet}\\,' },
  { phrase: 'mixed final states', replacement: '\\,\\text{mixed final states}\\,' },
];

// Prose words that frequently appear embedded inside otherwise-mathematical
// strings (e.g. "D^0-\bar D^0 mixing", "|V_cb|_inclusive vs ...") and would
// otherwise render as italic single-letter math (m-i-x-i-n-g, i-n-c-l-...).
// These are applied as word-boundary regex replacements so we never match an
// "inclusive" inside another token, and the leading/trailing thin-space keeps
// the surrounding math from colliding with the upright text.
const TEXT_WORD_REPLACEMENTS: Array<[RegExp, string]> = [
  [/\bmixing\b/g, '\\,\\text{mixing}\\,'],
  [/\binclusive\b/g, '\\,\\text{inclusive}\\,'],
  [/\bexclusive\b/g, '\\,\\text{exclusive}\\,'],
  [/\bvs\b/g, '\\,\\text{vs}\\,'],
  // " in " (with surrounding whitespace) only -- prevents matching "in" as
  // part of words like "spin", "string", "indirect".
  [/(?<=\s)in(?=\s)/g, '\\,\\text{in}\\,'],
];

function applySharedRules(value: string, isAlreadyLatex: boolean): string {
  let normalized = value.replace(/\s*(?<![<-])->\s*/g, ' \\to ');
  normalized = normalized.replace(/\s*(?:<-->|<->)\s*/g, ' \\leftrightarrow ');

  if (!isAlreadyLatex || !normalized.includes('\\gamma')) {
    normalized = normalized.replace(/\bgamma gamma\b/g, '\\gamma\\gamma');
  }

  return normalized;
}

function applyBarReplacements(value: string): string {
  return BAR_REPLACEMENTS.reduce(
    (current, [pattern, replacement]) => current.replace(pattern, replacement),
    value,
  );
}

function applyGreekReplacements(value: string): string {
  return GREEK_REPLACEMENTS.reduce(
    (current, [pattern, replacement]) => current.replace(pattern, replacement),
    value,
  );
}

function applyChargeSuperscripts(value: string): string {
  return value
    .replace(/(\\(?:mu|pi|tau|nu|gamma|ell))([+-])/g, '$1^$2')
    .replace(/\b(e|p|n|u|d|s|c|b|t|l|W)([+-])(?=\s|$|[,;)])/g, '$1^$2');
}

// Wrap multi-character subscripts and superscripts in braces so KaTeX does not
// truncate them to the first character.  Examples:
//   d_Hg     -> d_{Hg}
//   d_Ra     -> d_{Ra}
//   d_Xe     -> d_{Xe}
//   |V_cb|   -> |V_{cb}|
//   A_FB     -> A_{FB}
//   G^(1)    -> G^{(1)}    (KK level indices that authors write as ^(n))
// Already-braced groups ("_{...}"), single-character indices ("_b", "x_D"),
// and macro-style indices ("_\mu", "_\bar{b}") are left untouched.
function applyMultiCharSubAndSup(value: string): string {
  let out = value;
  // Parenthesised super/sub-scripts: ^(...) -> ^{(...)}, _(...) -> _{(...)}.
  // Without the braces KaTeX puts only "(" as the index and dumps the rest
  // (digits and closing paren) outside.
  out = out.replace(/([_^])\(([^()]*)\)/g, (_m, op, body) => `${op}{(${body})}`);
  // Multi-letter alphabetic indices: _<2+ letters> or ^<2+ letters>.
  out = out.replace(/([_^])([A-Za-z][A-Za-z]+)(?![A-Za-z{}])/g, (_m, op, idx) => `${op}{${idx}}`);
  return out;
}

function escapeRegExp(value: string): string {
  return value.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}

function applyTextWrappers(value: string): string {
  const replacements: string[] = [];
  let normalized = value;

  for (const { phrase, replacement, trimAdjacentSpaces } of TEXT_PHRASES) {
    const replacer = () => {
      const placeholder = `@@NOTATION_TEXT_${replacements.length}@@`;
      replacements.push(replacement);
      return placeholder;
    };

    if (trimAdjacentSpaces) {
      normalized = normalized.replace(new RegExp(`\\s*${escapeRegExp(phrase)}\\s*`, 'g'), replacer);
    } else {
      normalized = normalized.replaceAll(phrase, replacer);
    }
  }

  for (const [pattern, replacement] of TEXT_WORD_REPLACEMENTS) {
    normalized = normalized.replace(pattern, () => {
      const placeholder = `@@NOTATION_TEXT_${replacements.length}@@`;
      replacements.push(replacement);
      return placeholder;
    });
  }

  normalized = normalized.replace(/\band\b/g, () => {
    const placeholder = `@@NOTATION_TEXT_${replacements.length}@@`;
    replacements.push('\\,\\text{ and }\\,');
    return placeholder;
  });

  replacements.forEach((replacement, index) => {
    normalized = normalized.replaceAll(`@@NOTATION_TEXT_${index}@@`, replacement);
  });

  return normalized;
}

export function normalizeNotation(raw: string): string {
  const trimmed = raw.trim();
  const isAlreadyLatex = ALREADY_LATEX_PATTERN.test(trimmed);
  let normalized = applySharedRules(trimmed, isAlreadyLatex);

  if (!isAlreadyLatex) {
    normalized = applyBarReplacements(normalized);
    normalized = applyGreekReplacements(normalized);
    normalized = applyChargeSuperscripts(normalized);
  }

  // Run the text-phrase and multi-character subscript/superscript fixups for
  // every string, including ones that already contain LaTeX commands.  The
  // text phrases (mixing, in, vs, ...) are prose words that authors sometimes
  // embed in otherwise-math strings; the multi-char subscript pass corrects
  // unbraced indices such as d_Hg or |V_cb|.
  normalized = applyMultiCharSubAndSup(normalized);
  normalized = applyTextWrappers(normalized);

  return normalized;
}

export function notationToMath(raw: string): string {
  return `$${normalizeNotation(raw)}$`;
}

/**
 * processNameToHtml(raw)
 *
 * Format a process_name caption that mixes English prose with physics
 * shorthand (e.g. "B -> phi K_S penguin CP asymmetry").  Long stretches of
 * English are left as text; embedded physics tokens (arrows, subscripts,
 * particle names, bars) are normalised and wrapped in $...$ so KaTeX renders
 * them inline.  This is intentionally conservative: a string with no detected
 * physics tokens is returned unchanged.
 */
const PROSE_KEEP_WORDS = new Set([
  'and', 'or', 'in', 'of', 'the', 'a', 'an', 'to', 'with', 'over',
  'CP', 'EW', 'LFV', 'LFU', 'KK', 'NP', 'BSM', 'SM',
  'rare', 'inclusive', 'exclusive', 'radiative', 'leptonic', 'semileptonic',
  'baryonic', 'nonleptonic', 'invisible', 'electronic', 'tauonic', 'muonic',
  'charged', 'charmless', 'penguin', 'mixing', 'oscillation', 'frequency',
  'mass', 'splitting', 'difference', 'phase', 'asymmetry', 'decay', 'decays',
  'ratio', 'production', 'pair', 'resonance', 'partner', 'partners',
  'universality', 'flavor', 'flavour', 'gauge', 'electroweak', 'top',
  'bottom', 'down', 'up', 'strange', 'charm', 'beauty', 'kaon', 'neutral',
  'lepton', 'lepton-flavor-violating', 'lepton-flavor-universality',
  'fraction', 'branching', 'limit', 'bound', 'measurement', 'dipole',
  'high-mass', 'forward-backward', 'high', 'low', 'tail', 'long', 'short',
  'angular', 'polarization', 'helicity', 'angle', 'observable', 'observables',
  'integrated', 'modes', 'mode', 'channel', 'channels', 'state', 'states',
  'process', 'processes', 'transition', 'transitions', 'spectator', 'penguin-only',
  'is', 'are', 'be', 'been', 'we', 'this', 'that', 'these', 'those',
  'precision', 'oscillation', 'dependence', 'family',
]);

// Heuristic: this token looks like physics shorthand (uses subscripts,
// superscripts, particle letters paired with bar/_, etc.).  Used to decide
// whether to wrap a run of tokens in inline math.
function looksLikePhysics(token: string): boolean {
  if (!token) return false;
  // ASCII-arrow, fraction, particle symbols.
  if (/->|<->|\^|_/.test(token)) return true;
  if (/\b(?:bar|nubar|bbar|cbar|tbar|tibar|sbar|dbar|ubar|Dbar|Bbar|Tbar|Cbar|Kbar)\b/.test(token)) return true;
  if (/\b(?:pi|mu|nu|ell|gamma|tau|phi|psi|rho|omega|eta|chi|theta|kappa|sigma|lambda|xi|epsilon|varepsilon|Upsilon|Delta|Gamma|Sigma|Lambda|Omega|Phi|Psi|Theta|Xi)\b/.test(token)) return true;
  // Single capital letter followed by a charge sign / asterisk / digit at the
  // end of the token (e.g. "K+", "K*", "W-", "B0", "K*+").  Reject compounds
  // like "Z-pole" or "B-meson" where prose follows the dash.
  if (/^[A-Z][A-Za-z]?[+\-*][0-9]*$/.test(token)) return true;
  if (/^[A-Z][_^]/.test(token)) return true;
  // Pure mathematical fragments like "T_{5/3}" or "K_S".
  if (/\{|\}/.test(token)) return true;
  return false;
}

/**
 * Plain-text variant of processNameToHtml for places where math cannot
 * render (e.g. the <title> tag).  Replaces ASCII arrows and bar shorthand
 * with Unicode equivalents so the tab text reads naturally.
 */
export function processNameToPlain(raw: string): string {
  if (!raw) return '';
  return raw
    .replace(/(?<![\-<])->(?![>\-])/g, '→')
    .replace(/<->/g, '↔')
    .replace(/\bnubar\b/g, 'ν̄')
    .replace(/\bbbar\b/g, 'b̄')
    .replace(/\bcbar\b/g, 'c̄')
    .replace(/\btbar\b/g, 't̄')
    .replace(/\bsbar\b/g, 's̄')
    .replace(/\bdbar\b/g, 'd̄')
    .replace(/\bubar\b/g, 'ū')
    .replace(/\bBbar\b/g, 'B̄')
    .replace(/\bDbar\b/g, 'D̄')
    .replace(/\bKbar\b/g, 'K̄')
    .replace(/\bTbar\b/g, 'T̄')
    .replace(/\bCbar\b/g, 'C̄');
}

export function processNameToHtml(raw: string): string {
  if (!raw) return '';
  const trimmed = raw.trim();
  if (!trimmed) return '';
  // If author already wrapped the whole thing in math, defer to KaTeX.
  if (trimmed.startsWith('$') && trimmed.endsWith('$')) return trimmed;

  // Token grouping: split on whitespace but keep punctuation attached.
  const tokens = trimmed.split(/(\s+)/);
  const out: string[] = [];
  let buffer: string[] = [];
  let bufferTrailingSpace = '';

  function flushBuffer() {
    if (buffer.length === 0) {
      if (bufferTrailingSpace) {
        out.push(bufferTrailingSpace);
        bufferTrailingSpace = '';
      }
      return;
    }
    // Join, normalise (handles ->, bar, greek), wrap in math.  Drop only
    // surrounding whitespace that we collected as separators.
    const joined = buffer.join('').replace(/^\s+|\s+$/g, '');
    if (joined) {
      out.push(`$${normalizeNotation(joined)}$`);
    }
    buffer = [];
    if (bufferTrailingSpace) {
      out.push(bufferTrailingSpace);
      bufferTrailingSpace = '';
    }
  }

  for (const tok of tokens) {
    if (/^\s+$/.test(tok)) {
      if (buffer.length > 0) {
        // Save the trailing whitespace; we'll decide whether to keep it inside
        // the math span or emit it after the closing $ once we see the next
        // token.
        bufferTrailingSpace = tok;
      } else {
        out.push(tok);
      }
      continue;
    }
    // Strip trailing punctuation to test the bare word.
    const stripped = tok.replace(/[.,;:!?)]+$/u, '').replace(/^[(]/u, '');
    const tail = tok.slice(stripped.length + (tok.startsWith('(') ? 1 : 0));
    const head = tok.startsWith('(') ? '(' : '';
    const isPhysics = looksLikePhysics(stripped) || /[\\]/.test(stripped);
    const isProse = PROSE_KEEP_WORDS.has(stripped) || PROSE_KEEP_WORDS.has(stripped.toLowerCase());

    if (isPhysics && !isProse) {
      // Continue the math run: incorporate any pending whitespace into the math span.
      if (bufferTrailingSpace) {
        buffer.push(bufferTrailingSpace);
        bufferTrailingSpace = '';
      }
      buffer.push(head, stripped);
      // Trailing punctuation goes outside the math (so we don't render it inside $).
      if (tail) {
        flushBuffer();
        out.push(tail);
      }
    } else {
      flushBuffer();
      out.push(tok);
    }
  }
  flushBuffer();
  return out.join('');
}
