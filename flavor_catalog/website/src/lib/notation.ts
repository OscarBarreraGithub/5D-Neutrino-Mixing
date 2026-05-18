const ALREADY_LATEX_PATTERN =
  /\\(?:to|bar|mathrm|mathcal|varepsilon|pi|nu|mu|ell|rightarrow|frac)|\\[A-Za-z]/;

const BAR_REPLACEMENTS: Array<[RegExp, string]> = [
  [/\btibar(?=\b|_)/g, '\\bar{D}'],
  [/\bDbar(?=\b|_)/g, '\\bar{D}'],
  [/\bBbar(?=\b|_)/g, '\\bar{B}'],
  [/\bTbar(?=\b|_)/g, '\\bar{T}'],
  [/\bCbar(?=\b|_)/g, '\\bar{C}'],
  [/\bnubar(?=\b|_)/g, '\\bar{\\nu}'],
  [/\bbbar(?=\b|_)/gi, '\\bar{b}'],
  [/\btbar(?=\b|_)/g, '\\bar{t}'],
  [/\bcbar(?=\b|_)/g, '\\bar{c}'],
];

const GREEK_REPLACEMENTS: Array<[RegExp, string]> = [
  [/(?<!\\)\bpi\b/g, '\\pi'],
  [/(?<!\\)\bnu\b/g, '\\nu'],
  [/(?<!\\)\bmu\b/g, '\\mu'],
  [/(?<!\\)\bell\b/g, '\\ell'],
  [/(?<!\\)\bgamma\b/g, '\\gamma'],
  [/(?<!\\)\btau\b/g, '\\tau'],
  [/(?<!\\)\bphi\b/g, '\\phi'],
  [/(?<!\\)\bpsi\b/g, '\\psi'],
  [/(?<!\\)\bUpsilon\b/g, '\\Upsilon'],
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

  if (isAlreadyLatex) {
    return normalized;
  }

  normalized = applyBarReplacements(normalized);
  normalized = applyGreekReplacements(normalized);
  normalized = applyChargeSuperscripts(normalized);
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
  if (/\b(?:bar|nubar|bbar|cbar|tbar|tibar|Dbar|Bbar|Tbar|Cbar)\b/.test(token)) return true;
  if (/\b(?:pi|mu|nu|ell|gamma|tau|phi|psi|Upsilon|Delta|epsilon)\b/.test(token)) return true;
  // Single capital letter followed by digits, plus/minus, or asterisks.
  if (/^[A-Z][A-Za-z]?[+\-*]/.test(token)) return true;
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
    .replace(/\bBbar\b/g, 'B̄')
    .replace(/\bDbar\b/g, 'D̄');
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
