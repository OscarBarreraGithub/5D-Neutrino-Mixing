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
