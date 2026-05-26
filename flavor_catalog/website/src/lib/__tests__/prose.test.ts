/**
 * R22-I1 — snapshot suite for ``src/lib/prose.ts``.
 *
 * Pins ``texProseToHtml`` (LaTeX text-mode -> HTML, em-dash normalization,
 * jargon scrub, math-shorthand wrapping) for a small corpus of catalog
 * prose strings drawn from CLEANUP_PLAN.md §C C04.  We exercise both the
 * raw ``process_name`` fields (mostly plain prose) and short prose
 * fragments that contain the LaTeX commands the helper rewrites.
 *
 * Snapshots live in ``./__snapshots__/prose.test.ts.snap``.  Regenerate
 * with ``npm run test:update`` after a deliberate rule change.
 */
import { describe, it, expect } from 'vitest';
import { texProseToHtml } from '../prose.ts';

interface ProseCase {
  id: string;
  raw: string;
}

const PROCESS_NAME_CORPUS: readonly ProseCase[] = [
  { id: 'T010', raw: 'Z-pole bottom-quark partial-width ratio and asymmetries' },
  { id: 'CR002', raw: 'Pair production of exotic-charge T_{5/3} custodial top partner' },
  { id: 'K018', raw: 'Kaon semileptonic decays and source-level V_us input' },
  { id: 'B015', raw: 'Inclusive b -> s ell+ ell- rare B decay' },
  { id: 'K020', raw: 'LFV charged-kaon semileptonic decay' },
];

// A small set of prose-flavored strings that exercise the LaTeX text-mode
// rewriting + jargon scrub + shorthand wrapping paths.  These keep the
// snapshot suite useful even though ``process_name`` fields are mostly
// plain English.
const PROSE_CORPUS: readonly ProseCase[] = [
  {
    id: 'texttt-textbf-mix',
    raw: 'The \\texttt{deltaf2} module wraps \\textbf{Delta F = 2} Wilson coefficients.',
  },
  {
    id: 'em-dash-shorthand',
    raw: 'Inclusive b -> s gamma — the rate is dominated by the radiative penguin.',
  },
  {
    id: 'href-and-jargon',
    raw: 'See \\href{https://example.org/note}{the PKA note} for Wave-3 results.',
  },
  {
    id: 'shorthand-arrows-and-pairs',
    raw: 'The Z b bbar coupling and nu nubar pair-production constrain B -> K nubar.',
  },
];

describe('texProseToHtml — process_name corpus', () => {
  for (const { id, raw } of PROCESS_NAME_CORPUS) {
    it(`${id}`, () => {
      expect(texProseToHtml(raw)).toMatchSnapshot();
    });
  }
});

describe('texProseToHtml — prose fragments exercising rewrite rules', () => {
  for (const { id, raw } of PROSE_CORPUS) {
    it(`${id}`, () => {
      expect(texProseToHtml(raw)).toMatchSnapshot();
    });
  }
});
