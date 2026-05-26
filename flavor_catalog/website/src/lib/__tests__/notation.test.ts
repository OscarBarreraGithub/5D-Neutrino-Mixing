/**
 * R22-I1 — snapshot suite for ``src/lib/notation.ts``.
 *
 * Pins ``normalizeNotation``, ``notationToMath``, ``processNameToHtml``,
 * and ``processNameToPlain`` outputs for a small corpus of raw catalog
 * ``process_name`` strings drawn from CLEANUP_PLAN.md §C C04:
 *
 *   * T010 — Z-pole bottom-quark partial-width ratio and asymmetries
 *   * CR002 — Pair production of exotic-charge T_{5/3} custodial top partner
 *   * K018 — Kaon semileptonic decays and source-level V_us input
 *   * B015 — Inclusive b -> s ell+ ell- rare B decay
 *   * K020 — LFV charged-kaon semileptonic decay
 *
 * Snapshots live in ``./__snapshots__/notation.test.ts.snap``.  Regenerate
 * with ``npm run test:update`` after a deliberate notation-rule change.
 */
import { describe, it, expect } from 'vitest';
import {
  normalizeNotation,
  notationToMath,
  processNameToHtml,
  processNameToPlain,
} from '../notation.ts';

interface Entry {
  id: string;
  // Raw ``process_name`` field as committed in the catalog yaml.
  raw: string;
}

const CORPUS: readonly Entry[] = [
  { id: 'T010', raw: 'Z-pole bottom-quark partial-width ratio and asymmetries' },
  { id: 'CR002', raw: 'Pair production of exotic-charge T_{5/3} custodial top partner' },
  { id: 'K018', raw: 'Kaon semileptonic decays and source-level V_us input' },
  { id: 'B015', raw: 'Inclusive b -> s ell+ ell- rare B decay' },
  { id: 'K020', raw: 'LFV charged-kaon semileptonic decay' },
];

describe('notation helpers — catalog process_name corpus', () => {
  for (const { id, raw } of CORPUS) {
    describe(id, () => {
      it('normalizeNotation', () => {
        expect(normalizeNotation(raw)).toMatchSnapshot();
      });

      it('notationToMath', () => {
        expect(notationToMath(raw)).toMatchSnapshot();
      });

      it('processNameToHtml', () => {
        expect(processNameToHtml(raw)).toMatchSnapshot();
      });

      it('processNameToPlain', () => {
        expect(processNameToPlain(raw)).toMatchSnapshot();
      });
    });
  }
});
