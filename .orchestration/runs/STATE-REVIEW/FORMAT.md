# STATE REVIEW — shared format + rules (read first)

Goal: a complete, honest, code-grounded review of the project's CURRENT state so the
PI can see exactly what is rigorous, what is approximate, what choices were made, and
what remains. This feeds a master STATE_OF_PROJECT.md. Be thorough and precise; the PI
asked for "complete and total ... very thoroughly and closely."

## Grounding rules (NON-NEGOTIABLE)
- Read the ACTUAL CODE, not just the ledgers. Ledgers/memory can be stale — verify every
  claim against the real files and cite `path:line`.
- Distinguish what is IMPLEMENTED from what was actually SCANNED (the only scans run are
  quark-only: minimal `scan_outputs/wq_quarkonly_1M_20128400` and custodial
  `scan_outputs/wq_quarkonly_1M_custodial_20675555`; lepton sector was dropped in both).
- Be honest about approximations and untested paths. Do NOT call something rigorous if it
  is a proxy, an ideal-limit choice, a deferred term, or a benchmark stand-in.
- Where a number/convention matters, state it explicitly (e.g. M_KK physical = 2.45*Lambda_IR;
  floor threshold; MFV r weight; LFV C=0.02 diagnostic vs BR veto; epsilon_K normalization).

## Output format — write your section to the assigned file with EXACTLY these four buckets:
### IN (rigorous, fully implemented + tested; note if also SCANNED)
  - one bullet per item: what it is, the constraint/observable, file:line, test ref, scanned? (yes/no)
### APPROXIMATE (proxy / deferred / ideal-limit / benchmark stand-in)
  - one bullet per item: what it is, EXACTLY what approximation is made and why, file:line,
    how far it could be from rigorous, what would make it rigorous.
### CHOICES MADE (conventions, model decisions, defaults — with rationale + where set)
  - one bullet per item: the choice, the value/default, where it is set (file:line), the rationale,
    and any alternative that was considered/rejected.
### TODO / OPEN (what remains, and whether it needs human physics input vs is buildable)
  - one bullet per item: the work, blocker (human-input vs build), rough size.

End your file with `SECTION-READY`. Cite sources. Canonical docs to cross-check (but verify
against code): `.orchestration/PHASE2_PROGRAM_LEDGER.md`, `.orchestration/NEEDS_HUMAN_PHYSICS.md`,
`.orchestration/REBUILD_LEDGER.md` (if present), `docs/*.md`, `CLAUDE.md`, and the memory file
referenced in the ledger.
