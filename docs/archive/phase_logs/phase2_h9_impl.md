# Phase 2 Hole #9 Implementation

## Inventory summary

Created `docs/audits/scope_wording_inventory.md` with the required three-column
inventory: topic, current wording, and recommended fix. The inventory records
the methodology note's first-KK wording, the stale compact-note first-mode TODO,
the external EWPO-floor references, the spread of `g_s^*` conventions, and the
root/decision-log statements that make the repo broader than the present
quark-sector paper.

## New appendix subsection summary

KK tower truncation: added `Scope and approximations` to
`docs/quark_scan_methodology_note.tex`, stating that the scan includes only the
lightest colour-octet KK gluon (`n=1`). The subsection records the literature
full-tower estimate as an additional `+20-30%` systematic on the epsilon_K
ratio, separate from the BGS budget band and the LO BMU running uncertainty.

EWPO floor: added a paragraph that treats `M_KK >= 3 TeV` from electroweak
precision as an externally cited custodial-RS floor, not a result computed by
this repository. The same paragraph states that the quark-flavour p50 bound,
`47.26 TeV` at `g_s^*=3`, is independent and substantially stronger.

`g_s^*` convention: added a paragraph making `g_s^*=3` the paper headline
convention. The perturbative `g_s ~= 1.05` values are identified as the
code-default cross-reference, while `g_s^* ~= 6` is reserved for the CFW
boundary-term comparison only.

Lepton-sector scope: added a paragraph stating that the repository hosts both
quark and lepton 5D-flavour tools, but this paper is quark-sector only. The
`neutrinos/`, `flavorConstraints/`, and `yukawa/` directories are follow-up
work and were not audited or rerun under the present `deltaf2.py` updates.

## Inline body edits summary

The main text now cross-references the scope appendix from the reading guide and
from the first KK-gluon discussion. The result section now presents the
`g_s^*=3` values as the headline crossings and the perturbative values as the
code-default comparison. The previous inline NLO/systematics placeholder was
reduced to a cross-reference to the appendix, and the recommendation section now
points to the central BGS/FLAG inputs and scope conventions instead of
repeating loose "LO-only" wording.

## Verification

- Rebuilt `docs/quark_scan_methodology_note.pdf` with two `pdflatex
  -interaction=nonstopmode quark_scan_methodology_note.tex` passes from
  `docs/`.
- Confirmed the regenerated PDF has 19 pages with `pdfinfo`.
- Checked the build log for undefined references; none remain after the second
  pass.

Hole #9 ready for peer review.
