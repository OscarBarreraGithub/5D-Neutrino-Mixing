# Writer Worklog: wa_w4_charged_lepton

Family: `charged_lepton`
Processes: `L003`, `L004`, `L008`, `L009`, `L010`
Writer agent: `WA`
Timestamp: `2026-05-16T13:31:05-04:00`

## L003 -- \(\mu-e\) conversion in aluminum

- Normalized Section B headings and tightened the distinction between the
  published SINDRUM II gold benchmark and aluminum-only Mu2e/COMET prospects.
- Added sidecar field references next to the gold bound and each aluminum
  projection retained in the prose.
- Appended `WRITER-DONE` to `L003.yaml` and updated `last_updated_at`.

Open issues for CA:
- Verify whether a current PDG LFV summary has a dedicated muon-conversion
  table newer than the experiment/equivalent snapshots used by PKA.
- Check that no direct published aluminum-target limit should replace the
  prospective treatment.

## L004 -- \(\mu-e\) conversion in gold

- Normalized the PDG-or-equivalent and validity headings.
- Added sidecar references for the PDG live generation/access dates, SINDRUM II
  \(7\times10^{-13}\) limit, Crossref publication metadata, and aluminum-program
  comparison numbers.
- Appended `WRITER-DONE` to `L004.yaml` and updated `last_updated_at`.

Open issues for CA:
- Confirm that the PDG live NODE=S004 stream remains the preferred canonical
  source over a static annual listing snapshot.
- Confirm the operator/nuclear-overlap caveat is strong enough before any
  direct Au/Al/Ti comparison is made.

## L008 -- \(\tau\to e\gamma\)

- Normalized the entry to the plan-v1 Section B field names and tightened the
  PDG/BaBar/Belle/Belle II value hierarchy.
- Added sidecar references for the \(3.3\times10^{-8}\) PDG/BaBar bound,
  Belle \(5.6\times10^{-8}\) cross-check, and Belle II \(50\,{\rm ab}^{-1}\)
  prospect.
- Appended `WRITER-DONE` to `L008.yaml` and updated `last_updated_at`.

Open issues for CA:
- Confirm whether any Belle II \(\tau\to e\gamma\) measurement supersedes the
  BaBar/PDG value; PKA found only Belle II prospects.
- Decide later whether implementation should generalize `muToEGamma.py` or use
  a separate tau radiative-LFV module.

## L009 -- \(\tau\to3\mu\)

- Added explicit standard notation and normalized Section B headings.
- Added sidecar references for the PDG/Belle II \(1.9\times10^{-8}\) limit,
  the LHCb Run 2 \(1.9\,(2.3)\times10^{-8}\) preprint, and historical
  Belle/BaBar/LHCb limits.
- Appended `WRITER-DONE` to `L009.yaml` and updated `last_updated_at`.

Open issues for CA:
- Decide whether the 2026 LHCb preprint should be treated as co-primary with
  Belle II before PDG incorporates it.
- Confirm the future implementation basis should include dipole/contact
  interference rather than a pure contact-only bound.

## L010 -- \(\tau\to3e\)

- Added explicit standard notation and normalized Section B headings.
- Added sidecar references for the PDG/Belle \(2.7\times10^{-8}\) limit,
  BaBar \(2.9\times10^{-8}\) cross-check, LHCb context exclusion, and Belle II
  three-lepton prospects.
- Appended `WRITER-DONE` to `L010.yaml` and updated `last_updated_at`.

Open issues for CA:
- Independently check whether a post-2010 Belle II, LHCb, ATLAS, or CMS
  electron-mode result supersedes the PDG-selected Belle value.
- Confirm that LHCb is only contextual for L010 and not a numerical input.

## Batch notes

- No `\textbf{CHECK}` markers were introduced.
- No bibliography files were changed; process-local source keys remain listed
  pending catalog-wide bibliography consolidation.
- No PKA reference snapshots, PKA worklogs, catalog indexes, macros, or
  templates were modified.
