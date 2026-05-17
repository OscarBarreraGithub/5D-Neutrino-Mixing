# WA-WA-w9-vlq-4top Writer Worklog

**Agent:** WA-WA-w9-vlq-4top
**Batch:** Wave-9 PRIMARY `collider_rs` VLQ / four-top polish
**Cycle:** 1
**Timestamp:** 2026-05-17T17:08:10-04:00

## Scope

Touched only the assigned process entries and this worklog:

- `flavor_catalog/processes/collider_rs/CR008.{tex,yaml}`
- `flavor_catalog/processes/collider_rs/CR010.{tex,yaml}`
- `flavor_catalog/processes/collider_rs/CR014.{tex,yaml}`
- `flavor_catalog/worklogs/writer/wa_w9_vlq_4top.md`

## Required Reading Applied

- Read `flavor_catalog/AGENTIC_WORKFLOW.md` sections 3, 4, and 5.
- Read `flavor_catalog/signoff/by_process/L001.md` and `B001_B003.md`.
- Applied the CHK-1 carve-out uniformly for collider entries: observed mass exclusions, cross-section measurements, and signal/cross-section upper limits remain in `pdg_or_equivalent.values`; luminosity, collision energy, scan ranges, expected limits, projections, local RS-methodology numbers, and theory inputs stay in supporting or paper-era blocks.
- Used `flavor_catalog/processes/top_higgs_ew/T010.{tex,yaml}` for PRIMARY-entry style and sidecar transition shape.

## Per-Process Changes

### CR008: VLQ singlet T

- Tightened the Process and Current best limit(s) prose to distinguish the ATLAS singlet benchmark \(m_T>1.36\) TeV from the CMS all-branching-envelope \(m_T>1.48\) TeV context.
- Kept both measured mass exclusions in `pdg_or_equivalent.values`.
- Kept ATLAS/CMS luminosities in `supporting_measurements`, not `pdg_or_equivalent`.
- Clarified code-coverage evidence with explicit local `file:line` anchors: `scanParams/scan.py:399`, `scanParams/scan.py:523`, `scanParams/scan.py:535`, `quarkConstraints/scan.py:377`, and `quarkConstraints/scan.py:380`.
- Appended `WRITER-DONE` cycle-1 status transition.

### CR010: VLQ doublet (T,B)

- Added an explicit Collider RS partner-search framing sentence in the Process section.
- Preserved the PDG/ATLAS simultaneous doublet limits \(m_T,m_B>1.37\) TeV as the canonical assignment values.
- Kept simplified endpoint exclusions in `pdg_or_equivalent.values` because they are observed mass-exclusion limits, while noting in prose that they are not simultaneous doublet exclusions.
- Left luminosity and \(\sqrt{s}\) descriptors in `supporting_measurements`.
- Clarified code-coverage evidence with explicit `quarkConstraints/scan.py:359`, `quarkConstraints/scan.py:441`, `quarkConstraints/validation.py:108`, and `quarkConstraints/validation.py:145` anchors.
- Appended `WRITER-DONE` cycle-1 status transition.

### CR014: Four-top production

- Tightened the inclusive four-top prose to avoid treating the ATLAS/CMS cross-section measurements as a BSM excess.
- Kept inclusive cross-section measurements, the CMS observed \(850\) GeV benchmark exclusion, and the ATLAS \(21\)--\(119\) fb observed upper-limit range in `pdg_or_equivalent.values`.
- Kept CMS-B2G-25-005 luminosities, collision energies, scan ranges, and expected exclusion in `supporting_measurements`.
- Left the local quark-scan \(M_{\rm KK}^{\min}=47.26\) TeV comparison in `paper_era_reference`.
- Clarified code-coverage evidence with explicit `quarkConstraints/couplings.py:96`, `quarkConstraints/deltaf2.py:449`, `quarkConstraints/scan.py:377`, and `quarkConstraints/scan.py:441` anchors.
- Appended `WRITER-DONE` cycle-1 status transition.

## Source SHA Confirmation

Ran `sha256sum` over all tracked `CR008`, `CR010`, and `CR014` `.txt` source snapshots. The computed hashes match the sidecar `source_shas` entries for all cited snapshots:

- CR008: 11 source snapshots verified.
- CR010: 8 source snapshots verified.
- CR014: 19 source snapshots verified.

The `sha256sums.txt` files themselves are auxiliary checksum manifests and were not used as source keys in the sidecars.

## Code-Coverage Checks

Reran the collider/VLQ/four-top greps over:

`quarkConstraints qcd flavorConstraints neutrinos yukawa warpConfig solvers scanParams tests`

Results:

- VLQ/top-partner/bottom-partner/recast keyword grep returned no implementation matches.
- CR014/four-top/top-philic keyword grep returned no implementation matches.
- Collider-term grep returned only `tests/test_alpha_s.py:89`, an unrelated CMS/RunDec alpha-s test comment.
- Adjacent code remains low-energy only: `scanParams/scan.py:399`, `scanParams/scan.py:523`, `quarkConstraints/scan.py:359`, `quarkConstraints/scan.py:377`, `quarkConstraints/couplings.py:96`, and `quarkConstraints/deltaf2.py:449`.

## Bibliography Proposals

The requested path `flavor_catalog/references/catalog.bib` is not present in this checkout, and the batch instruction forbids touching files outside the assigned process paths plus this worklog. Proposed unresolved bibliography keys for a future bib merge:

- CR008: `PDG2025_TprimeVLQ`, `ATLAS2024_TsingletPair`, `CMS2023_VLQPairLeptonic`, `ATLAS2023_VLQPairMET`, `CMS2018_VLQPairLeptonic`, `CMS2018_TtoBW`, `CMS2017_VLQPairSingleLepton`, `ATLAS2015_VLQPair8TeV`, `AguilarSaavedra2009_TopPartners`, `AguilarSaavedraEtAl2013_VLQHandbook`, `LariEtAl2008_ColliderFlavourHighQ`.
- CR010: `PDG2025_TprimeListing`, `PDG2025_BprimeListing`, `ATLAS2018_TBCombination`, `ATLAS2024_VLQLeptonJets`, `CMS2023_VLQLeptonic`, `CMSReview2025_VLQ`, `CMS2017_TBBoosted`, `GarbersonGolling2013_ExoticQuarks`, `AguilarSaavedra2009_TopPartners`.
- CR014: `PDG2025_TopFourTop`, `ATLAS2023_FourTopObservation`, `CMS2023_FourTopObservation`, `CMSB2G25005_PublicResult`, `CMS2026_FourTopResonance`, `ATLASTopPhilic2024`, `CMS2014_FourTop8TeV`, `CMS2017_FourTop13TeV`, `CMS2018_FourTopMultilepton`, `ATLAS2019_FourTopLJDL`, `CMS2020_FourTopFullRun2Multilepton`, `ATLAS2021_FourTopCrossSection`, `CMS2023_FourTopEvidence`, `CaoChenLiu2016`, `DeBlasEberhardtKrause2018`, `BanelliEtAl2021_FourTopOperators`, `DarmeFuksMaltoni2021_TopPhilic`, `DarmeEtAl2025_BoostedTopPhilic`, `QuarkScanMethodology2026`.

## Open Issues

- No `catalog.bib` file exists at the requested path, so bibliography merge remains unresolved.
- CR014 uses the current CMS-B2G-25-005 public result / arXiv:2604.14058 as the direct resonance benchmark. Checker may decide to prefer only fully published ATLAS/CMS results for final signoff, but the observed CMS mass exclusion is explicitly sourced and stored as a measured collider limit.
- None of these entries should be implemented as a single hard-coded RS mass cut; all need a future collider recast or simplified-likelihood interface if promoted into live code.

## Status Transitions Appended

- CR008: `WRITER-INITIATED -> WRITER-DONE`, `agent_id: "WA"`, `cycle: 1`, `at: "2026-05-17T17:08:10-04:00"`.
- CR010: `WRITER-INITIATED -> WRITER-DONE`, `agent_id: "WA"`, `cycle: 1`, `at: "2026-05-17T17:08:10-04:00"`.
- CR014: `WRITER-INITIATED -> WRITER-DONE`, `agent_id: "WA"`, `cycle: 1`, `at: "2026-05-17T17:08:10-04:00"`.
