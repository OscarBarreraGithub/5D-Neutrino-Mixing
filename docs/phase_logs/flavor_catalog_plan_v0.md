# Flavor Catalog Planner v0

Date: 2026-05-16  
Branch for this planning commit: `paper/quark-scan-2026q2`  
Planning scope: design only. Do not create the catalog subdirectory in this round.

This plan assumes the catalog is a discovery-mode companion to the quark-scan paper, not a new constraint backend. The current repo already has an audited quark Delta F=2 lane and a lepton-sector mu->e gamma filter, but most of the PI's requested flavor-sensitive processes are not implemented and should be cataloged before any integration decision. Evidence for the current implemented surface:

- The methodology note says the scan uses five Delta F=2 observables, `epsilon_K`, `Delta m_K`, `Delta m_Bd`, `Delta m_Bs`, and `Delta m_D`, at [docs/quark_scan_methodology_note.tex](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/docs/quark_scan_methodology_note.tex:115).
- The repo v1 Delta F=2 default inputs include `epsilon_k`, `b_d`, `b_s`, and `d` systems at [quarkConstraints/deltaf2.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/deltaf2.py:209), with Delta m_K available as a helper at [quarkConstraints/deltaf2.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/deltaf2.py:755).
- The modern policy surface enumerates `epsilon_K`, `K`, `B_d`, `B_s`, and `D0`, but explicitly says it is policy-only and not a full backend at [quarkConstraints/modern/phenomenology.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/modern/phenomenology.py:23) and [quarkConstraints/modern/phenomenology.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/modern/phenomenology.py:165).
- The lepton-sector LFV implementation is mu->e gamma, documented and checked in [flavorConstraints/muToEGamma.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavorConstraints/muToEGamma.py:1), with the scan default using the MEG II 2025 limit at [scanParams/scan.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/scanParams/scan.py:32).
- The current paper scope explicitly excludes `neutrinos/`, `flavorConstraints/`, and `yukawa/` from the quark paper audit at [docs/quark_scan_methodology_note.tex](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/docs/quark_scan_methodology_note.tex:1173).

## A. Recommended Subdirectory Layout

Create the catalog later at repo root:

```text
flavor_catalog/
  README.md
  catalog_index.tex
  catalog_index.yaml
  catalog_master.tex
  latex/
    macros.tex
    process_template.tex
  processes/
    kaon/
      k_epsilon_k.tex
      k_epsilon_k.yaml
      ...
    charm/
    beauty/
    top_higgs_ew/
    charged_lepton/
    edm/
    neutrino_universality/
  references/
    catalog.bib
    by_process/
      k_epsilon_k/
        source_manifest.yaml
        sha256sums.txt
        pdg/
          pdg_snapshot_accessed_YYYYMMDD.txt
        hflav/
        flag/
        arxiv/
          0804.1954v2.pdf
          0804.1954v2.txt
        notes/
          extraction_notes.md
  figures/
    by_process/
      <process_id>/
        README.md
        source/
        generated/
  worklogs/
    pka/
    writer/
    checker/
    discovery/
    opus/
  signoff/
    round_000_index.md
    by_process/
      <process_id>.md
  tools/
    README.md
    snapshot_sources.sh
    validate_metadata.py
```

Rationale:

- `flavor_catalog/` should be intentionally separate from `quarkConstraints/`, `flavorConstraints/`, `neutrinos/`, and `docs/quark_scan_methodology_note.tex`. The present repo state marks lepton-sector tooling as follow-up scope, not audited quark-paper scope, at [docs/quark_scan_methodology_note.tex](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/docs/quark_scan_methodology_note.tex:1173).
- `processes/<family>/<process_id>.tex` keeps process drafts navigable and prevents one enormous file from hiding omissions.
- A sidecar `.yaml` next to every `.tex` is preferred over LaTeX-only comments because agents and validators can parse it without TeX.
- `catalog_index.yaml` is the machine-readable table of contents and status ledger. `catalog_index.tex` is the human-facing TOC. `catalog_master.tex` should compile all process `.tex` files after writer/checker approval.
- `references/by_process/<process_id>/` stores process-local raw materials and immutable snapshots. Shared citations live in `references/catalog.bib`.
- `worklogs/` records what each agent did, including failed source checks. `signoff/` records durable status decisions.

Branch recommendation:

- This planning commit stays on `paper/quark-scan-2026q2`, as requested.
- Actual catalog construction should happen on a separate branch, `flavor-catalog/2026q2`, cut from `paper/quark-scan-2026q2`. Rationale: the catalog will add many PDFs/text snapshots and long-lived draft files. Keeping it off the paper branch avoids churn while the rc paper artifacts remain frozen. Merge or cherry-pick only after Opus and PI sign-off.

Figures and external files:

- Store generated figures only under `flavor_catalog/figures/by_process/<process_id>/generated/`, with source code or notebooks under `source/`. No figures are required in v0 unless a process has a canonical plot that materially clarifies a constraint.
- Track at most 200 MB total under `flavor_catalog/references/`. Individual tracked source files should normally be under 5 MB. If the total would exceed 200 MB, prefer text snapshots, arXiv IDs, DOI metadata, and sha256 records in `source_manifest.yaml`; keep oversized PDFs external-only.
- Track arXiv PDFs and text extracts when useful. Do not track publisher PDFs unless the license explicitly allows redistribution. For PDG, HFLAV, FLAG, CKMfitter, UTfit, and experiment pages, track a minimal text snapshot of the relevant table/value, plus URL, access date, and sha256. Full PDG review PDFs should be external-only unless the license is explicitly cleared.

## B. Per-Process TeX Template

Recommendation: each process has both `processes/<family>/<process_id>.tex` and `processes/<family>/<process_id>.yaml`. The `.tex` should begin with a short comment pointing to the sidecar, but the sidecar is canonical metadata.

YAML sidecar skeleton:

```yaml
schema: flavor_catalog.process.v1
process_id: k_epsilon_k
family: kaon
standard_notation: "\\varepsilon_K"
process_name: "Indirect CP violation in neutral-kaon mixing"
owner_agent_id: "PKA-K-001"
writer_agent_id: null
checker_agent_id: null
opus_signoff_id: null
status: "DRAFT"   # DRAFT | WRITER-DONE | CHECKER-DONE | OPUS-APPROVED
version: "0.1.0"
last_updated: "2026-05-16"
source_shas:
  pdg_snapshot_accessed_YYYYMMDD.txt: "sha256:<fill>"
  1911.06822v2.pdf: "sha256:<fill>"
pdg_or_equivalent:
  source: "PDG Review of Particle Physics"
  access_date: "YYYY-MM-DD"
  value_summary: "<agent-filled exact value with uncertainty and year>"
code_coverage:
  status: "YES"
  evidence:
    - "quarkConstraints/deltaf2.py:209"
    - "quarkConstraints/deltaf2.py:729"
implementation_difficulty_if_missing: null
open_issues: []
```

LaTeX skeleton:

```tex
% flavor_catalog.process_id = k_epsilon_k
% Canonical metadata sidecar: k_epsilon_k.yaml
% Status is controlled by the YAML sidecar, not by this comment block.

\section*{Process}
\textbf{Standard notation:} \(\varepsilon_K\).  % fill exact process name

\subsection*{PDG value}
% Quote the canonical current PDG or equivalent value, uncertainty, CL if applicable,
% year/review edition, access date, and local snapshot path.

\subsection*{Relevance to RS / anarchic flavor}
% Explain tree-level KK gauge exchange, fermion localization, chirality enhancement,
% dipoles, Z/Higgs flavor, or neutrino-sector mechanism as applicable.

\subsection*{Post-2008 developments}
% What changed since the CFW/Perez-Randall era: experiment, lattice, SM theory,
% global fits, Belle II/LHCb/NA62/KOTO/MEG/ACME/Mu2e/COMET prospects.

\subsection*{Constraint validity / model dependence}
% Classify: robust tree-level, loop/dipole, SM-theory-limited, long-distance-limited,
% global-fit-dependent, lepton-extension-only, EDM-adjacent, etc.

\subsection*{Code coverage in this repo}
% YES/PARTIAL/NO. Include exact file:line evidence for YES/PARTIAL.

\subsection*{Implementation difficulty}
% LOW/MEDIUM/HIGH/BLOCKED with reason and blockers.

\subsection*{Key references}
% Use BibTeX keys from references/catalog.bib or process-local raw refs before merge.
```

Bibliography convention:

- During parallel PKA work, allow process-local `references/by_process/<process_id>/refs.bib` files to avoid merge conflicts.
- During writer/checker consolidation, merge accepted entries into `references/catalog.bib` with stable keys such as `PDG2026:Kepsilon`, `HFLAV2025:CharmMixing`, `NA62:KpPipNunu`, `BrodGorbahnStamou2020:EpsilonK`.
- `catalog_master.tex` should use only `references/catalog.bib`. Process-local `.bib` files are raw deposits, not final compile inputs.

## C. Initial Process List

Coverage aliases used below:

- `YES-D2-KCP`: implemented for kaon CP violation via `epsilon_k` default input at [quarkConstraints/deltaf2.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/deltaf2.py:209), constants at [quarkConstraints/deltaf2.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/deltaf2.py:613), and evaluator at [quarkConstraints/deltaf2.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/deltaf2.py:729).
- `PARTIAL-D2-KM`: Delta m_K helper exists at [quarkConstraints/deltaf2.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/deltaf2.py:755), and the methodology lists it as one of five scan observables at [docs/quark_scan_methodology_note.tex](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/docs/quark_scan_methodology_note.tex:115).
- `YES-D2-BD`: B_d mixing default input at [quarkConstraints/deltaf2.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/deltaf2.py:224), constants at [quarkConstraints/deltaf2.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/deltaf2.py:631), evaluator at [quarkConstraints/deltaf2.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/deltaf2.py:903).
- `YES-D2-BS`: B_s mixing default input at [quarkConstraints/deltaf2.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/deltaf2.py:238), constants at [quarkConstraints/deltaf2.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/deltaf2.py:645), evaluator at [quarkConstraints/deltaf2.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/deltaf2.py:922).
- `YES-D2-D0`: D0 mixing default input at [quarkConstraints/deltaf2.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/deltaf2.py:252), constants at [quarkConstraints/deltaf2.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/deltaf2.py:655), evaluator at [quarkConstraints/deltaf2.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/deltaf2.py:941).
- `YES-LFV-MUEG`: mu->e gamma implemented in [flavorConstraints/muToEGamma.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavorConstraints/muToEGamma.py:1) and [flavorConstraints/muToEGamma.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavorConstraints/muToEGamma.py:75), with scan default at [scanParams/scan.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/scanParams/scan.py:32).
- `NOT-SEEN`: no current implementation found by the planning grep across `quarkConstraints/`, `flavorConstraints/`, `scanParams/`, `neutrinos/`, `yukawa/`, `qcd/`, `scripts/`, `tests/`, and `docs/`, excluding notebooks. Later agents must rerun and refine this.

### Kaon sector

| ID | Process | One-line target / likely canonical source | Initial code coverage |
|---|---|---|---|
| K001 | `epsilon_K` | Indirect CP violation in K0-K0bar mixing; PDG kaon mixing plus SM theory inputs. | YES-D2-KCP |
| K002 | `Delta m_K` | K_L-K_S mass difference / neutral-kaon mixing magnitude; PDG. | PARTIAL-D2-KM |
| K003 | `epsilon'/epsilon` | Direct CP violation in K->pi pi; PDG/NA48/KTeV plus lattice theory. | NOT-SEEN |
| K004 | `K^+ -> pi^+ nu nubar` | Ultra-rare s->d nu nubar; NA62/PDG. | NOT-SEEN |
| K005 | `K_L -> pi^0 nu nubar` | CP-violating s->d nu nubar; KOTO/PDG. | NOT-SEEN |
| K006 | `K_L -> mu^+ mu^-` | Rare FCNC with long-distance component; PDG. | NOT-SEEN |
| K007 | `K_L -> e^+ e^-` | Rare helicity-suppressed neutral-kaon decay; PDG. | NOT-SEEN |
| K008 | `K_L -> pi^0 e^+ e^-` | Rare semileptonic FCNC with CP components; PDG. | NOT-SEEN |
| K009 | `K_L -> pi^0 mu^+ mu^-` | Muon counterpart of K_L -> pi0 ll; PDG. | NOT-SEEN |
| K010 | `K_S -> pi^0 e^+ e^-` | Short-distance/indirect input for K_L -> pi0 ee; PDG. | NOT-SEEN |
| K011 | `K_S -> pi^0 mu^+ mu^-` | Short-distance/indirect input for K_L -> pi0 mumu; PDG. | NOT-SEEN |
| K012 | `K_S -> mu^+ mu^-` | Rare short-lived neutral-kaon dimuon decay; LHCb/PDG. | NOT-SEEN |
| K013 | `K_L -> pi^0 gamma gamma` | Radiative rare kaon decay, likely related to PI's `K1 -> pi gamma` note; PDG. | NOT-SEEN |
| K014 | `K_S -> pi^0 gamma gamma` | Radiative short-lived kaon channel; clarify if this is the PI's intended `K1 -> pi gamma`. | NOT-SEEN |
| K015 | `K^+ -> pi^+ l^+ l^-` | Charged-kaon rare semileptonic modes, e and mu; PDG. | NOT-SEEN |
| K016 | `K_L -> pi^+ pi^- e^+ e^-` | CP/T-sensitive radiative four-body mode; PDG. | NOT-SEEN |
| K017 | `K^+ -> l^+ nu` / `R_K` | Lepton-universality ratio K->e nu / K->mu nu; PDG/FlaviaNet. | NOT-SEEN |
| K018 | `K_{l3}` | Semileptonic K->pi l nu and V_us / CKM unitarity; PDG/FLAG. | NOT-SEEN |
| K019 | `K_L -> e^\pm mu^\mp` | Charged-lepton flavor violation in kaon decay; PDG. | NOT-SEEN |
| K020 | `K^+ -> pi^+ e^\pm mu^\mp` | LFV charged-kaon semileptonic decay; PDG. | NOT-SEEN |
| K021 | `K_L -> pi^0 e^\pm mu^\mp` | LFV neutral-kaon semileptonic decay; PDG. | NOT-SEEN |
| K022 | `K -> pi pi` isospin amplitudes | Re A0/Re A2 and Delta I=1/2 inputs relevant to epsilon prime; PDG/lattice. | NOT-SEEN |

### Charm sector

| ID | Process | One-line target / likely canonical source | Initial code coverage |
|---|---|---|---|
| C001 | `D^0-D0bar` mixing, `x`, `y`, `Delta m_D` | Neutral charm mixing; HFLAV/PDG. | YES-D2-D0 |
| C002 | `|q/p|`, `phi_D` in D mixing | CP violation in neutral-charm mixing; HFLAV. | NOT-SEEN |
| C003 | `Delta A_CP(D^0 -> K^+K^-, pi^+pi^-)` | Direct charm CP violation; LHCb/HFLAV/PDG. | NOT-SEEN |
| C004 | `D^0 -> mu^+ mu^-` | Rare c->u ll FCNC; PDG/LHCb. | NOT-SEEN |
| C005 | `D^0 -> e^+ e^-` | Rare dilepton charm decay; PDG. | NOT-SEEN |
| C006 | `D^0 -> e^\pm mu^\mp` | LFV charm decay; PDG. | NOT-SEEN |
| C007 | `D^+ -> pi^+ mu^+ mu^-` | Rare c->u ll semileptonic mode; PDG/LHCb. | NOT-SEEN |
| C008 | `D^+ -> pi^+ e^\pm mu^\mp` | LFV semileptonic charm decay; PDG. | NOT-SEEN |
| C009 | `D^0 -> gamma gamma` | Radiative rare charm decay; PDG. | NOT-SEEN |
| C010 | `D -> rho gamma`, `D -> phi gamma` | Long-distance dominated radiative charm FCNC probes; PDG. | NOT-SEEN |
| C011 | `D_s -> l nu`, `D^+ -> l nu` | Leptonic charm decays, CKM/decay-constant consistency; HFLAV/FLAG/PDG. | NOT-SEEN |
| C012 | `D -> pi l nu`, `D -> K l nu` | Semileptonic charm decays, V_cd/V_cs and lepton universality; HFLAV/FLAG. | NOT-SEEN |

### Beauty sector

| ID | Process | One-line target / likely canonical source | Initial code coverage |
|---|---|---|---|
| B001 | `Delta m_d` / B_d mixing | Neutral B_d mixing magnitude; HFLAV/PDG/FLAG. | YES-D2-BD |
| B002 | `S_{psi K_S}` / `sin 2 beta` | CP phase in B_d -> J/psi K_S; HFLAV/PDG. | NOT-SEEN |
| B003 | `Delta m_s` / B_s mixing | Neutral B_s mixing magnitude; HFLAV/PDG/FLAG. | YES-D2-BS |
| B004 | `phi_s` in B_s -> J/psi phi | CP phase in B_s mixing; HFLAV/LHCb. | NOT-SEEN |
| B005 | `B_s -> mu^+ mu^-` | Clean rare b->s ll decay; HFLAV/PDG/LHCb/CMS/ATLAS. | NOT-SEEN |
| B006 | `B_d -> mu^+ mu^-` | Rare b->d ll decay; HFLAV/PDG. | NOT-SEEN |
| B007 | `B_s -> e^+ e^-`, `B_d -> e^+ e^-` | Electron rare B decays; PDG. | NOT-SEEN |
| B008 | `B_s -> tau^+ tau^-`, `B_d -> tau^+ tau^-` | Tau rare B decays; PDG/LHCb/Belle. | NOT-SEEN |
| B009 | `B^+ -> tau^+ nu_tau` | Charged-Higgs/charged-current sensitive seed; HFLAV/PDG. | NOT-SEEN |
| B010 | `B^+ -> mu^+ nu_mu` | Leptonic B decay; HFLAV/PDG. | NOT-SEEN |
| B011 | `B -> X_s gamma` / `b -> s gamma` | Inclusive radiative penguin; HFLAV/PDG. | NOT-SEEN |
| B012 | `B -> K^* gamma` | Exclusive b->s gamma; HFLAV/PDG. | NOT-SEEN |
| B013 | `B_s -> phi gamma` | Exclusive b->s gamma in B_s; HFLAV/PDG. | NOT-SEEN |
| B014 | `B -> rho gamma`, `B -> omega gamma` | b->d gamma exclusive modes; HFLAV/PDG. | NOT-SEEN |
| B015 | `B -> X_s l^+ l^-` | Inclusive b->s ll; HFLAV/PDG. | NOT-SEEN |
| B016 | `B -> K l^+ l^-` | Exclusive b->s ll branching fractions; HFLAV/LHCb/Belle. | NOT-SEEN |
| B017 | `B -> K^* l^+ l^-` | Angular observables and branching fractions, including P5 prime; HFLAV/LHCb. | NOT-SEEN |
| B018 | `R_K` | Lepton-universality ratio in B->K ll; HFLAV/LHCb/Belle. | NOT-SEEN |
| B019 | `R_{K^*}` | Lepton-universality ratio in B->K* ll; HFLAV/LHCb/Belle. | NOT-SEEN |
| B020 | `B_s -> phi l^+ l^-` | b->s ll in B_s; HFLAV/LHCb. | NOT-SEEN |
| B021 | `Lambda_b -> Lambda l^+ l^-` | Baryonic b->s ll channel; HFLAV/PDG/LHCb. | NOT-SEEN |
| B022 | `B -> K nu nubar` | b->s invisible rare decay; Belle II/PDG. | NOT-SEEN |
| B023 | `B -> K^* nu nubar` | Vector invisible rare decay; Belle/BaBar/PDG. | NOT-SEEN |
| B024 | `B -> pi nu nubar` | b->d invisible rare decay; PDG. | NOT-SEEN |
| B025 | `R_D` | B->D tau nu over light leptons; HFLAV. | NOT-SEEN |
| B026 | `R_{D^*}` | B->D* tau nu over light leptons; HFLAV. | NOT-SEEN |
| B027 | `R_{J/psi}` | B_c -> J/psi tau nu ratio; LHCb/PDG. | NOT-SEEN |
| B028 | `B_c -> tau nu` | B_c lifetime / charged-current tau constraint; PDG/theory. | NOT-SEEN |
| B029 | `B -> pi l nu` | V_ub exclusive semileptonic channel; HFLAV/FLAG. | NOT-SEEN |
| B030 | `B -> D l nu`, `B -> D^* l nu` | V_cb exclusive semileptonic channels; HFLAV/FLAG. | NOT-SEEN |
| B031 | `B -> X_u l nu`, `B -> X_c l nu` | Inclusive V_ub/V_cb channels; HFLAV. | NOT-SEEN |
| B032 | `Bbar -> pi Kbar` | Nonleptonic charmless B->pi K amplitudes and CP asymmetries; HFLAV/PDG. | NOT-SEEN |
| B033 | `B -> phi K_S` | Penguin CP asymmetry, seed related to `Bbar^0 -> phi`; HFLAV/PDG. | NOT-SEEN |
| B034 | `B_s -> phi phi` | b->s sbar s penguin and CP phase; LHCb/HFLAV. | NOT-SEEN |
| B035 | `B -> eta' K`, `B -> K K` | Additional charmless penguin modes; HFLAV/PDG. | NOT-SEEN |
| B036 | `B -> K mu tau`, `B -> K e tau`, `B -> K e mu` | LFV B decays; PDG/LHCb/Belle. | NOT-SEEN |
| B037 | `B_s -> mu tau`, `B_s -> e mu`, `B_d -> e mu` | LFV neutral-B dilepton modes; PDG/LHCb. | NOT-SEEN |

### Top, Higgs, and electroweak flavor

| ID | Process | One-line target / likely canonical source | Initial code coverage |
|---|---|---|---|
| T001 | `t -> c Z` | Top FCNC Z decay; ATLAS/CMS/PDG. | NOT-SEEN |
| T002 | `t -> u Z` | Top FCNC Z decay; ATLAS/CMS/PDG. | NOT-SEEN |
| T003 | `t -> c gamma` | Top FCNC photon decay; ATLAS/CMS/PDG. | NOT-SEEN |
| T004 | `t -> u gamma` | Top FCNC photon decay; ATLAS/CMS/PDG. | NOT-SEEN |
| T005 | `t -> c g` | Top FCNC gluon decay / production; ATLAS/CMS/PDG. | NOT-SEEN |
| T006 | `t -> u g` | Top FCNC gluon decay / production; ATLAS/CMS/PDG. | NOT-SEEN |
| T007 | `t -> c h` | Top-Higgs FCNC decay; ATLAS/CMS/PDG. | NOT-SEEN |
| T008 | `t -> u h` | Top-Higgs FCNC decay; ATLAS/CMS/PDG. | NOT-SEEN |
| T009 | single-top FCNC `qg -> t` | Production interpretation of tqg FCNC; ATLAS/CMS. | NOT-SEEN |
| T010 | `Z -> b bbar`, `R_b` | Z-pole bottom coupling constraint; PDG/LEP/SLC. | NOT-SEEN |
| T011 | `A_FB^b`, `A_b` | Z-pole bottom asymmetries; PDG/LEP/SLC. | NOT-SEEN |
| T012 | `Z -> c cbar`, `R_c` | Z-pole charm coupling constraint; PDG/LEP/SLC. | NOT-SEEN |
| T013 | `A_FB^c`, `A_c` | Z-pole charm asymmetries; PDG/LEP/SLC. | NOT-SEEN |
| T014 | `Z -> b s`, `Z -> b d`, `Z -> s d` | Flavor-changing Z decays; PDG/LEP/LHC searches. | NOT-SEEN |
| T015 | `Z -> e mu` | LFV Z decay; ATLAS/CMS/PDG. | NOT-SEEN |
| T016 | `Z -> e tau` | LFV Z decay; ATLAS/CMS/PDG. | NOT-SEEN |
| T017 | `Z -> mu tau` | LFV Z decay; ATLAS/CMS/PDG. | NOT-SEEN |
| T018 | `h -> mu tau` | Higgs LFV decay; ATLAS/CMS. | NOT-SEEN |
| T019 | `h -> e tau` | Higgs LFV decay; ATLAS/CMS. | NOT-SEEN |
| T020 | `h -> e mu` | Higgs LFV decay; ATLAS/CMS. | NOT-SEEN |
| T021 | `h -> b s`, `h -> b d`, `h -> s d` | Higgs quark-flavor violation; LHC/PDG limits. | NOT-SEEN |
| T022 | W/Z charged-current universality | Electroweak lepton/quark flavor universality tests; PDG. | NOT-SEEN |

### Charged-lepton flavor violation and lepton universality

| ID | Process | One-line target / likely canonical source | Initial code coverage |
|---|---|---|---|
| L001 | `mu -> e gamma` | Dipole LFV; MEG/MEG II/PDG. | YES-LFV-MUEG |
| L002 | `mu -> 3e` | Four-lepton/contact LFV; SINDRUM/Mu3e/PDG. | NOT-SEEN |
| L003 | `mu-e` conversion in Al | Coherent conversion; SINDRUM II/Mu2e/COMET projections. | NOT-SEEN |
| L004 | `mu-e` conversion in Au | Existing SINDRUM II strong limit; PDG. | NOT-SEEN |
| L005 | `mu-e` conversion in Ti | Existing conversion limits; PDG. | NOT-SEEN |
| L006 | muonium-antimuonium conversion | Delta L_mu - Delta L_e = 2 LFV; PDG. | NOT-SEEN |
| L007 | `tau -> mu gamma` | Tau dipole LFV; Belle/BaBar/Belle II/PDG. | NOT-SEEN |
| L008 | `tau -> e gamma` | Tau dipole LFV; Belle/BaBar/Belle II/PDG. | NOT-SEEN |
| L009 | `tau -> 3 mu` | Four-lepton tau LFV; Belle/BaBar/LHCb/PDG. | NOT-SEEN |
| L010 | `tau -> 3 e` | Four-lepton tau LFV; Belle/BaBar/PDG. | NOT-SEEN |
| L011 | `tau -> mu e e` | Mixed-flavor tau LFV; PDG. | NOT-SEEN |
| L012 | `tau -> e mu mu` | Mixed-flavor tau LFV; PDG. | NOT-SEEN |
| L013 | `tau -> mu eta` | Hadronic tau LFV; Belle/BaBar/PDG. | NOT-SEEN |
| L014 | `tau -> e eta` | Hadronic tau LFV; Belle/BaBar/PDG. | NOT-SEEN |
| L015 | `tau -> mu eta'` | Hadronic tau LFV; Belle/BaBar/PDG. | NOT-SEEN |
| L016 | `tau -> e eta'` | Hadronic tau LFV; Belle/BaBar/PDG. | NOT-SEEN |
| L017 | `tau -> mu pi^0`, `tau -> e pi^0` | Hadronic tau LFV; PDG. | NOT-SEEN |
| L018 | `tau -> mu rho`, `tau -> e rho` | Vector-meson tau LFV; PDG. | NOT-SEEN |
| L019 | `tau -> mu phi`, `tau -> e phi` | Strange vector-meson tau LFV; PDG. | NOT-SEEN |
| L020 | `tau -> l K_S` | Tau LFV with neutral kaon; PDG. | NOT-SEEN |
| L021 | pion lepton-universality `R_pi` | pi->e nu / pi->mu nu; PDG. | NOT-SEEN |
| L022 | tau lepton-universality ratios | tau->l nu nu and tau->h nu universality; HFLAV/PDG. | NOT-SEEN |
| L023 | neutrino trident production | Muon-neutrino trident constraint on lepton gauge/contact operators; CCFR/CHARM-II. | NOT-SEEN |
| L024 | PMNS oscillation observables | Lepton mixing angles, mass splittings, delta_CP; NuFIT/PDG. | NOT-SEEN |
| L025 | neutrinoless double beta decay | Lepton-number violation adjacent to neutrino extension; KamLAND-Zen/GERDA/LEGEND/PDG. | NOT-SEEN |

### EDM-adjacent flavor and CP

| ID | Process | One-line target / likely canonical source | Initial code coverage |
|---|---|---|---|
| E001 | electron EDM | CP/flavor alignment via lepton dipoles; ACME/PDG. | NOT-SEEN |
| E002 | muon EDM | Lepton dipole CP phase; PDG. | NOT-SEEN |
| E003 | tau EDM | Tau dipole CP phase; Belle/PDG. | NOT-SEEN |
| E004 | neutron EDM | Quark EDM/chromo-EDM/Weinberg constraints; nEDM/PDG. | NOT-SEEN |
| E005 | proton/deuteron EDM prospects | Hadronic EDM future sensitivity; storage-ring proposals. | NOT-SEEN |
| E006 | mercury EDM | Diamagnetic atom CP constraint; PDG. | NOT-SEEN |
| E007 | radium/xenon EDM | Complementary hadronic/nuclear EDM probes; PDG. | NOT-SEEN |
| E008 | quark chromo-EDM bounds | EFT-level translation relevant to RS loops; lattice/QCD sum rules. | NOT-SEEN |
| E009 | Weinberg three-gluon operator | CP-odd gluonic operator relevant to KK loops; EDM global fits. | NOT-SEEN |
| E010 | CPV in D/B radiative decays | EDM-adjacent dipole phases in flavor transitions; HFLAV/PDG. | NOT-SEEN |

This initial list has 128 entries. It is intentionally overcomplete. Discovery agents should add missing baryon, hyperon, and nuclear processes if they find credible RS-flavor relevance.

## D. Agent Role Specification

### Process-Knowledge Agent (PKA)

Recommended model: Codex `gpt-5.5` with `xhigh` reasoning. Run in background. Assign one process by default; allow up to three tightly coupled processes only when the same literature and PDG table covers all of them, e.g. `t->uZ` and `t->cZ`, or `tau->mu eta` and `tau->e eta`.

Expected wall time: 60-150 minutes per process. Parallelism budget: 8 concurrent PKAs initially, rising to 12 only if repository writes are process-disjoint and source downloads are rate-limited politely.

Exact spawn prompt template:

```text
You are PKA-{agent_id} for the 5D-Neutrino-Mixing flavor catalog.
Work only under flavor_catalog/ once that directory exists. Do not edit quarkConstraints/,
flavorConstraints/, neutrinos/, yukawa/, qcd/, scan outputs, or paper docs.

Assigned process(es):
- {process_id}: {standard_notation} -- {plain_name}

For each process:
1. Identify the canonical current experimental value or limit from PDG, HFLAV, FLAG,
   experiment papers, or another accepted source. Record value, uncertainty or CL,
   source edition/year, URL, and access date.
2. Save a tracked minimal source snapshot under references/by_process/{process_id}/:
   PDG/HFLAV/FLAG text extract where allowed; arXiv PDF and pdftotext extract for
   arXiv preprints; source_manifest.yaml and sha256sums.txt for all files.
3. Find post-2008 experimental/theory developments relative to CFW/Perez-Randall-era
   literature. Include arXiv IDs and explain why each source matters.
4. Run the catalog coverage grep from the plan and record whether the process is
   implemented, partially implemented, or absent. Cite exact file:line evidence.
5. Draft processes/<family>/{process_id}.tex using the process template and write
   processes/<family>/{process_id}.yaml with source SHAs.
6. Write worklogs/pka/{process_id}.md with search queries, sources rejected, and
   unresolved issues.

Success criteria:
- The PDG/equivalent value can be traced to a local snapshot or manifest entry.
- The .tex has all required sections filled, even if some say "unknown; checker must verify".
- Code coverage has exact file:line evidence for YES/PARTIAL and explicit grep commands for NO.
- No publisher PDFs are tracked unless redistribution is explicitly licensed.
```

Deliverables:

- `processes/<family>/<process_id>.tex`
- `processes/<family>/<process_id>.yaml`
- `references/by_process/<process_id>/source_manifest.yaml`
- `references/by_process/<process_id>/sha256sums.txt`
- source snapshots/PDFs/text extracts subject to the license policy
- `worklogs/pka/<process_id>.md`

### Writer Agent (WA)

Recommended model: Codex `gpt-5.5` `xhigh`. Run in background. Assign batches of 3-5 processes from the same family. The WA should not be the same agent as the PKA for paper-quality batches. A PKA may draft raw text, but WA must independently structure, tighten, and normalize it.

Expected wall time: 2-4 hours per 3-5 process batch. Parallelism budget: 4 concurrent WAs.

WA task:

- Convert PKA raw deposits into coherent paper-quality process entries.
- Normalize constraint classifications and implementation difficulty.
- Merge raw process-local refs into a proposed shared `catalog.bib` patch or a batch-local merge file.
- Leave disputed or weak claims marked with `\textbf{CHECK}` rather than smoothing over uncertainty.

### Checker Agent (CA)

Recommended model: Codex `gpt-5.5` `xhigh`. Run in background. A CA must never check a batch it wrote. Expected wall time: 2-5 hours per batch. Parallelism budget: 4 concurrent CAs.

Checker checklist:

- Reopen every PDG/HFLAV/FLAG/experiment source and confirm the value, uncertainty, CL, year/edition, and access date.
- Verify every arXiv citation exists, matches the claim, and is not confused with a published-only version.
- Confirm local snapshot files exist and match `source_shas`.
- Rerun code coverage grep and inspect exact file:line citations.
- Confirm RS relevance is model-specific and not overstated. Mark whether it applies to quark RS scan, lepton extension, both, or EDM-adjacent loops only.
- Confirm implementation difficulty rubric is applied consistently.
- Check that no process silently depends on `quarkConstraints/` integration.
- Record PASS/FAIL with actionable fixes in `worklogs/checker/<batch_id>.md`.

### Discovery Agent (DA)

Recommended model: Codex `gpt-5.5` `xhigh`. Run after the first WA+CA round, and again after fixes. Expected wall time: 3-6 hours per discovery round. Parallelism budget: 3 concurrent DAs with non-overlapping scopes:

- DA-1: PDG rare/forbidden decays and kaon/charm coverage.
- DA-2: HFLAV/Belle II/LHCb B and tau coverage.
- DA-3: RS/extra-dimensional flavor literature, EDMs, top/Higgs/Z flavor.

DA success criteria:

- Produce a proposed additions table with process id, source, one-line rationale, and whether it is already covered by an existing entry.
- Identify duplicate or over-split entries that should be merged.
- Explain search coverage, not just final additions.

Run at least two DA rounds. Converge only after two consecutive discovery rounds find at most one genuinely new process each and every existing process has CA pass plus Opus approval.

### Opus Sign-Off Agent

Recommended model: Opus, foreground for final batch decisions. Expected wall time: 30-60 minutes per 10-15 process batch. Parallelism budget: 1-2 sign-offs at a time because sign-off should preserve global consistency.

Opus responsibilities:

- Arbitrate unresolved WA/CA disagreements.
- Check tone, overclaiming, and paper-readiness.
- Approve status transition to `OPUS-APPROVED`.
- Escalate physics judgments to PI when model dependence cannot be resolved from sources.

## E. Orchestration Sequence

1. PI green-light is already given for planning. Before execution, Claude should confirm that creating `flavor_catalog/` on `flavor-catalog/2026q2` is approved.
2. Claude creates `flavor-catalog/2026q2` from `paper/quark-scan-2026q2` and scaffolds only the layout, template, and index files.
3. Claude spawns PKAs. With 128 initial entries, use 44 PKA tasks with M=3 maximum, but force M=1 for high-risk seed processes: `epsilon'/epsilon`, `K->pi nu nubar`, `B_s->mumu`, `b->s gamma`, `R_K(*)`, `R_D(*)`, `mu-e conversion`, electron/neutron EDM, and all PI-ambiguous kaon radiative modes.
4. PKAs deposit drafts and raw materials. Claude validates only file existence, metadata schema, sha256sums, and no writes outside `flavor_catalog/`.
5. Claude groups deposits into WA batches of 3-5 processes by family. WAs produce writer-quality drafts.
6. Claude assigns different CAs to those batches. A CA cannot check its own WA batch.
7. Mismatches go back to the WA with the CA worklog attached. Iterate WA->CA until CA passes.
8. DA round 1 starts after at least 50% of the initial process list has CA pass. Add missing processes and send them through PKA->WA->CA.
9. DA round 2 starts after all initial and DA-1 additions have CA pass. Repeat if more than one genuinely new process appears.
10. Opus signs off each approved batch. Disagreements unresolved by Opus go to the PI.
11. Build `catalog_master.tex`. The master compile is a consistency check, not a physics approval substitute.
12. Commit and push the catalog branch. Merge to paper branch only after PI approval.

Convergence rule:

- Two discovery rounds in a row find at most one new genuine process each.
- Every process has `CHECKER-DONE` and `OPUS-APPROVED`.
- `catalog_master.tex` compiles.
- `catalog_index.yaml` has no `DRAFT`, `WRITER-DONE`, or stale source snapshot statuses.

## F. Coding-Coverage Annotation Plan

The current planning grep indicates the implemented process surface is limited to Delta F=2 neutral meson mixing plus mu->e gamma. The catalog agents must rerun coverage, because new commits may change this.

Required grep commands:

```bash
rg -n "epsilon_K|epsilon_k|Delta m_K|DELTA_M_K|evaluate_delta_mk|B_d|b_d|B_s|b_s|D0|d_mix|mu.?e.?gamma|muToEGamma|BR_LIMIT|MEG" quarkConstraints flavorConstraints scanParams neutrinos yukawa qcd scripts tests docs -g '!*.ipynb'
rg -n "b.?->.?s|X_s|K\\*|R_K|R_D|B_s.*mu|pi.*nu|nu.*nubar|epsilon.?prime|eps.?prime|t.?->|Z.?->|h.?->|EDM|dipole|conversion|tau.?->|mu.?->.?3e" quarkConstraints flavorConstraints scanParams neutrinos yukawa qcd scripts tests docs -g '!*.ipynb'
rg -n "MODERN_PHENOMENOLOGY_SYSTEM_IDS|DEFAULT_DELTA_F2_INPUTS|evaluate_bd_mixing|evaluate_bs_mixing|evaluate_d0_mixing|check_mu_to_e_gamma" quarkConstraints flavorConstraints scanParams tests docs -g '!*.ipynb'
```

Modules to inspect manually for YES/PARTIAL claims:

- `quarkConstraints/deltaf2.py`
- `quarkConstraints/modern/phenomenology.py`
- `quarkConstraints/modern/matching.py`
- `quarkConstraints/paper_0710_1869/eft_deltaf2/`
- `flavorConstraints/muToEGamma.py`
- `scanParams/scan.py`
- `neutrinos/`, `yukawa/`, and `qcd/` for supporting inputs only

Implementation difficulty rubric:

- LOW: same Delta F=2 or mu->e gamma pattern already exists; needs only a new input value, sidecar metadata, or a direct ratio wrapper.
- MEDIUM: needs a new observable formula or new Wilson/operator mapping, but no new RG machinery and no global fit.
- HIGH: needs new RG running, lattice/hadronic matrix elements not currently represented, loop matching, angular observable likelihoods, or nontrivial SM subtraction.
- BLOCKED: canonical theory prediction or experimental likelihood is unavailable, licensing prevents source snapshotting, or PI must decide model scope.

Initial difficulty expectations:

- Existing Delta F=2 entries: LOW for cataloging, MEDIUM for production-grade integration because endpoint/systematic issues remain documented in the methodology note at [docs/quark_scan_methodology_note.tex](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/docs/quark_scan_methodology_note.tex:1055).
- `mu->e gamma`: LOW for cataloging, MEDIUM for integration update because the lepton-sector scan uses the MEG II 2025 default at [scanParams/scan.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/scanParams/scan.py:32), while the implementation docstring still records the Perez-Randall paper-era coefficient at [flavorConstraints/muToEGamma.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavorConstraints/muToEGamma.py:20).
- Rare K/B/charm decays: MEDIUM to HIGH, depending on whether a clean short-distance Wilson basis exists.
- `b->s ll` anomalies and angular observables: HIGH because they require likelihoods, correlated bins, and SM form-factor treatment.
- EDMs and chromo-EDMs: HIGH because they require loop matching and hadronic/nuclear translation.
- Top/Higgs/Z FCNC: MEDIUM for limit cataloging, HIGH for RS matching and integration.

## G. Iteration and Convergence

Expected rounds:

- PKA raw deposit: one round, with source-fix micro-iterations only.
- WA/CA: two full rounds expected for most batches. Rare kaon, B anomalies, and EDM batches may need three.
- DA: two rounds minimum after the first CA pass wave.
- Opus: one sign-off pass per final batch, with one possible corrective loop.

Disagreement resolution:

- WA and CA disagreements first go to Claude for routing only, not physics judgment.
- Opus arbitrates source interpretation, wording, and status if the disagreement is resolvable from literature.
- PI arbitrates scope/model-dependence questions, including whether a process should be included as a hard constraint, a future projection, or contextual only.

Versioning:

- Process sidecar version starts at `0.1.0`.
- PKA deposit increments patch: `0.1.1`, `0.1.2`.
- WA major draft increments minor: `0.2.0`.
- CA-approved draft increments: `0.3.0`.
- Opus-approved draft becomes `1.0.0`.
- Any post-approval PDG/source update increments minor: `1.1.0`; typo-only changes increment patch.

Sign-off logs:

- One process log at `signoff/by_process/<process_id>.md`.
- One round index at `signoff/round_<NNN>_index.md`.
- Each log entry records process id, version, status transition, agent id, date, source manifest sha, and unresolved caveats.

## H. Budget and Risk

Agent-hour estimate for the 128-entry initial catalog:

- PKA round: 190-320 agent-hours, assuming 1.5-2.5 hours per process, with grouped trivial pairs.
- WA round: 90-140 agent-hours, assuming 32-40 batches at 2-4 hours each.
- CA round: 120-180 agent-hours, assuming 32-40 batches at 3-5 hours each.
- DA rounds: 24-54 agent-hours total for two rounds across three DAs.
- Opus approval: 12-24 agent-hours.
- Orchestrator overhead: 30-50 human/Claude hours for routing, schema validation, conflict handling, and final compile.

Risk register:

| Risk | Impact | Mitigation |
|---|---|---|
| Stale PDG value passes because PKA and CA read the same old source | Wrong catalog values | Require current PDG/equivalent access date, local snapshot, and CA independent reopening of source. |
| PDG/HFLAV page changes after snapshot | Non-reproducible citation | Store text snapshot, URL, access date, and sha256 in `source_manifest.yaml`. |
| Publisher PDFs get committed | License problem | Track arXiv preprints and minimal text extracts only; publisher PDFs external-only unless license-cleared. |
| One agent shallowly covers too many processes | Missing constraints | Cap PKA at one process by default, three only for tightly coupled channels. |
| Checker focuses on prose rather than values | False approval | CA checklist starts with value/source verification before writing quality. |
| Rare-decay theory overclaimed as robust RS bound | Misleading physics | Mandatory validity/model-dependence classification in every process file. |
| Duplicate process IDs or bibliography keys | Compile and maintenance failures | `catalog_index.yaml` and `validate_metadata.py` enforce unique IDs and keys. |
| Large tracked snapshots bloat repo | Unusable branch | 200 MB tracked cap, 5 MB per file guidance, external-only oversized sources. |
| Existing quark scan accidentally integrated with draft catalog | Paper branch instability | Work on `flavor-catalog/2026q2`, keep all files under `flavor_catalog/`, no code changes. |

Acceptance criteria from the PI's view:

- Every seed process from the PI appears in `catalog_index.yaml` or is explicitly mapped to a clarified/renamed process.
- At least two discovery rounds have converged.
- Every process has a source snapshot or manifest, exact value/limit, RS relevance, post-2008 update, validity classification, code coverage note, and implementation difficulty.
- All `YES/PARTIAL` code coverage claims cite file:line evidence.
- `catalog_master.tex` compiles and the bibliography resolves.
- Opus-approved status exists for every process.
- The PI can decide, process by process, which entries should become future code constraints.

## I. Open Questions for the PI

1. Clarify the seed `K1 -> pi gamma`. Likely interpretations are `K_S -> pi^0 gamma gamma`, `K_L -> pi^0 gamma gamma`, or another radiative kaon shorthand. A single-photon `K -> pi gamma` two-body mode is not the usual neutral-kaon rare-decay constraint.
2. Confirm scope. This plan includes quark, charged-lepton, neutrino, electroweak, Higgs/top FCNC, and EDM-adjacent processes by default. The PI can prune later, but discovery should not start quark-only.
3. Confirm branch policy. Recommendation is execution on `flavor-catalog/2026q2`, not directly on `paper/quark-scan-2026q2`.
4. Confirm whether full PDG PDFs may be cached, or whether minimal text snapshots plus source manifests are preferred for licensing.
5. Confirm whether the catalog is a companion/follow-up artifact rather than part of rc1.1/rc2. This plan treats it as a companion unless the PI says otherwise.
6. Confirm whether future implementation should target `quarkConstraints/modern/` only, or whether the older `quarkConstraints/deltaf2.py` v1 path remains acceptable for new constraints.

===FLAVOR_CATALOG_PLAN_V0_END===
