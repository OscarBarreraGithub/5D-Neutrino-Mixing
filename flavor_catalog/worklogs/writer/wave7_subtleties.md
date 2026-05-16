# Wave 7 Subtlety Writer Worklog

Agent: subtlety-writer-wave7
Branch: flavor-catalog/2026q2
Timestamp used for sidecar updates: 2026-05-16T18:26:12-04:00

This pass threaded only the six PI-approved cross-cutting subtleties from the
May 15 and May 16 external Deep Research review summaries into existing process
entries.  The edits were limited to the existing constraint-validity/model-
dependence subsections, with one short paragraph or two compact sentences per
entry.  No PDG/equivalent values, uncertainties, citations, source snapshots,
or subsection ordering were changed.  Each touched process sidecar received a
single appended `SUBTLETY-ADDED` status-history item with agent
`subtlety-writer-wave7`; no other sidecar fields were intentionally modified.

## Files Touched

- `flavor_catalog/processes/beauty/B011.tex` and
  `flavor_catalog/processes/beauty/B011.yaml`: added subtlety D, photon
  polarization in radiative B.  The TeX now distinguishes the inclusive
  \(B\to X_s\gamma\) branching fraction as a \(|C_7,C_8|\) dipole constraint
  from the RS right-handed-dipole signature, which needs exclusive polarization
  and time-dependent CP observables sensitive to \(C_7'\), with B012 noted as
  the future exclusive entry.

- `flavor_catalog/processes/beauty/B015.tex` and
  `flavor_catalog/processes/beauty/B015.yaml`: added subtlety A, custodial-
  protection dependence for rare-b constraints.  The entry now states that
  reduced \(Z b_L b_L\) pressure in custodial RS can make inclusive
  \(b\to s\ell\ell\) relatively more discriminating, and that RS bounds require
  specifying custodial protection, fermion embeddings, and brane kinetic terms.

- `flavor_catalog/processes/beauty/B017.tex` and
  `flavor_catalog/processes/beauty/B017.yaml`: added subtleties A and C.  The
  b -> s l+l- umbrella row now carries the custodial-RS caveat for semileptonic
  channels and explicitly frames \(R_K/R_{K^*}\) as precision LFU null tests,
  not as inputs with a pre-2023 anomaly prior.

- `flavor_catalog/processes/beauty/B018.tex` and
  `flavor_catalog/processes/beauty/B018.yaml`: added subtleties A and C.  The
  \(R_K\) row now notes custodial-dependence of its discriminating power in a
  combined rare-b picture and repeats the modern LFU-null-test framing.

- `flavor_catalog/processes/beauty/B019.tex` and
  `flavor_catalog/processes/beauty/B019.yaml`: added subtleties A and C.  The
  \(R_{K^*}\) row now carries the same custodial-RS model-dependence warning
  and the instruction not to use the pre-2023 anomaly narrative as a prior.

- `flavor_catalog/processes/charm/C001.tex` and
  `flavor_catalog/processes/charm/C001.yaml`: added subtlety B.  Neutral charm
  mixing is now explicitly promoted to a leading up-sector diagnostic in
  down-aligned or kaon-protected RS variants.

- `flavor_catalog/processes/charm/C002.tex` and
  `flavor_catalog/processes/charm/C002.yaml`: added subtlety B.  CP violation
  in charm mixing now carries the same down-alignment/kaon-protection triage
  language.

- `flavor_catalog/processes/charm/C003.tex` and
  `flavor_catalog/processes/charm/C003.yaml`: added subtlety B.  Direct charm
  CP violation is now marked as a leading up-sector observable when down-sector
  alignment suppresses kaon and B pressure.

- `flavor_catalog/processes/charm/C004.tex` and
  `flavor_catalog/processes/charm/C004.yaml`: added subtlety B.  Rare
  \(D^0\to\mu^+\mu^-\) is now covered by the same down-aligned RS caveat.

- `flavor_catalog/processes/top_higgs_ew/T001.tex`,
  `T001.yaml`, `T002.tex`, `T002.yaml`, `T005.tex`, `T005.yaml`, `T006.tex`,
  and `T006.yaml`: added subtlety B to the top FCNC rows \(t\to cZ\),
  \(t\to uZ\), \(t\to cg\), and \(t\to ug\).  Each now says top FCNCs become
  leading, not secondary, diagnostics in down-aligned or kaon-protected RS.

- `flavor_catalog/processes/top_higgs_ew/T007.tex` and
  `flavor_catalog/processes/top_higgs_ew/T007.yaml`: added subtleties B and F.
  The row now promotes top-Higgs FCNCs in down-aligned variants and adds the
  kinematics nit that on-shell \(h\to tc\) is forbidden for a 125 GeV Higgs;
  \(Htc\) constraints enter through \(t\to Hc\) and associated single-top + H
  production.

- `flavor_catalog/processes/top_higgs_ew/T010.tex` and
  `flavor_catalog/processes/top_higgs_ew/T010.yaml`: added subtlety A.  The
  Zbb row now explicitly warns that custodial RS can make the observable less
  direct as a mass setter or more diagnostic of embeddings, and that RS bounds
  require the custodial/embedding/brane-kinetic-term assumptions.

- `flavor_catalog/processes/top_higgs_ew/EW001.tex`,
  `EW001.yaml`, `EW002.tex`, and `EW002.yaml`: added subtlety A.  The oblique
  and CKM-unitarity EW rows now tie their RS mass-limit interpretation to the
  same custodial-protection, fermion-embedding, and brane-kinetic-term
  assumptions.

- `flavor_catalog/processes/edm_neutrino/E001.tex`, `E001.yaml`, `E004.tex`,
  `E004.yaml`, `E006.tex`, `E006.yaml`, `E008.tex`, `E008.yaml`, `E009.tex`,
  and `E009.yaml`: added subtlety E.  The electron, neutron, mercury, qCEDM,
  and Weinberg-operator EDM rows now say anarchic CP phases should be reviewed
  jointly with flavor rows sharing Wilson coefficients, with cross-references
  K001, K003, B011, B033, and B034, rather than isolated in an appendix.

- `flavor_catalog/worklogs/writer/wave7_subtleties.md`: this worklog.
