# External Review: deepresearch_may15 Compared With Flavor Catalog v0.1

## Section 1: Processes / observables in the external research that we are MISSING

The strongest genuine miss is the radiative-B chirality block beyond the inclusive
`B -> X_s gamma` rate. Deep Research explicitly says the modern `b -> s gamma`
block should include `B -> X_s gamma`, `B -> K* gamma`, `B_s -> phi gamma`,
CP observables, polarization observables, `b -> d gamma`, `C7`, `C7'`, and
time-dependent radiative CP asymmetries (`deepresearch_may15.txt:691-722`);
it later singles out photon polarization because RS can generate sizeable
right-handed dipoles (`deepresearch_may15.txt:897-901`). The catalog has B011
for inclusive `B -> X_s gamma`, and B011 itself says it is separated from
exclusive `B -> K* gamma` and `B_s -> phi gamma` (`B011.tex:13-19`). The plan
listed B012/B013/B014 for exclusive `b -> s gamma` and `b -> d gamma`
(`flavor_catalog_plan_v1.md:273-276`), but DA-4 deferred "exclusive radiative
variants after inclusive b->s gamma" (`round_004_convergence.md:35`). My
judgment: the deferral is only partly sound. It is defensible for a first
rate-only catalog, but not for an RS chirality catalog: `C7'` and photon
polarization are exactly the observables that distinguish right-handed dipoles
from a simple inclusive-rate constraint.

Top FCNC photon modes are also missing. Deep Research ranks top FCNCs as
`t -> cZ`, `t -> cH`, `t -> c gamma` plus production analogues
(`deepresearch_may15.txt:426-467`) and says the seed list should add `t -> qH`,
`t -> q gamma`, and `tqg` production (`deepresearch_may15.txt:832-838`). The
catalog covers `t -> cZ`/`t -> uZ` (T001/T002), `t -> cg`/`t -> ug`
(T005/T006), and `t -> Hc` (T007), with values fact-checked in the inventory
(`factcheck_status.md:74-79`). There is no T003/T004 process file for
`t -> c gamma` or `t -> u gamma`, even though plan v1 listed them
(`flavor_catalog_plan_v1.md:307-308`). Judgment: this is a real catalog miss,
not just a harmless tail. The DA-4 deferred list does not explicitly close
T003/T004 (`round_004_convergence.md:36-38`), so the current 75-row closure is
ambiguous here.

Deep Research includes direct KK-gluon, vector-like-quark, and custodial-partner
searches as a ranked item (`deepresearch_may15.txt:500-555`) and later adds
high-pT dilepton/dijet tails as the cleaner handle once resonances are broad or
off-shell (`deepresearch_may15.txt:914-935`). The catalog does not have a direct
collider-search process family; direct KK-gluon scale statements appear only as
CFW theory context in process files such as K001 (`K001.tex:40-47`) and not as
current LHC reinterpretations. Judgment: I would not force these into the 75
process entries, because Deep Research itself warns that old resonance limits
are benchmark-, width-, and coupling-dependent (`deepresearch_may15.txt:1006-1009`).
But the catalog should probably add a synthesis note or deferred-scope record,
because RS readers will expect collider tails and VLQ/custodial-partner limits
to be named explicitly.

The electroweak/global-fit block is broader in Deep Research than in the
catalog. It asks for CKM/W universality, `mW`, `S,T,U,W,Y`, and charged-current
universality in one correlated fit (`deepresearch_may15.txt:575-616`), and in
the Z-pole seed audit it says to include `Z -> cc`, `Gamma_Z`, `S,T,U,W,Y`, and
`mW` in custodial analyses (`deepresearch_may15.txt:844-855`). The catalog has
EW001 for `S,T,U` (`EW001.tex:13-18`), EW002 for first-row CKM unitarity
(`EW002.tex:13-19`), and T010 for `R_b`, `A_FB^b`, and `A_b` (`T010.tex:13-29`).
Plan v1 had T012/T013 for Z-charm pole observables and T022 for W/Z universality
(`flavor_catalog_plan_v1.md:315-326`), but the 75-row fact-check inventory has
no Z-charm or W/Y row (`factcheck_status.md:74-88`). Judgment: W/Y can be
deferred to a synthesis-level EW fit, but Z-charm, `Gamma_Z`, and W/Y should be
explicitly marked deferred; otherwise the EW scope looks narrower than advertised.

Two lower-priority B-sector tails are absent but reasonably deferred. Deep
Research points to `B_s -> tau tau` and `B -> K(*) tau tau`
(`deepresearch_may15.txt:912-932`), and says `B -> tau nu` should be interpreted
with `B+ -> mu nu`, `R(D)`, `R(D*)`, and `B_c -> tau nu` consistency
(`deepresearch_may15.txt:786-790`). Plan v1 listed B008/B010/B027/B028
(`flavor_catalog_plan_v1.md:269-290`), while DA-4 deferred "electron/tau rare
tails" and charged-current constituents already covered at overview level
(`round_004_convergence.md:35`). Judgment: that rationale is sound for v0.1,
but B009/B025/B026 should cross-reference these tails so the omission does not
look accidental.

Deep Research also flags semileptonic asymmetries `a_sl^d,s` as cross-checks of
Delta B = 2 phases (`deepresearch_may15.txt:885-895`). I found no matching
process by grep. This should not be a new standalone first-tier process, but it
belongs as a nuance in the B-mixing rows because it probes a phase combination
not identical to `S_psiKS` or `phi_s`.

## Section 2: Processes / constraints we have that the external research disagrees with

I found no BLOCKER-level value, sign, or interpretation disagreement where the
catalog is plainly wrong. The catalog's current value-bearing claims were already
fact-checked with 74 VERIFIED rows and one non-value PARTIAL row
(`factcheck_status.md:108-155`), and Deep Research is mostly a ranking memo
rather than an independent numerical audit.

The clearest WARNING is granularity in B017. The catalog says B017 is
`B^0 -> K* l+l-` with linked `b -> s l+l-` observables and records branching and
LFU-ratio context (`B017.tex:11-20`, `B017.tex:22-44`). It correctly states that
angular observables require form factors, nonlocal charm, bin definitions, and
correlated likelihoods (`B017.tex:69-75`). Deep Research, however, specifically
asks for the modern angular basis `P_i`, `F_L`, `A_FB`, `S_i`, `A_i` in standard
`q^2` bins (`deepresearch_may15.txt:774-780`) and says giant blocks are only
family-complete, not bin-exhaustive (`deepresearch_may15.txt:1013-1019`). Type:
granularity / validity range. Severity: WARNING, because B017 is not wrong, but
the catalog should make the missing angular-basis/covariance import more explicit.

T010/EW001/EW002 have a scope mismatch with Deep Research's custodial-EW framing.
T010 says it packages `R_b^0`, `A_FB^{0,b}`, and `A_b` (`T010.tex:13-29`) and
already emphasizes custodial protection and bottom embedding choices
(`T010.tex:31-40`, `T010.tex:54-62`). EW001 covers only `S,T,U`
(`EW001.tex:20-33`). Deep Research asks for the "complete set" including
`delta g_Lb`, `delta g_Rb`, `Z -> cc`, `Gamma_Z`, `W,Y`, and `mW`
(`deepresearch_may15.txt:844-855`). Type: model-dependence / global-fit scope.
Severity: WARNING, not because the recorded T010/EW001 values are wrong, but
because an RS custodial analysis needs a coupling-shift representation, not only
the pseudo-observable values.

There is one NIT where Deep Research is sloppier than the catalog: it lists
"h -> tc and t -> cH as a correlated pair" (`deepresearch_may15.txt:909-910`).
For a 125 GeV Higgs, on-shell `h -> t c` is kinematically forbidden. Catalog T007
correctly treats `t -> Hc` and associated single-top-plus-H production as the
relevant direct probes (`T007.tex:10-15`, `T007.tex:63-70`). Type:
interpretation / wording. Severity: NIT; no catalog physics fix is needed, though
T007 could mention that "Htc" coupling constraints often enter through top decay
and associated production, not Higgs decay to an on-shell top.

## Section 3: Subtleties / nuances the external research surfaces that we did NOT make explicit

Photon chirality in radiative B decays: RS right-handed dipoles make polarization
more diagnostic than the inclusive branching ratio (`deepresearch_may15.txt:897-901`).
Relevant entries: B011, plus deferred B012/B013/B014. B011 already says the
inclusive rate tests `C7,C8` dipoles (`B011.tex:46-53`) and that using it as a
hard RS constraint needs a dedicated Wilson pipeline (`B011.tex:66-73`).
Suggested wording: "Inclusive `B -> X_s gamma` constrains `C7,C8`, but RS
chirality tests require exclusive polarization and time-dependent CP observables
sensitive to `C7'`."

Custodial dependence of EW constraints: Deep Research repeatedly ties Zbb/EW
precision to custodial/non-custodial architecture (`deepresearch_may15.txt:575-616`).
Relevant entries: T010, EW001, EW002. T010 already states that custodial
protection determines whether Zbb is leading or diagnostic (`T010.tex:31-40`);
EW001 similarly warns that custodial protection and embeddings change the mass
translation (`EW001.tex:55-62`). Suggested wording: "When translating these
EWPOs to RS mass limits, quote results only after specifying custodial protection,
fermion embeddings, brane kinetic terms, and whether W/Y or vertex corrections
are included."

Up-aligned model hierarchy: Deep Research says D mixing and charm CPV become
central when kaons are engineered safe (`deepresearch_may15.txt:13-18`,
`deepresearch_may15.txt:420-467`). Relevant entries: C001, C002, C003, T001,
T002, T005, T006, T007. C001 frames D mixing as a test of whether alignment
suppresses both down- and up-sector neutral currents (`C001.tex:43-49`), and
T001 frames top FCNCs as complementary direct top-sector probes (`T001.tex:36-43`).
Suggested wording: "In down-aligned or kaon-protected RS variants, up-sector
FCNCs and top FCNCs can become leading rather than secondary diagnostics."

RK/RK* anomaly status: Deep Research says RK and RK* are now broadly SM-like and
should not be treated as an anomaly prior (`deepresearch_may15.txt:400-407`,
`deepresearch_may15.txt:994-1011`). Relevant entries: B017, B018, B019. The
catalog already says this in B018 and B019 (`B018.tex:50-61`,
`B019.tex:48-64`), but the synthesis should repeat it. Suggested wording:
"Use RK/RK* as precision LFU null tests; do not impose the pre-2023 anomaly
narrative as a prior."

Long-distance limitations: Deep Research cautions against Delta m_K standalone
fits, KL dimuon/electron modes as primary bounds, B -> piK as hard mass setters,
and integrated rare charm decays (`deepresearch_may15.txt:945-981`,
`deepresearch_may15.txt:997-1002`). Relevant entries: K002, K006, B032, C007.
The catalog has the same caveats in K002 (`K002.tex:58-65`), K006
(`K006.tex:51-58`), B032 (`B032.tex:73-79`), and C007 (`C007.tex:59-66`).
Suggested wording: "These rows are diagnostic or conservative-envelope inputs
unless a short-distance subtraction, amplitude model, or bin-specific recast is
supplied."

EDMs inside a flavor fit: Deep Research says leaving EDMs out can give a
misleading viability picture (`deepresearch_may15.txt:936-938`). Relevant
entries: E001, E004, E006, E008, E009, and cross-references from K/B/D/top
dipole rows. E004 and E008 already connect EDMs to CP-odd quark/gluon dipoles
outside the Delta F = 2 basis (`E004.tex:41-48`, `E008.tex:57-62`). Suggested
wording: "For anarchic CP phases, EDM rows should be reviewed alongside flavor
branching ratios, not as an optional appendix."

## Section 4: Overall assessment of the external research itself

Deep Research gets the broad hierarchy right. Its top tier matches the catalog's
strongest themes: `epsilon_K` and LR Delta S = 2 pressure (K001/K002),
`epsilon'/epsilon` (K003), rare `K -> pi nu nubar` (K004/K005), `b -> s gamma`
(B011), B mixing phases (B001-B004), Zbb/EW precision (T010/EW001/EW002),
muon CLFV (L001-L005), and EDMs (E001/E004/E006/E008/E009). It also agrees with
the catalog's conservative validity framing: K006 is long-distance limited
(`K006.tex:51-58`), B032 is a penguin diagnostic rather than a hard mass setter
(`B032.tex:73-79`), and B017 needs a global-fit interface rather than a
single-number cut (`B017.tex:69-75`).

What it gets wrong or states too loosely is smaller but real. The `h -> tc`
wording is kinematically wrong for the observed Higgs. It also uses broad source
domain labels instead of exact papers/values, so it is not a substitute for the
catalog's fact-check process. It sometimes treats direct searches and high-pT
tails as if they belong naturally beside low-energy observables; for RS this is
useful context, but benchmark dependence is severe and Deep Research itself
admits no single mass number is universal (`deepresearch_may15.txt:1006-1009`).

What it gets right that the catalog under-covered: photon-polarization/radiative
B chirality, top photon FCNCs, EW W/Y and Z-charm/global-width pieces, explicit
direct-search/high-pT context, and semileptonic asymmetries in B mixing. The top
photon omission is the most concrete process-list miss. The radiative-B
polarization omission is the most RS-specific physics miss.

Both documents still have gaps. Neither is a machine-readable likelihood. The
catalog is process-complete in prose but not bin/covariance complete for
`b -> s l+l-`, not a global EDM fit with cancellations, and not a direct-search
reinterpretation framework. Deep Research is useful as an external checklist,
but it is too coarse to arbitrate numerical values or implementation policy.
The practical next step is not to reopen all discovery, but to add a short
deferred-scope/synthesis addendum covering the concrete misses above, with T003/T004
top-photon FCNCs treated as the only clear candidate for a new process row.

===EXTERNAL_REVIEW_deepresearch_may15_END===
