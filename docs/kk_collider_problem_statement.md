# Making the KK collider bounds rigorous in our warped (RS) flavor scan — problem statement

## Setup
We run a large parameter scan of a Randall-Sundrum warped extra dimension with
anarchic 5D Yukawas and bulk fermions localized by dimensionless masses c, a
brane-localized Higgs, and an IR/KK scale M_KK. For each point we impose flavor
and electroweak constraints, and we also want collider limits on the Kaluza-Klein
(KK) resonances. Right now the collider bounds are a proxy: we compare each point
against a generic, rescaled cross-section-times-branching limit rather than
computing the rate from that point's own RS couplings. We want to replace this
with a defensible, per-point recast of real LHC searches.

## The states we care about
- KK gluon (color octet, first KK level), the dominant strong resonance, decaying
  mostly to top pairs and broad in width.
- KK electroweak gauge bosons (the KK photon, KK Z, and KK W), decaying to top
  pairs, dijets, and diboson final states.
- KK graviton, decaying to diphoton, ZZ, WW.
- Vector-like quark partners (the custodial top and bottom partners), decaying to
  a third-generation quark plus a W, Z, or Higgs.

We are quark-sector-focused, so dilepton channels are out of scope for now; the
relevant final states are top pairs, dijets, diboson, and the vector-like-quark
channels.

## What is point-dependent (and why a generic limit is not enough)
In this framework the couplings of the KK resonances to light quarks are set by
the same bulk localization (the c parameters) that controls the quark masses and
mixings. Light quarks are localized away from the IR brane, so their coupling to
the KK states is suppressed and roughly universal, while the third generation is
IR-localized and strongly coupled. As a result both the production cross section
(driven by the light-quark couplings and the parton luminosities) and the decay
branching fractions (dominated by tops) vary point to point. The KK gluon in
particular is broad, with a width that can be tens of percent of its mass, so the
narrow-width approximation that most published limits assume does not directly
apply.

## The questions we need answered
1. Channel scope. Which searches should we recast to bound these states: top-pair
   resonances for the KK gluon and KK electroweak states, dijet resonances for the
   light-quark channel, diboson and diphoton for the KK graviton, and
   vector-like-quark pair and single production for the partners? Which are the
   binding ones, and which can we safely drop?

2. Production cross section. How should we compute the production rate of each KK
   state from the point's actual RS couplings (light-quark coupling squared folded
   with the parton luminosity), rather than scaling a single benchmark cross
   section? What is the right level of approximation: leading order with a K-factor,
   or is more needed?

3. Width and line shape. The KK gluon (and to a lesser extent the others) is broad.
   How do we model the broad-width line shape, and how does a width of tens of
   percent change the acceptance relative to the narrow benchmark the experiments
   quote? When does the narrow-width approximation fail badly enough that we must
   use the full line shape?

4. Branching fractions. How do we compute the per-point branching into the searched
   final state (for the KK gluon, the split between top pairs, light jets, and
   bottom pairs), given the point's couplings?

5. Acceptance and efficiency. Published limits fold in an analysis-specific
   acceptance times efficiency tuned to a benchmark model. How do we obtain the
   acceptance appropriate to our coupling structure and width so the recast is
   honest, rather than borrowing the benchmark acceptance?

6. Interference. How important is interference, both between the KK state and its
   Standard Model partner (for example KK-gluon and gluon interference in the
   top-pair spectrum) and between signal and background for a broad resonance, and
   do we need to include it to get the bound right?

7. Which limits, and combination. Which specific published searches (collider
   energy, integrated luminosity, final state) should we map each point onto, and
   how do we combine channels into a single exclusion while stating the
   model-dependence honestly?

## What has changed since 2007 (current state, as of 2026)
The foundational reference is Lillie, Randall, Wang, "The Bulk RS KK-gluon at the
LHC" (hep-ph/0701166, JHEP 0709:074). That was a pre-data discovery study; three
things have changed since.

1. Boosted-top reconstruction is now standard. LRW had to invent techniques for
   tagging the highly boosted tops from a heavy resonance; jet substructure and
   boosted-top tagging are now routine, so the experiments directly measure the
   acceptance for a heavy, broad g_KK to boosted top pairs.

2. We now have real, width-binned limits instead of a projection. Full Run-2
   (13 TeV, 140 fb^-1) top-pair resonance searches quote limits as a function of
   resonance width, which is exactly the broad-width issue here. Representative
   numbers: CMS excludes a leptophobic topcolor Z' up to about 3.8, 5.25, and
   6.65 TeV for widths of 1, 10, and 30 percent, and an RS KK gluon up to about
   4.55 TeV; the latest ATLAS full-Run-2 search excludes a 30-percent-width g_KK
   below about 4.1 TeV and a spin-2 KK graviton below about 1.3 TeV. HL-LHC
   projections reach about 5.7 TeV (discovery) to 6.6 TeV (exclusion).

3. The model framing shifted. It is now understood that minimal anarchic RS is
   pushed to tens of TeV by flavor and electroweak data, not by colliders. Our
   corrected scan finds a TYPICAL ~30 TeV minimal floor from `epsilon_K` (flavor)
   and an irreducible EXISTENCE floor ~18-20 TeV from oblique S,T,U; the
   post-B1-fix Z->bb constraint is only ~5 TeV. Collider limits (~4 TeV) sit far
   below all of these. The field moved to custodial RS,
   composite Higgs / partial compositeness, and "little RS" with brane-localized
   kinetic terms, which can bring the scale down toward 5 TeV.

Why this matters here: for minimal RS the collider limit (about 4 TeV) is far below
the flavor/electroweak floor (tens of TeV), so it never binds. But our custodial
scan puts the floor at about 7 TeV (inclusive) and 2-3 TeV (rigorous), right in the
window where the 13 TeV top-pair limits bite. So the collider recast is a
co-leading constraint precisely in the custodial regime, which is the reason to
make it rigorous now rather than later.

## RS KK-gluon model inputs (from Lillie-Randall-Wang)
These specify the per-point g_KK rate and width (questions 2-4) in terms of the
bulk localization c (our convention c = M_5/k, c > 1/2 UV-peaked / light,
c < 1/2 IR-peaked):
- Coupling of the first KK gluon to a bulk quark, in units of the QCD coupling g_s:
  about -0.2 g_s for UV-localized light quarks (c > 1/2, roughly universal and
  suppressed); about 1 g_s for the third-generation doublet (c approximately 0.4);
  about 4 g_s for the IR-localized right-handed top (c approximately 0.3). The
  coupling is a smooth function of c interpolating between these.
- Total width Gamma/M is about 0.17 for a typical configuration, with branching to
  top pairs about 92.5 percent (the rest mostly to light jets and bottom pairs).
- Production is q qbar to g_KK through the suppressed light-quark coupling, so the
  cross section scales as the sum over light quarks of (g_q/g_s)^2 times the q qbar
  parton luminosity, and is point-dependent through the c values.
- Mass relation: the first KK gluon mass is m ≈ x1 * Lambda with x1 ≈ 2.45 (our
  physical-mass convention, distinct from the geometric Lambda_IR).

Per-point recipe: from each point's c values compute g_q(c) for every quark, hence
the production cross section (light-quark-coupling-squared times parton luminosity),
the width Gamma/M, and the top-pair branching fraction; then test cross-section-
times-branching against the width-binned Run-2 boosted top-pair limit at the
matching width bin (the broad ~15-30 percent bin for the KK gluon), using the
experiment's measured boosted-top acceptance.

## What we will do with the answers
Wire a per-point collider module that, for each scan point, computes the
production cross section and branching fractions of the KK states from its RS
couplings, applies the appropriate width and acceptance, and tests the resulting
cross-section-times-branching against the chosen published limits, replacing the
current generic proxy with a recast we can defend as rigorous.
