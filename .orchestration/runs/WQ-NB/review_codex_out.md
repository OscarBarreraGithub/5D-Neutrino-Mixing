nbconvert: exit 0; 33 cells, 12 PNG outputs/fig files, 0 traceback/error outputs.
4-TeV floor: FAIL. A4 input panel samples `df[~df.skipped]`, leaking 264,568 non-skipped 1/2/3 TeV rows; C2 y-axis lower limit is 3.5 TeV.
No-scale-mixing: FAIL. A1/A2/A3/A4 pool over `M_KK`; A4/C3b pool over `r`; C3a bins RH-compositeness across `r`.
Tags/advisories: per-row tag fields are invariant, not hardcoded: 9 rigorous / 12 proxy / 16 partial; inclusive survival matches HARD rigorous/proxy veto accounting exactly, mismatch 0.
Advisory isolation: partials are excluded from carving/survival and isolated in C1c “DOES NOT CARVE M_KK”; wording fix: T010 is `partial/HARD`, not INFO/SOFT.
Fit spot-check: r=0.05,5TeV cQ/up_sv/dn_sv calc=stored, maxdiff 0; r=0.25,10TeV maxdiff 0; r=1,50TeV maxdiff 0.
Headline recompute: eligible peak median ratios EW001=1.67805, B012=0.369326 next; EW001 veto=1.0 at 5TeV and no eligible veto at >=7TeV.
Survival: strict=1.0 for all retained tiles; inclusive=0 at 5TeV and 1.0 at 7/10/15/20/30/50TeV.
Hierarchy/top column: cQ medians descend gen1>gen2>gen3 and cross 0.5 in all exact cells; top-column norm medians decrease 4.842->4.130->3.470->3.153->2.954.
Quark-only labeling: present; 37 allowlisted constraints listed by physics name, and lepton sector dropped/not fully rigorous is stated.
Required fixes: filter A4 input to retained `M_KK`, facet/group all distributions by exact `(r,M_KK)` or remove pooled panels, set M_KK axes floor to 4, and correct T010 INFO/SOFT wording.
WQ-NB-REVIEW: NEEDS-FIXES