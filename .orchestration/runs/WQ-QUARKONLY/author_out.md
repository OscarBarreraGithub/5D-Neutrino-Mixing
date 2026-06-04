Plan done: verified candidates, added quark-only path, added allowlist tallies/provenance, added tests, ran smoke and full pytest.

Changed: `--quark-only` / `ScanConfig.quark_only`; full-mode config payload omits false flag to preserve hash surface. No commit made.

Allowlist table: `K001,K002,B001-B004,C001,C002 -> qmbc -> N -> IN -> ΔF=2 quark couplings only`.
`B011-B014,T001-T004,T007,T008 -> qmbc+kk_ew -> N -> IN -> quark/top radiative/EW/H proxy inputs only`.
`T005,T006 -> qmbc+kk_gluon -> N -> IN -> top gluon inputs only`.
`T010,T011,T012,T014 -> rs_ew_couplings+qmbc -> N -> IN -> quark Z shifts only; no lepton Z fields read`.
`EW001 -> kk_ew/kk_gluon/qmbc -> N -> IN -> oblique mass proxy, no swept leptons`.
`EW003 -> optional rs_charged_current diagnostic -> N -> IN -> scalar is data-only; CC not built/read for veto`.
`B032-B034,C003,K003,K013,E004,E006-E009 -> none -> N -> IN -> INFO/data stubs use no point extras`.
`EW002 -> rs_charged_current -> Y -> OUT -> needs charged-current e/mu epsilon/G_F object; tagged deferred_lepton_followup`.

Build flags in quark-only: `include_charged_current=False`, `include_fermion_kk_mixing=True`, `include_higgs_yukawas=False`, `lepton_yukawa_result=None`.

Smoke 1k: 999/1000 evaluated, 37/37 allowlist evaluated+active per point, 0 constraint exceptions, 0.2786 s/evaluated point, 1e5 extrapolation 7.74 core-hours.

Veto ranking: `EW001 750 proxy`, `B003 250 rigorous`, `B004 250 rigorous`, `B011 250 proxy`, `B012 250 proxy`, `K001 250 rigorous`, `B013 138 proxy`.

M_KK sanity: same quark fit `K001` ratio 1.659 at 1.5 TeV to 0.0179 at 15 TeV; `K002` 0.0391 to 0.000413.

Tests: `tests/test_full_catalog_scan_harness.py` 11 passed; full `python -m pytest tests/ -q` 1711 passed, 1 skipped.

WQ-AUTHOR-DONE.