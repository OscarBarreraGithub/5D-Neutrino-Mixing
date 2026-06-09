# W7-MUEGAMMA implementation plan

Scope: implementation plan only. No production code is written in this run. The locked model is Perez-Randall LMFV, using the spurion `(Y_N_bar Y_N_bar^dagger)_{12}` and the existing `flavorConstraints.muToEGamma` dipole machinery.

## Hard-veto contract

The catalog L001 veto is driven entirely by the predicted branching fraction:

```text
passes = (BR_NP <= br_limit)
BR_NP = prefactor_br * |(Y_N_bar Y_N_bar^dagger)_{e mu}|^2 * (reference_scale_gev / M_KK)^4
br_limit = 1.5e-13
```

The injected `C=0.02` (`C_PAPER`) is an inert diagnostic for `dipole_rhs` and `dipole_ratio_to_bound`; it must not decide pass/fail in the catalog path. `C=0.02` is paired with the paper-era Perez-Randall/MEGA limit `1.2e-11`, while the MEG-II-consistent coefficient is `sqrt(1.5e-13 / 4.0e-8) = 0.001936491673...`, already documented as the live/default coefficient in `flavor_catalog/processes/charged_lepton/L001.tex:40-48` and `L001.yaml:254-275`.

Keep these three ratios separate because they measure different quantities:

- `BR_NP / br_limit = 1033.8851067579139`: the catalog hard-veto ratio.
- `lhs / rhs = 3.1133057793694858` with `C=0.02`: a paper-C dipole-bound diagnostic.
- `lhs / rhs = 32.154083827064859` with `C=0.001936491673...`: a MEG-II-consistent dipole-bound diagnostic.

Do not edit `flavor_catalog/processes/charged_lepton/L001.tex` or `L001.yaml` for W7. The existing catalog text describes the generic bound and the legacy/live coefficient relationship. W7 adds an explicit LMFV-prediction path in code while keeping that catalog prose separate. If the final catalog should report only the BR-form veto, both the generic C-bound text and the LMFV prediction path, that remains a catalog-presentation question to flag rather than silently resolving in this implementation.

## Re-verified repository facts

- `flavorConstraints/muToEGamma.py` already has the reusable dipole core: `PREFAC_BR=4.0e-8` at line 31, paper-era `BR_LIMIT_PAPER=1.2e-11` at line 33, `C_PAPER=0.02` at line 36, `coefficient_from_br_limit` at lines 61-73, `check_mu_to_e_gamma` at lines 75-146, and `check_mu_to_e_gamma_raw` at lines 149-198. The core returns `core["passes"] = lhs <= rhs`; W7 must treat that as diagnostic only in the catalog adapter.
- `flavor_catalog_constraints/physics_adapters/lepton.py` already converts the core result into a branching-fraction result. `_result_from_core` computes `branching_fraction` and `ratio_to_limit`, then sets `passes=ratio_to_limit <= 1.0` at lines 176-196. Current line 259 derives `c_lfv` from the BR limit; W7 needs an optional fixed `c_lfv` so L001 can emit the paper-C diagnostic without changing the BR veto.
- `flavor_catalog_constraints/primary/charged_lepton/L001.py` currently requires `"lepton_mass_basis_couplings"` at lines 50-52 and calls the adapter at lines 280-288 without passing `c_lfv`. W7 changes L001 to require the new LMFV extra and to pass `c_lfv=self.anchor.c_paper.value`, while keeping pass/fail on `BR_NP <= br_limit`.
- `flavor_catalog/processes/charged_lepton/L001.yaml` records the MEG II 2025 limit `1.5e-13` at lines 151-165, the paper-era `c_paper=0.02` at lines 218-229, and the repo-default `lfv_C=0.0019364916731037085` plus `c_paper=0.02` at lines 230-275. `L001.tex:40-48` states the same live/default-vs-paper distinction; do not overwrite that text.
- The legacy lepton scan in `scanParams/scan.py` is separate from catalog W7: it computes `lfv_C = coefficient_from_br_limit(config.br_limit, prefactor=config.prefac_br)` at line 595 and passes that coefficient into `check_mu_to_e_gamma` at lines 523-534. W7 should not change this legacy behavior.
- `YukawaResult` carries `Y_N`, `Y_N_bar`, `Y_N_matrix`, `epsilon`, and `params` in `yukawa/compute_yukawas.py:32-82`. `compute_all_yukawas` constructs PMNS with `get_pmns(...)`, computes `Y_N_matrix = V_pmns @ diag(Y_N)` at lines 283-287, and stores `Lambda_IR`, `c_L`, `c_E`, `c_N`, `M_N`, phases, `k`, and `v` at lines 301-313.
- `RSLeptonMassBasisCouplings` already carries many LMFV ingredients in `quarkConstraints/rs_ew_couplings.py:98-126`, including `kk_ew_mass_gev`, `Y_N_bar_vector`, `Y_N_matrix`, `Y_N_bar_matrix`, `pmns`, `lfv_dipole_spurion`, and `params`. Its builder validates `Y_N_bar_matrix == 2*k*Y_N_matrix == PMNS@diag(Y_N_bar)` at lines 388-417 and records `params["M_KK"] = kk_ew_mass_gev` at lines 422-424.
- `build_rs_ew_extras` resolves optional lepton Yukawas at `flavor_catalog_constraints/rs_ew_builder.py:102-107`, builds `lepton_mass_basis_couplings` at lines 108-115, and returns that extra only when lepton Yukawas are available at lines 167-179. Full scan passes `lepton_yukawa_result` at `scripts/run_full_catalog_scan.py:555-585`; quark-only mode passes `lepton_yukawa_result=None` at lines 465-479 and evaluates only `QUARK_ONLY_ALLOWLIST_IDS`.
- `scripts/run_full_catalog_scan.py:75-82` centralizes quark-only forbidden extras. Add the new LMFV key there for documentation and regression coverage. Do not add `ScanConfig` fields or alter `_config_payload`/`_config_hash` at lines 1486-1502.
- The harness tags any evaluated result with truthy diagnostic keys containing `"proxy"` or `"recast"` as proxy at `scripts/run_full_catalog_scan.py:955-989` and `_proxy_flags` at lines 1045-1057. The generated LMFV carrier path must emit no truthy proxy/recast-named diagnostic; `used_proxy=False` is safe because `_proxy_flags` ignores false values.

## Carrier decision

Use a dedicated `LMFVLeptonParameters` extra under `"lepton_lmfv_parameters"` rather than making L001 read `RSLeptonMassBasisCouplings` directly.

Reason: `RSLeptonMassBasisCouplings` is a broader EW/tree-current carrier in `quarkConstraints/rs_ew_couplings.py`, which is shared with W8. W7 only needs a stable dipole-spurion contract and should avoid extending or depending semantically on that W8 carrier. The duplication is intentional and bounded: the new carrier is a catalog-side, read-only LMFV view derived from the same `YukawaResult` and physical `base_spectrum.kk_ew_mass_gev`; it records `epsilon`, which the existing carrier does not expose, and validates only the fields needed for L001.

Do not add an L001 fallback from `"lepton_lmfv_parameters"` to `"lepton_mass_basis_couplings"`. A fallback with `used_proxy=False` can be classified rigorous by the harness even when it is only transitional. If adapter-level compatibility for `RSLeptonMassBasisCouplings` is kept, route it through the same validation helper and record `extra_used`/`input_kind`, but L001 itself should require the new extra so the rigorous path is unambiguous.

## Implementation steps

1. Add the LMFV carrier and declared extra key.

   - Add `"lepton_lmfv_parameters"` to `KNOWN_EXTRA_KEYS` in `flavor_catalog_constraints/point_builder.py`.
   - Define immutable `LMFVLeptonParameters` in `flavor_catalog_constraints/physics_adapters/lepton.py`.
   - Required fields: `Y_N`, `Y_N_bar`, `Y_N_matrix`, `Y_N_bar_matrix`, `pmns`, `lmfv_spurion`, `M_KK_gev`, `M_N_gev`, `c_L`, `c_E`, `c_N`, `v_gev`, `k_gev`, `epsilon`, `Lambda_IR_gev`, `ordering`, `majorana_alpha`, `majorana_beta`, `max_abs_ybar`, `source`, and a safe `matching_status` string such as `"locked_lmfv_nda_carrier"`.
   - The constructor must copy arrays to read-only `np.complex128`/`np.float64`, require shapes `(3,)` or `(3,3)`, require positive finite masses/scales and finite `epsilon`, verify PMNS unitarity, verify `Y_N_bar_matrix == 2*k*Y_N_matrix == pmns @ diag(Y_N_bar)`, and verify `lmfv_spurion == Y_N_bar_matrix @ Y_N_bar_matrix.conj().T`.
   - The carrier must not impose a new physics prior. Scan perturbativity remains `_require_perturbative_leptons` in `scripts/run_full_catalog_scan.py:832-838`; the carrier only validates finiteness and records diagnostics.

2. Attach the carrier in the existing builder.

   - In `flavor_catalog_constraints/rs_ew_builder.py`, after `resolved_lepton_yukawas` and `base_spectrum` are available, build `lmfv_lepton_parameters_from_yukawa_result(resolved_lepton_yukawas, m_kk_gev=base_spectrum.kk_ew_mass_gev)`.
   - Return `{"lepton_lmfv_parameters": lmfv}` only when `resolved_lepton_yukawas is not None`, parallel to the existing `lepton_mass_basis_couplings` conditional.
   - Add `"lepton_lmfv_parameters"` to `QUARK_ONLY_FORBIDDEN_EXTRAS` in `scripts/run_full_catalog_scan.py:75-82`. Quark-only mode remains unchanged because it passes `lepton_yukawa_result=None`; this is a documentation/guardrail edit, not a config-hash edit.

3. Update the adapter without changing `flavorConstraints/muToEGamma.py`.

   - Extend `mu_to_e_gamma_from_lepton_input(...)` with `c_lfv: float | None = None`.
   - Preserve compatibility: if `c_lfv is None`, keep the current derived `sqrt(br_limit/prefactor_br)` behavior so existing L002/L003/L004 reuse paths and legacy proxy tests continue to work.
   - Add a coercion path for `LMFVLeptonParameters`, using `Y_N_bar`, `pmns`, and `M_KK_gev` with `check_mu_to_e_gamma_raw(..., C=c_lfv, reference_scale=3000.0)`.
   - Always compute catalog pass/fail from `BR_NP <= br_limit`, independent of `core["passes"]` and independent of `c_lfv`. Keep `core_passes`, `dipole_rhs`, and `dipole_ratio_to_bound` only as diagnostics.
   - For the generated carrier path set `input_kind="LMFVLeptonParameters"`, `extra_used="lepton_lmfv_parameters"`, and `used_proxy=False`. Do not emit truthy keys or values with `"proxy"`/`"recast"` in their names; avoid `needs_human_physics` on this locked-model path.
   - Keep explicit `MuToEGammaProxyInput` and loose mapping inputs marked `used_proxy=True` with the existing `needs_human_physics` diagnostic.

4. Wire L001 to the locked LMFV prediction path.

   - Change `_REQUIRED_EXTRA` from `"lepton_mass_basis_couplings"` to `"lepton_lmfv_parameters"`.
   - If the new extra is absent, preserve the non-vetoing unevaluated semantics: `passes=True`, `predicted=None`, `ratio=None`, `diagnostics["evaluated"]=False`, and `diagnostics["missing_extra"] == "lepton_lmfv_parameters"`.
   - Do not fall back to `"lepton_mass_basis_couplings"` in L001. Existing tree/contact charged-lepton constraints can continue using the legacy extra.
   - Call the adapter with `br_limit=self.anchor.budget`, `prefactor_br=self.anchor.prefactor_br.value`, `c_lfv=self.anchor.c_paper.value`, and `reference_scale_gev=3000.0`.
   - Prefer the carrier's `M_KK_gev`. If `kk_ew_mass_gev` is also present on the point, check/report consistency rather than silently overriding the carrier mass.
   - Diagnostics for the generated carrier path should include `evaluated=True`, `lmfv_model="Perez-Randall LMFV NDA"`, `lfv_coefficient=0.02`, `c_lfv_role="dipole_rhs_diagnostic_only"`, `br_limit=1.5e-13`, `br_ratio_to_limit`, `dipole_lhs`, `dipole_rhs`, `dipole_ratio_to_bound`, `core_passes`, `branching_formula`, `off_diagonal_12`, `product_matrix`, `input_kind`, `extra_used`, and `used_proxy=False`.
   - Result classification should be `rigorous` under the locked LMFV model. Do not place `needs_human_physics`, truthy `proxy` diagnostics, truthy `recast` diagnostics, or `matching_status` text containing `partial`, `deferred`, `proxy`, or `recast` on this path.

## Numeric oracle

Use the repo-local benchmark in `yukawa/compute_yukawas.py:220-238` and `tests/test_mu_to_e_gamma.py:58-99` for a direct adapter/carrier oracle. This oracle is not a builder mass-convention test.

```text
compute_all_yukawas(
  Lambda_IR=3000.0,
  c_L=0.58,
  c_E=[0.75, 0.60, 0.50],
  c_N=0.27,
  M_N=1.22e18,
  lightest_nu_mass=0.002,
  ordering="normal",
)
M_KK = 3000.0 GeV  # explicit adapter/carrier oracle mass
C = 0.02           # diagnostic only
br_limit = 1.5e-13

Y_N_bar = [0.20416916, 0.43091265, 1.02237364]
off_diagonal_12 = -0.034773700046830024 + 0.05165132075170267j
lhs = 0.062266115587389717
rhs(C=0.02) = 0.02
dipole_ratio_to_bound(C=0.02) = 3.1133057793694858
BR_NP = 1.5508276601368708e-10
BR_NP / br_limit = 1033.8851067579139
passes = False

C_consistent = 0.0019364916731037084
dipole_ratio_to_bound(C_consistent) = 32.154083827064859
```

The oracle must assert that `predicted == BR_NP`, `ratio == BR_NP/br_limit`, and `passes is False` because the BR ratio is above one. It may also assert the two dipole diagnostic ratios, but it must not describe those diagnostics as gating the catalog result.

Add a separate builder mass-convention test. Build a point with `Lambda_IR` derived from a target physical mass (`lambda_ir = mkk / xi`, using the existing exact-root helper where available, e.g. `tests/rs_ew_phase3b_helpers._scales_for_mkk`). Assert the carrier mass equals `point.extras["kk_ew_mass_gev"] == base_spectrum.kk_ew_mass_gev`, not raw `Lambda_IR`. Do not reuse the direct oracle's `Lambda_IR=3000.0, M_KK=3000.0` assumptions for this builder test.

## Tests to add or update

- `tests/constraints/primary/charged_lepton/test_L001.py`
  - Update absent-extra and invalid-extra tests to expect `"lepton_lmfv_parameters"`.
  - Add an end-to-end generated-carrier test using the direct oracle above. Assert `predicted`, `ratio`, `passes`, `dipole_lhs`, `dipole_rhs`, `dipole_ratio_to_bound`, `lfv_coefficient == 0.02`, `c_lfv_role`, `used_proxy is False`, and `tag_result(result)[0] == "rigorous"`.
  - Add a pass/fail independence test: evaluate the same carrier with `c_lfv=0.02` and a very large `c_lfv` (large enough to flip `core_passes` if desired). Assert `dipole_rhs` changes, but `branching_fraction`, `ratio_to_limit`, and `passes` are identical and pinned to the BR ratio.
  - Add a no-proxy-flag test for the generated carrier diagnostics: every diagnostic key containing `proxy` or `recast` must have value `False` or `None`; `tag_result` must return no proxy flags.
  - Move legacy proxy/mapping coverage to adapter-level assertions or keep it clearly separate from the L001 rigorous path; those inputs must remain `proxy` if `used_proxy=True`.
  - Add malformed-carrier tests for non-finite arrays, non-unitary PMNS, zero/negative `M_KK_gev`, and inconsistent `Y_N_bar_matrix`; L001 should return unevaluated with `invalid_extra`.
- `tests/test_mu_to_e_gamma.py`
  - Keep the existing core tests for paper `C=0.02`, derived MEG-II `C`, and explicit `M_KK` override.
  - Add or tighten the benchmark assertions above for raw/core values. Phrase them as core diagnostics plus BR oracle, not as two independent catalog gates.
- `tests/test_lmfv_lepton_parameters.py` or `tests/test_rs_ew_phase4a.py`
  - Assert `build_from_rs_ew_inputs(..., lepton_sweep_inputs=...)` returns both `lepton_mass_basis_couplings` and `lepton_lmfv_parameters`.
  - Assert arrays are read-only, `epsilon` is finite, `M_KK_gev` is the physical `kk_ew_mass_gev`, and `lmfv_spurion == Y_N_bar_matrix @ Y_N_bar_matrix.conj().T`.
  - Add the mass-convention split test described above.
- `tests/constraints/test_contract.py`
  - Add declared-key coverage for `point_builder.make_point(lepton_lmfv_parameters=marker)`.
- `tests/test_full_catalog_scan_harness.py`
  - Add `"lepton_lmfv_parameters"` to the expected quark-only forbidden extras behavior.
  - Strengthen the deterministic quark-only draw test to compare exact serialized rows, using the harness writer settings: `json.dumps(harness._json_sanitize(row), sort_keys=True, separators=(",", ":")).encode()`. Compare two deterministic draws byte-for-byte and against an inline/golden expected serialized row, not just config hash, `params`, or selected fields.
  - Assert quark-only rows still have `lepton_sector == "dropped (not rigorous)"`, no `lepton_inputs`, and no `lepton_lmfv_parameters` in builder extras.
  - Preserve the current config-hash surface; the reference hashes for the two-mass test config remain `full=45e21a07585f7489` and `quark_only=d96cb734f724aedb` unless unrelated existing tests prove those values are already stale.
- Perturbativity/finite guards:
  - Unit-test carrier rejection of non-finite `Y_N_bar`, non-unitary PMNS, zero/negative `M_KK_gev`, and inconsistent `Y_N_bar_matrix`.
  - Keep nonperturbative lepton rejection in `_require_perturbative_leptons`; add a scan-style regression showing nonperturbative `Y_N_bar` is skipped before L001 evaluation.

Suggested verification commands after implementation:

```bash
pytest tests/constraints/primary/charged_lepton/test_L001.py \
       tests/test_mu_to_e_gamma.py \
       tests/test_lmfv_lepton_parameters.py \
       tests/test_rs_ew_phase4a.py \
       tests/constraints/test_contract.py \
       tests/test_full_catalog_scan_harness.py
```

## Isolation and serialization

Files planned to touch in the future implementation:

- `flavor_catalog_constraints/physics_adapters/lepton.py`: add `LMFVLeptonParameters`, builder/coercion support, optional fixed `c_lfv`, and BR-pinned pass/fail tests.
- `flavor_catalog_constraints/point_builder.py`: declare `"lepton_lmfv_parameters"`.
- `flavor_catalog_constraints/rs_ew_builder.py`: attach the new extra when lepton Yukawas are available.
- `flavor_catalog_constraints/primary/charged_lepton/L001.py`: require the new extra and call the adapter with `C=0.02` as diagnostic-only input plus `br_limit=1.5e-13` as the hard veto.
- `scripts/run_full_catalog_scan.py`: add only the new forbidden-extra key for quark-only documentation/guardrails; do not alter `ScanConfig`, `_config_payload`, `_config_hash`, or row schema except through tests that prove byte identity is preserved.
- Test files listed above.

Files not planned to touch:

- `flavorConstraints/muToEGamma.py`: core already implements the reusable dipole diagnostic and must remain unchanged.
- `flavor_catalog/processes/charged_lepton/L001.yaml` and `L001.tex`: keep the existing generic-bound/live-default description separate from the new LMFV prediction path.
- `scanParams/scan.py`: legacy lepton scan keeps its derived-`C` behavior.
- W8 shared files `quarkConstraints/rs_ew_couplings.py` and `quarkConstraints/oblique_stu.py`: do not edit for W7.

Serialization notes:

- `flavor_catalog_constraints/point_builder.py`, `flavor_catalog_constraints/rs_ew_builder.py`, and `scripts/run_full_catalog_scan.py` are used by the W9 scan harness, so W7 should serialize with W9 if W9 is editing builder integration or scan serialization at the same time.
- The new extra must not alter quark-only rows, quark-only config hashes, or quark-only allowlist semantics.

## Deferred / NEEDS-HUMAN scope

- Full loop-level RS dipole matching is not implemented. W7 uses the locked Perez-Randall LMFV NDA formula already present in `flavorConstraints/muToEGamma.py`; it does not claim a new loop calculation.
- Non-LMFV/anarchic/Chen-Yu/gauged-LFV models remain out of scope by locked decision.
- Arbitrary off-diagonal charged-lepton rotation effects beyond the LMFV spurion carrier remain out of scope for L001. Tree-level LFV contacts for other charged-lepton constraints continue to live in the existing `rs_ew_couplings` path.
- Heavy neutral exchange, boxes, scalar pieces, and relative phases remain deferred for the other LFV observables as documented in their adapters; W7 should not alter them.

PLAN-READY
