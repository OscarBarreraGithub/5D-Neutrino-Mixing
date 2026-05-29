1. NIT: No amplitude defect. E007 computes no `M12^NP` and no RS EDM amplitude; it records direct `|d_A|` limits only, with `predicted=None` and no RS CP-odd matching flags. Future live EDM physics would need CP-odd Wilsons/imaginary phases before Schiff+atomic response. `flavor_catalog_constraints/primary/edm_neutrino/E007.py:403`

2. NIT: QCD-running check is N/A, not failed. No low-energy CP-odd Wilsons are produced, so there is no non-running verdict path; the stub explicitly refuses RS CP-odd matching. `E007.py:420`, `flavor_catalog_constraints/physics_adapters/atomic_edm.py:48`

3. NIT: Budget/anchors are defensible for non-vetoing bookkeeping: scalar budget is Xe direct limit `1.4e-27 e cm`; diagnostics also retain Ra `1.4e-23`, Ra 2015 `5.0e-22`, independent Xe `1.5e-27`, and Xe central/uncertainties `1.4 +/- 6.6_stat +/- 2.0_syst e-28`. Not a central residual. `E007.py:147`, `flavor_catalog/processes/edm_neutrino/E007.yaml:99`

4. NIT: Severity/units are appropriate: `Severity.INFO`, `passes=True`, `non_vetoing=True`, units consistently `e cm`, and both nuclear/atomic plus RS CP-odd sides are flagged `NEEDS-HUMAN-PHYSICS`. `E007.py:376`, `E007.py:420`

PHYSICS-OK