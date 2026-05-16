# WA Worklog: wa_w5b_kaon

Family: kaon

Processes: K009, K010

## K009 -- \(K_L \to \pi^0\mu^+\mu^-\)

- Tightened the process, RS-relevance, post-2008, and validity language while keeping the PKA-documented physics content.
- Normalized Section B headings and added explicit source-file tags for the PDG/KTeV limit, NA48 \(K_S\) support input, lattice/Snowmass context, and Isidori--Smith--Unterdorfer numerical decomposition.
- Appended `WRITER-DONE` to the YAML `status_history` and updated `last_updated_at`.

Checker notes:
- Re-query PDG `S013.16` before freeze, as already flagged by the PKA.
- Keep the supporting \(K_S\to\pi^0\mu^+\mu^-\) value synchronized with the separate process entry that carries it canonically.

## K010 -- \(K_S \to \pi^0 e^+e^-\)

- Tightened the process and relevance wording, removing filler while preserving the measured partial-rate versus extrapolated-total distinction.
- Normalized Section B headings and made the numerical claims point back to the YAML `pdg_or_equivalent` and supporting NA48 source snapshot entries.
- Appended `WRITER-DONE` to the YAML `status_history` and updated `last_updated_at`.

Checker notes:
- Confirm whether downstream catalog plotting should default to the directly measured \(m_{ee}>0.165~\mathrm{GeV}\) partial branching fraction or the model-extrapolated full-rate value.

No `\textbf{CHECK}` markers were introduced.  No shared `catalog.bib` patch was made because this batch was restricted to process files and the writer worklog; process-local source keys remain listed in the entries.
