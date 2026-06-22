# Constraint Explorer notebook — regenerate & open

Interactive view of how each constraint bites on the RS minimal-model M_KK reach,
toggleable on/off over a PRECOMPUTED per-draw pass/fail matrix (no physics
recompute on toggle).

## Files

- `scripts/build_constraint_matrix.py` — JSONL scan dir -> compact `constraint_matrix.parquet`
- `notebooks/constraint_explorer.ipynb` — the interactive notebook (executed)
- `notebooks/_constraint_explorer_src.py` — editable source for the notebook
- `notebooks/_build_notebook.py` — assembles the .ipynb from the source

## 1. Build the matrix

```bash
cd /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing
# real fix100k run (once it finishes):
python scripts/build_constraint_matrix.py scan_outputs/fix100k_minimal_<TS>
# -> writes scan_outputs/fix100k_minimal_<TS>/constraint_matrix.parquet (+ .manifest.json)
```

Options: `--collider-gev 5500` (collider M_KK cut, TeV*1000), `--out PATH`,
`--max-tiles N` / `--limit N` (dev subsampling). Robust to missing IDs
(L001 is full-catalog-only; quark-only runs warn and fill it with NA).
Falls back to a numpy `.npz` if pyarrow is unavailable.

Matrix columns: `M_KK_TeV`, `M_KK_GeV`, `r`, `fit_success` (the SM-mass+CKM
precondition), `pass_<ID>` + `ratio_<ID>` for each of
{K001, C001, C002, B001, B003, L001, T010, T011, EW001}, and derived
`pass_COLLIDER` (M_KK >= threshold).

## 2. Open / re-run the notebook

```bash
# rebuild .ipynb from source after editing _constraint_explorer_src.py:
python notebooks/_build_notebook.py
# execute headless to refresh outputs:
jupyter nbconvert --to notebook --execute --inplace \
  --ExecutePreprocessor.timeout=600 notebooks/constraint_explorer.ipynb
```

The notebook auto-resolves its input to the newest `fix100k_minimal_*` run that
has a built matrix, else the schema-dev quark-only run. Override with
`CONSTRAINT_MATRIX_PATH=/abs/path/constraint_matrix.parquet`. **Only the input
path changes when the real fix100k data swaps in.**

Interactive use (needs `pip install ipywidgets` in the kernel): launch Jupyter,
toggle the per-constraint checkboxes. Without ipywidgets the notebook still runs
and renders the static all-ON view plus the comparison renders.

## Design: precompute-then-toggle (the key requirement)

Physics is evaluated ONCE by the build script. On any checkbox toggle the
notebook does only `survivor = AND over selected pass_<ID> columns` — a single
boolean reduce over precomputed columns. No constraint is re-evaluated and no
point is re-run; we only ask which already-evaluated draws survive the selected
cut subset.

## What it shows (validated against the schema-dev quark-only run, 240k draws)

- Per-constraint own M_KK floor (lowest M_KK where its veto fraction < 0.5):
  **T010 (Z->bb) = 20 TeV (rank 1)**, EW001 (S,T,U) = 7 TeV.
- Combined floor (survival > 0.5), all ON = **20 TeV**, set by Z->bb.
- Drop Z->bb -> combined floor falls to **7 TeV**.
- So at ~20-25 TeV the binding constraint is Z->bb, NOT S,T,U — the ranking
  inversion is explicit in panel (b) and the summary table.
