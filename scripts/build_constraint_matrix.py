#!/usr/bin/env python3
"""Build a compact per-draw pass/fail matrix from a scan output dir.

Reads the per-draw JSONL rows produced by the quark/full-catalog scan and
extracts, for every draw, a small fixed set of columns:

    M_KK_TeV   : first KK resonance mass in TeV  (params["M_KK"] / 1000)
    r          : quark-fit residual tolerance    (quark_fit_r / params["quark_fit_r"])
    fit_success: bool, the per-draw quark-fit precondition
                 (fit_diagnostics["success"]; "Reproduce SM masses" + "Reproduce CKM")
    skipped    : bool, row["skipped"] -- the point builder bailed on this draw
                 BEFORE evaluating constraints (e.g. non-perturbative lepton
                 seesaw Yukawa). Such rows carry an EMPTY constraints dict.
    skip_reason: str, row["skip_reason"] (e.g. "nonperturbative_lepton_yukawa",
                 "quark_fit_failed"); empty string when not skipped.
    evaluated  : bool, DERIVED -- True iff constraints were actually run for this
                 draw, i.e. (constraints dict is non-empty) AND NOT skipped. This
                 is the CRUCIAL flag: ~95% of full-catalog draws are SKIPPED with
                 an empty constraints dict, and a skipped draw must NOT be counted
                 as "passes every constraint". All survival/veto/floor statistics
                 in the explorer notebook must restrict to `evaluated` rows.
    pass_<ID>  : bool, constraints[ID]["passes"]   for each ID in CONSTRAINT_SET
                 (None/NA for skipped rows and for IDs absent from the row)
    ratio_<ID> : float, constraints[ID]["ratio"]   for each ID in CONSTRAINT_SET
    pass_COLLIDER : DERIVED bool, M_KK >= COLLIDER_THRESHOLD_GEV (geometric;
                 defined for every row regardless of `evaluated`)

The whole point: physics is evaluated ONCE here. The explorer notebook then
only does boolean AND-masks over these precomputed columns when constraints are
toggled -- it never re-runs any physics. ALL rows are kept (including skipped
ones); the `evaluated` flag is how downstream code restricts to the physically
correct population.

Output: <scan_dir>/constraint_matrix.parquet  (parquet if pyarrow present,
else a numpy .npz fallback at constraint_matrix.npz).

CLI:
    python scripts/build_constraint_matrix.py <scan_dir>
        [--collider-gev 5500] [--out PATH] [--max-tiles N] [--limit N]

Robust to missing constraint IDs: any ID in CONSTRAINT_SET that never appears in
the data is reported as a warning and its pass_/ratio_ columns are filled with
NaN/pd.NA (and excluded from the discovered-ID manifest).
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

# --------------------------------------------------------------------------
# CONFIG -- edit here. Maps the user's "review_local" constraint set to the
# scan's internal constraint IDs. Auto-discovery warns about any that are
# absent from the actual data (e.g. L001 is full-catalog-only).
# --------------------------------------------------------------------------
CONSTRAINT_SET = {
    "K001": "epsilon_K",
    "C001": "D0-D0bar mixing",
    "C002": "D0-D0bar CPV",
    "B001": "Delta m_d",
    "B003": "Delta m_s",
    "L001": "mu -> e gamma  (full-catalog only)",
    "T010": "Z -> bb (coupling)",
    "T011": "Z -> bb (aux)",
    "EW001": "S, T, U oblique",
}

# Collider cut: pass iff M_KK >= this many GeV. m1 = 5.5 TeV by default; the
# scan's M_KK variable already IS the physical first KK resonance mass.
COLLIDER_THRESHOLD_GEV = 5500.0


# --------------------------------------------------------------------------
def find_jsonl_files(scan_dir):
    """Find per-draw JSONL tiles under a scan dir.

    Supports two layouts seen in the wild:
      <dir>/r*/shard-*/tile-*.jsonl        (completed runs)
      <dir>/tile-*.jsonl[.tmp.*]           (in-progress flat layout)
    """
    patterns = [
        os.path.join(scan_dir, "r*", "shard-*", "tile-*.jsonl"),
        os.path.join(scan_dir, "r*", "tile-*.jsonl"),
        os.path.join(scan_dir, "tile-*.jsonl"),
        os.path.join(scan_dir, "tile-*.jsonl.tmp.*"),
    ]
    files = []
    for p in patterns:
        files.extend(glob.glob(p))
    # de-dup, keep deterministic order
    return sorted(set(files))


def iter_rows(files, limit=None):
    n = 0
    for f in files:
        try:
            with open(f) as fh:
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        yield json.loads(line)
                    except json.JSONDecodeError:
                        # truncated final line of an in-progress tmp file
                        continue
                    n += 1
                    if limit is not None and n >= limit:
                        return
        except OSError:
            continue


def build_records(files, collider_gev, limit=None):
    ids = list(CONSTRAINT_SET.keys())
    seen_ids = set()
    records = []
    n_total = 0
    for row in iter_rows(files, limit=limit):
        n_total += 1
        params = row.get("params", {})
        mkk_gev = params.get("M_KK")
        if mkk_gev is None:
            continue
        fd = row.get("fit_diagnostics") or {}
        cons = row.get("constraints", {}) or {}
        skipped = bool(row.get("skipped", False))
        skip_reason = row.get("skip_reason") or ""
        # A draw is "evaluated" iff constraints actually ran for it: the
        # constraints dict is non-empty AND the point builder did not bail
        # (skipped is not True). Skipped draws carry an empty constraints dict
        # and must be excluded from all survival/veto/floor statistics.
        evaluated = bool(cons) and not skipped
        rec = {
            "M_KK_TeV": float(mkk_gev) / 1000.0,
            "M_KK_GeV": float(mkk_gev),
            "r": row.get("quark_fit_r", params.get("quark_fit_r")),
            "fit_success": bool(fd.get("success", True)),
            "skipped": skipped,
            "skip_reason": str(skip_reason),
            "evaluated": evaluated,
        }
        for cid in ids:
            c = cons.get(cid)
            if c is None:
                rec[f"pass_{cid}"] = None
                rec[f"ratio_{cid}"] = float("nan")
            else:
                seen_ids.add(cid)
                rec[f"pass_{cid}"] = bool(c.get("passes"))
                ratio = c.get("ratio")
                rec[f"ratio_{cid}"] = float(ratio) if ratio is not None else float("nan")
        # derived collider cut
        rec["pass_COLLIDER"] = bool(rec["M_KK_GeV"] >= collider_gev)
        records.append(rec)
    return records, seen_ids, n_total


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("scan_dir", help="scan output directory")
    ap.add_argument("--collider-gev", type=float, default=COLLIDER_THRESHOLD_GEV,
                    help=f"collider M_KK threshold in GeV (default {COLLIDER_THRESHOLD_GEV})")
    ap.add_argument("--out", default=None, help="output path (default <scan_dir>/constraint_matrix.parquet)")
    ap.add_argument("--max-tiles", type=int, default=None, help="cap number of tile files read (for dev)")
    ap.add_argument("--limit", type=int, default=None, help="cap number of rows read (for dev)")
    args = ap.parse_args()

    scan_dir = os.path.abspath(args.scan_dir)
    if not os.path.isdir(scan_dir):
        print(f"ERROR: not a directory: {scan_dir}", file=sys.stderr)
        sys.exit(1)

    files = find_jsonl_files(scan_dir)
    if args.max_tiles is not None:
        files = files[: args.max_tiles]
    if not files:
        print(f"ERROR: no JSONL tiles found under {scan_dir}", file=sys.stderr)
        sys.exit(1)
    print(f"[build] {len(files)} tile file(s) under {scan_dir}")

    records, seen_ids, n_total = build_records(files, args.collider_gev, limit=args.limit)
    print(f"[build] read {n_total} rows, kept {len(records)} draws with M_KK")

    # auto-discovery report
    requested = set(CONSTRAINT_SET.keys())
    missing = sorted(requested - seen_ids)
    present = sorted(seen_ids)
    print(f"[discover] present IDs ({len(present)}): {present}")
    if missing:
        print(f"[discover] WARNING missing IDs (filled with NA): {missing}")
        for cid in missing:
            print(f"            - {cid}: {CONSTRAINT_SET[cid]}")

    out = args.out or os.path.join(scan_dir, "constraint_matrix.parquet")

    # try pandas+pyarrow, else numpy npz fallback
    wrote = None
    try:
        import pandas as pd
        df = pd.DataFrame.from_records(records)
        # pass_ columns -> nullable boolean so missing IDs stay distinguishable
        for cid in CONSTRAINT_SET:
            col = f"pass_{cid}"
            if col in df.columns:
                df[col] = df[col].astype("boolean")
        df["pass_COLLIDER"] = df["pass_COLLIDER"].astype("boolean")
        df["fit_success"] = df["fit_success"].astype("boolean")
        for bcol in ("skipped", "evaluated"):
            if bcol in df.columns:
                df[bcol] = df[bcol].astype("boolean")
        # attach metadata via a sidecar (parquet schema metadata is awkward in pandas 3)
        try:
            df.to_parquet(out, index=False)
            wrote = out
        except Exception as e:  # pyarrow missing or write error
            print(f"[write] parquet failed ({e}); falling back to .npz")
            wrote = None
        if wrote is None:
            out = os.path.splitext(out)[0] + ".npz"
            _write_npz(out, df)
            wrote = out
    except ImportError:
        print("[write] pandas unavailable; writing numpy .npz fallback")
        out = os.path.splitext(out)[0] + ".npz"
        _write_npz_from_records(out, records)
        wrote = out

    # sidecar manifest (machine-readable discovery report)
    manifest = {
        "scan_dir": scan_dir,
        "n_draws": len(records),
        "collider_threshold_gev": args.collider_gev,
        "constraint_set": CONSTRAINT_SET,
        "present_ids": present,
        "missing_ids": missing,
        "columns": sorted(records[0].keys()) if records else [],
    }
    man_path = os.path.splitext(wrote)[0] + ".manifest.json"
    with open(man_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"[write] wrote matrix -> {wrote}")
    print(f"[write] wrote manifest -> {man_path}")


def _write_npz(out, df):
    import numpy as np
    arrays = {}
    for col in df.columns:
        s = df[col]
        if str(s.dtype) == "boolean":
            # store as int8 with -1 for NA
            v = s.astype("object")
            arrays[col] = np.array([(-1 if x is None or (x is not False and x is not True and x != x)
                                     else int(bool(x))) for x in v], dtype=np.int8)
        else:
            arrays[col] = s.to_numpy()
    np.savez_compressed(out, **arrays)


def _write_npz_from_records(out, records):
    import numpy as np
    cols = sorted(records[0].keys())
    bool_cols = {"fit_success", "skipped", "evaluated"}
    arrays = {}
    for col in cols:
        vals = [r.get(col) for r in records]
        if col.startswith("pass_") or col in bool_cols:
            arrays[col] = np.array([(-1 if v is None else int(bool(v))) for v in vals], dtype=np.int8)
        elif col == "skip_reason":
            arrays[col] = np.array(["" if v is None else str(v) for v in vals], dtype=object)
        else:
            arrays[col] = np.array([(float("nan") if v is None else v) for v in vals], dtype=float)
    np.savez_compressed(out, **arrays)


if __name__ == "__main__":
    main()
