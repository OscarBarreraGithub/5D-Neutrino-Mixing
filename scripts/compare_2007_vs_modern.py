#!/usr/bin/env python3
"""
Compare 2007-era (Fitzpatrick-Perez-Randall, 0710.1869) vs modern (2024+)
exclusion bounds.

SUPERSEDED/CORRECTED: the former implementation post-hoc rescaled modern scan
ratios to "2007-era" ratios. That is invalid for B_d/B_s/D0: the scan rows
store hadronic |M12|/budget ratios in GeV, while the old rescale factors used
legacy dimensionless operator-weight bounds. The prior "9.4x D0 tightening",
"tightening driver", and acceptance-delta figures were units artifacts.

Usage:
    python scripts/compare_2007_vs_modern.py \
        scan_outputs/dense_20260414T125549/merged/results.jsonl \
        --output-dir results/figures/analysis
"""

import argparse
import json
import os
import textwrap
from collections import Counter, defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Supersession metadata
# ---------------------------------------------------------------------------

SUPERSEDED_TITLE = 'SUPERSEDED/CORRECTED: 2007-vs-modern all-system rescale'

SUPERSEDED_NOTICE = (
    'The former 2007-vs-modern comparison is retracted. Modern scan rows store '
    'ratios against hadronic |M12| budgets in GeV for B_d, B_s, and D0. The '
    'old rescale used legacy dimensionless operator-weight bounds for those '
    'systems, so it multiplied incompatible quantities. The prior 9.4x D0 '
    'tightening claim, tightening-driver plot, and 2007-vs-modern acceptance '
    'delta were units artifacts. No all-system 2007 ratio or acceptance '
    'comparison is computed by this script.'
)

INVALID_LEGACY_OUTPUTS = (
    'fig_2007_vs_modern_exclusion.pdf',
    'fig_2007_vs_modern_mkk_bound.pdf',
    'fig_2007_vs_modern_driver.pdf',
    'fig_2007_vs_modern_acceptance.pdf',
)

SYSTEM_IDS = ['epsilon_K', 'K', 'B_d', 'B_s', 'D0']

SYSTEM_LABELS = {
    'epsilon_K': r'$\epsilon_K$',
    'K': r'$\Delta M_K$',
    'B_d': r'$B_d$ mixing',
    'B_s': r'$B_s$ mixing',
    'D0': r'$D^0$ mixing',
}

SYSTEM_COLORS = {
    'epsilon_K': '#e41a1c',
    'K': '#377eb8',
    'B_d': '#4daf4a',
    'B_s': '#984ea3',
    'D0': '#ff7f00',
}


def load_data(path):
    """Load JSONL scan results, return list of dicts with parsed fields."""
    records = []
    with open(path) as f:
        for line in f:
            d = json.loads(line)
            ratios = d['ratio_to_bound_by_system']
            ratios_modern = {s: float(ratios[s]) for s in SYSTEM_IDS}
            max_ratio_modern = max(ratios_modern.values())
            ratio_pass_modern = max_ratio_modern <= 1.0
            verifier_ok = all(
                bool(d.get(field, True))
                for field in (
                    'verifier_ok',
                    'bridge_verifier_ok',
                    'phenomenology_verifier_ok',
                )
            )
            fit_success = bool(d.get('fit_success', d.get('fit_converged', True)))
            rec = {
                'r': d['r'],
                'M_KK': d['M_KK'],
                'overall_scale': d['overall_scale'],
                'accepted_modern': bool(
                    d.get(
                        'accepted',
                        fit_success and verifier_ok and ratio_pass_modern,
                    )
                ),
                'ratio_pass_modern': ratio_pass_modern,
                'fit_success': fit_success,
                'verifier_ok': verifier_ok,
                'fit_score': d.get('fit_score', None),
                'ratios_modern': ratios_modern,
                'max_ratio_modern': max_ratio_modern,
                'binding_modern': max(ratios_modern, key=ratios_modern.get),
            }
            records.append(rec)
    return records


def aggregate_best_scale(records):
    """
    At each (r, M_KK) grid point, pick the overall_scale that minimizes
    the max ratio (the 'best' point). Return arrays for plotting.
    """
    grid = defaultdict(list)
    for rec in records:
        key = (rec['r'], rec['M_KK'])
        grid[key].append(rec)

    r_vals, mkk_vals = [], []
    max_ratio_modern = []
    accepted_modern = []
    ratio_pass_modern = []
    binding_modern = []

    for (r, mkk), recs in grid.items():
        # Best overall_scale for modern
        best_mod = min(recs, key=lambda x: x['max_ratio_modern'])

        r_vals.append(r)
        mkk_vals.append(mkk)
        max_ratio_modern.append(best_mod['max_ratio_modern'])
        accepted_modern.append(any(rec['accepted_modern'] for rec in recs))
        ratio_pass_modern.append(any(rec['ratio_pass_modern'] for rec in recs))
        binding_modern.append(best_mod['binding_modern'])

    return {
        'r': np.array(r_vals),
        'mkk': np.array(mkk_vals),
        'max_ratio_modern': np.array(max_ratio_modern),
        'accepted_modern': np.array(accepted_modern),
        'ratio_pass_modern': np.array(ratio_pass_modern),
        'binding_modern': binding_modern,
    }


def write_superseded_markdown(records, agg, outdir):
    """Write an explicit replacement report for the retracted comparison."""
    path = os.path.join(outdir, 'compare_2007_vs_modern_SUPERSEDED.md')
    n = len(records)
    n_mod = sum(1 for r in records if r['accepted_modern'])
    n_ratio = sum(1 for r in records if r['ratio_pass_modern'])
    n_grid = len(agg['r'])
    n_grid_mod = int(agg['accepted_modern'].sum())

    lines = [
        f'# {SUPERSEDED_TITLE}',
        '',
        f'> {SUPERSEDED_NOTICE}',
        '',
        '## What remains valid',
        '',
        '- The input scan rows can be summarized under the modern scan convention.',
        '- Modern acceptance is the stored scan gate: fit success, verifier gates, and ratio gates.',
        '- The stored B_d, B_s, and D0 ratios are hadronic |M12|/budget quantities in GeV.',
        '',
        '## What is not computed',
        '',
        '- No all-system 2007 ratio is reconstructed from these scan rows.',
        '- No 2007-vs-modern acceptance delta is reported.',
        '- No tightening-driver figure is reported.',
        '',
        '## Modern scan summary',
        '',
        f'- Total scan points: {n:,}',
        f'- Accepted by stored modern gate: {n_mod:,} ({100*n_mod/max(n, 1):.2f}%)',
        f'- Ratio-only pass under modern ratios: {n_ratio:,} ({100*n_ratio/max(n, 1):.2f}%)',
        f'- Unique (r, M_KK) grid points: {n_grid:,}',
        f'- Grid points with at least one stored-modern accepted row: {n_grid_mod:,}',
        '',
        (
            'The difference between ratio-only pass and stored acceptance is expected: '
            'the modern scan acceptance also includes fit and verifier gates.'
        ),
        '',
    ]
    with open(path, 'w') as fh:
        fh.write('\n'.join(lines))
    print(f'  Saved: {path}')


def plot_superseded_notices(outdir):
    """Overwrite legacy figure names with an explicit supersession banner."""
    wrapped = textwrap.fill(SUPERSEDED_NOTICE, width=86)
    for filename in INVALID_LEGACY_OUTPUTS:
        fig, ax = plt.subplots(figsize=(9.0, 5.2))
        ax.axis('off')
        ax.text(
            0.5,
            0.72,
            SUPERSEDED_TITLE,
            ha='center',
            va='center',
            fontsize=15,
            fontweight='bold',
        )
        ax.text(
            0.5,
            0.47,
            wrapped,
            ha='center',
            va='center',
            fontsize=10.5,
            linespacing=1.35,
        )
        ax.text(
            0.5,
            0.18,
            'This file intentionally contains no 2007-vs-modern comparison plot.',
            ha='center',
            va='center',
            fontsize=10,
            style='italic',
        )
        path = os.path.join(outdir, filename)
        fig.savefig(path, bbox_inches='tight', dpi=150)
        print(f'  Saved: {path}')
        plt.close(fig)


def print_summary(records, agg):
    """Print summary statistics to stdout."""
    n = len(records)
    n_mod = sum(1 for r in records if r['accepted_modern'])
    n_ratio = sum(1 for r in records if r['ratio_pass_modern'])

    print('\n' + '=' * 70)
    print(f'  {SUPERSEDED_TITLE}')
    print('=' * 70)
    print('  ' + textwrap.fill(SUPERSEDED_NOTICE, width=66).replace('\n', '\n  '))
    print()
    print(f'  Total scan points:       {n:>8,}')
    print(f'  Accepted (modern gate):  {n_mod:>8,}  ({100*n_mod/max(n, 1):.2f}%)')
    print(f'  Ratio-only pass:         {n_ratio:>8,}  ({100*n_ratio/max(n, 1):.2f}%)')
    print()

    binding_modern = Counter()
    for rec in records:
        if not rec['accepted_modern']:
            binding_modern[rec['binding_modern']] += 1
    print('  Binding system among modern-gate rejected rows:')
    for s in SYSTEM_IDS:
        if binding_modern[s] > 0:
            print(f'    {s:>12s}: {binding_modern[s]:>6,}')
    print()

    n_grid = len(agg['r'])
    n_grid_mod = agg['accepted_modern'].sum()
    n_grid_ratio = agg['ratio_pass_modern'].sum()
    print(f'  Grid points (best overall_scale): {n_grid:,}')
    print(f'    Modern accepted:     {n_grid_mod:>5}')
    print(f'    Modern ratio-pass:   {n_grid_ratio:>5}')
    print('=' * 70)
    print()


def main():
    parser = argparse.ArgumentParser(
        description='Supersede the invalid 2007-era vs modern quark comparison',
    )
    parser.add_argument('results_file', help='Path to merged results.jsonl')
    parser.add_argument('--output-dir', default='results/figures/analysis',
                        help='Directory for output figures')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print(f'Loading data from {args.results_file} ...')
    records = load_data(args.results_file)
    print(f'  Loaded {len(records):,} points')

    print('Aggregating best overall_scale per (r, M_KK) grid point ...')
    agg = aggregate_best_scale(records)
    print(f'  {len(agg["r"]):,} unique (r, M_KK) grid points')

    print_summary(records, agg)

    print('Writing superseded/corrected outputs ...')
    write_superseded_markdown(records, agg, args.output_dir)
    plot_superseded_notices(args.output_dir)

    print('\nDone.')


if __name__ == '__main__':
    main()
