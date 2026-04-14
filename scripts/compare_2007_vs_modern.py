#!/usr/bin/env python3
"""
Compare 2007-era (Fitzpatrick-Perez-Randall, 0710.1869) vs modern (2024+)
exclusion bounds by post-hoc rescaling of modern scan results.

The key insight: ratio_to_bound = effective_amplitude / bound.  If the bound
was looser in 2007, the same amplitude gives a smaller ratio:

    ratio_2007 = ratio_modern * (bound_modern / bound_2007)

Usage:
    python scripts/compare_2007_vs_modern.py \
        scan_outputs/dense_20260414T125549/merged/results.jsonl \
        --output-dir results/figures/analysis
"""

import argparse
import json
import os
import sys
from collections import Counter, defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from scipy.interpolate import griddata

# ---------------------------------------------------------------------------
# 2007 vs modern bound rescaling factors
# ---------------------------------------------------------------------------
# ratio = bound_modern / bound_2007
# If < 1, modern bound is tighter (the usual case).
RESCALE = {
    'epsilon_K': 4.18e-4 / 6.0e-4,    # ~0.697
    'K':         1.742e-15 / 2.5e-15,  # ~0.697
    'B_d':       4.0e-7 / 6.0e-7,      # ~0.667
    'B_s':       5.5e-6 / 1.5e-5,      # ~0.367
    'D0':        8.5e-9 / 8.0e-8,      # ~0.106
}

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
            rec = {
                'r': d['r'],
                'M_KK': d['M_KK'],
                'overall_scale': d['overall_scale'],
                'accepted_modern': d['accepted'],
                'fit_score': d.get('fit_score', None),
            }
            ratios = d['ratio_to_bound_by_system']
            rec['ratios_modern'] = {s: ratios[s] for s in SYSTEM_IDS}

            # 2007 rescaling
            ratios_2007 = {}
            for s in SYSTEM_IDS:
                ratios_2007[s] = ratios[s] * RESCALE[s]
            rec['ratios_2007'] = ratios_2007

            rec['max_ratio_modern'] = max(rec['ratios_modern'].values())
            rec['max_ratio_2007'] = max(rec['ratios_2007'].values())
            rec['accepted_2007'] = rec['max_ratio_2007'] <= 1.0

            # Binding system (which has the largest ratio)
            rec['binding_modern'] = max(rec['ratios_modern'], key=rec['ratios_modern'].get)
            rec['binding_2007'] = max(rec['ratios_2007'], key=rec['ratios_2007'].get)

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
    max_ratio_modern, max_ratio_2007 = [], []
    accepted_modern, accepted_2007 = [], []
    binding_modern, binding_2007 = [], []
    # For the "driver" figure: which system tightened the most
    tightening_driver = []

    for (r, mkk), recs in grid.items():
        # Best overall_scale for modern
        best_mod = min(recs, key=lambda x: x['max_ratio_modern'])
        # Best overall_scale for 2007
        best_2007 = min(recs, key=lambda x: x['max_ratio_2007'])

        r_vals.append(r)
        mkk_vals.append(mkk)
        max_ratio_modern.append(best_mod['max_ratio_modern'])
        max_ratio_2007.append(best_2007['max_ratio_2007'])
        accepted_modern.append(best_mod['accepted_modern'])
        accepted_2007.append(best_2007['accepted_2007'])
        binding_modern.append(best_mod['binding_modern'])
        binding_2007.append(best_2007['binding_2007'])

        # Driver: at the best modern point, which system's ratio increased
        # the most going modern -> 2007?  I.e., which system's tightening
        # contributed the most?  We compare ratio_modern / ratio_2007 = 1/RESCALE.
        # But for the actual point, we look at the modern best point and ask
        # which system was most helped by the 2007 loosening:
        # tightening_factor[s] = ratio_modern[s] / ratio_2007[s] = 1 / RESCALE[s]
        # But since this is constant, let's instead look at which system's
        # absolute ratio change (modern - 2007) is largest at this grid point.
        # That captures which system drives the actual exclusion change.
        rec = best_mod
        biggest_sys = None
        biggest_diff = -np.inf
        for s in SYSTEM_IDS:
            diff = rec['ratios_modern'][s] - rec['ratios_2007'][s]
            if diff > biggest_diff:
                biggest_diff = diff
                biggest_sys = s
        tightening_driver.append(biggest_sys)

    return {
        'r': np.array(r_vals),
        'mkk': np.array(mkk_vals),
        'max_ratio_modern': np.array(max_ratio_modern),
        'max_ratio_2007': np.array(max_ratio_2007),
        'accepted_modern': np.array(accepted_modern),
        'accepted_2007': np.array(accepted_2007),
        'binding_modern': binding_modern,
        'binding_2007': binding_2007,
        'tightening_driver': tightening_driver,
    }


def plot_exclusion_comparison(agg, outdir):
    """Side-by-side exclusion contour: 2007 vs modern."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True)

    log_r = np.log10(agg['r'])
    log_mkk = np.log10(agg['mkk'] / 1e3)  # M_KK in TeV

    # Create regular grid for interpolation
    ri = np.linspace(log_r.min(), log_r.max(), 200)
    mi = np.linspace(log_mkk.min(), log_mkk.max(), 200)
    RI, MI = np.meshgrid(ri, mi)

    vmin, vmax = -1.5, 2.5

    for ax, ratio_key, title in [
        (ax1, 'max_ratio_2007', '2007-era bounds (FPR 0710.1869)'),
        (ax2, 'max_ratio_modern', 'Modern bounds (2024+)'),
    ]:
        log_ratio = np.log10(np.clip(agg[ratio_key], 1e-10, None))
        ZI = griddata(
            np.column_stack([log_r, log_mkk]),
            log_ratio,
            (RI, MI),
            method='cubic',
        )

        cmap = plt.cm.RdYlGn_r
        im = ax.pcolormesh(
            ri, mi, ZI,
            cmap=cmap, vmin=vmin, vmax=vmax, shading='auto',
        )
        # Contour at ratio = 1 (log10 = 0)
        cs = ax.contour(ri, mi, ZI, levels=[0.0], colors='black', linewidths=2)
        ax.clabel(cs, fmt='ratio=1', fontsize=9)

        ax.set_xlabel(r'$\log_{10}\,r$', fontsize=14)
        ax.set_title(title, fontsize=13)

    ax1.set_ylabel(r'$\log_{10}(M_{\mathrm{KK}}\,/\,\mathrm{TeV})$', fontsize=14)

    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    cb = fig.colorbar(im, cax=cbar_ax)
    cb.set_label(r'$\log_{10}(\max\; \mathrm{ratio\_to\_bound})$', fontsize=12)

    fig.suptitle('Allowed Region: 2007 vs 2024+ Bounds', fontsize=15, y=1.02)

    path = os.path.join(outdir, 'fig_2007_vs_modern_exclusion.pdf')
    fig.savefig(path, bbox_inches='tight', dpi=150)
    print(f'  Saved: {path}')
    plt.close(fig)


def plot_mkk_bound(agg, outdir):
    """M_KK lower bound vs r for 2007 and modern."""
    r_unique = np.sort(np.unique(agg['r']))

    mkk_min_modern = []
    mkk_min_2007 = []

    for r_val in r_unique:
        mask = agg['r'] == r_val

        # Modern: find min M_KK where accepted
        acc_mask = mask & agg['accepted_modern']
        if acc_mask.any():
            mkk_min_modern.append(agg['mkk'][acc_mask].min() / 1e3)
        else:
            mkk_min_modern.append(np.nan)

        # 2007
        acc_mask = mask & agg['accepted_2007']
        if acc_mask.any():
            mkk_min_2007.append(agg['mkk'][acc_mask].min() / 1e3)
        else:
            mkk_min_2007.append(np.nan)

    mkk_min_modern = np.array(mkk_min_modern)
    mkk_min_2007 = np.array(mkk_min_2007)

    fig, ax = plt.subplots(figsize=(8, 5.5))

    ax.plot(r_unique, mkk_min_modern, 'b-', lw=2.5, label='Modern (2024+)', zorder=3)
    ax.plot(r_unique, mkk_min_2007, 'r--', lw=2.5, label='2007-era (FPR)', zorder=3)

    # Shade between them where both exist
    valid = np.isfinite(mkk_min_modern) & np.isfinite(mkk_min_2007)
    if valid.any():
        ax.fill_between(
            r_unique[valid],
            mkk_min_2007[valid],
            mkk_min_modern[valid],
            alpha=0.25, color='purple',
            label='Constraint tightening',
        )

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$r = m_{\mathrm{KK}}^{(1)} / M_{\mathrm{KK}}$', fontsize=14)
    ax.set_ylabel(r'Minimum $M_{\mathrm{KK}}$ [TeV]', fontsize=14)
    ax.set_title(r'$M_{\mathrm{KK}}$ Lower Bound: 2007 vs Modern', fontsize=14)
    ax.legend(fontsize=12, loc='upper right')
    ax.grid(True, alpha=0.3, which='both')
    ax.tick_params(labelsize=12)

    path = os.path.join(outdir, 'fig_2007_vs_modern_mkk_bound.pdf')
    fig.savefig(path, bbox_inches='tight', dpi=150)
    print(f'  Saved: {path}')
    plt.close(fig)


def plot_tightening_driver(agg, outdir):
    """Color-code (r, M_KK) plane by which system drove the biggest tightening."""
    fig, ax = plt.subplots(figsize=(8, 5.5))

    log_r = np.log10(agg['r'])
    log_mkk = np.log10(agg['mkk'] / 1e3)

    # Map system names to integers
    sys_to_int = {s: i for i, s in enumerate(SYSTEM_IDS)}
    driver_int = np.array([sys_to_int[s] for s in agg['tightening_driver']])

    # Create regular grid
    ri = np.linspace(log_r.min(), log_r.max(), 200)
    mi = np.linspace(log_mkk.min(), log_mkk.max(), 200)
    RI, MI = np.meshgrid(ri, mi)

    ZI = griddata(
        np.column_stack([log_r, log_mkk]),
        driver_int.astype(float),
        (RI, MI),
        method='nearest',
    )

    colors = [SYSTEM_COLORS[s] for s in SYSTEM_IDS]
    cmap = mcolors.ListedColormap(colors)
    bounds = np.arange(-0.5, len(SYSTEM_IDS) + 0.5, 1)
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    im = ax.pcolormesh(ri, mi, ZI, cmap=cmap, norm=norm, shading='auto')

    cb = fig.colorbar(im, ax=ax, ticks=range(len(SYSTEM_IDS)))
    cb.ax.set_yticklabels([SYSTEM_LABELS[s] for s in SYSTEM_IDS], fontsize=11)
    cb.set_label('Dominant tightening driver', fontsize=12)

    ax.set_xlabel(r'$\log_{10}\,r$', fontsize=14)
    ax.set_ylabel(r'$\log_{10}(M_{\mathrm{KK}}\,/\,\mathrm{TeV})$', fontsize=14)
    ax.set_title('Which system drove the largest constraint tightening?', fontsize=13)
    ax.tick_params(labelsize=12)

    path = os.path.join(outdir, 'fig_2007_vs_modern_driver.pdf')
    fig.savefig(path, bbox_inches='tight', dpi=150)
    print(f'  Saved: {path}')
    plt.close(fig)


def plot_acceptance_bar(records, outdir):
    """Bar chart of acceptance rates and binding-system breakdown."""
    n_total = len(records)
    n_accepted_modern = sum(1 for r in records if r['accepted_modern'])
    n_accepted_2007 = sum(1 for r in records if r['accepted_2007'])

    # For rejected points, count which system is binding
    binding_modern = Counter()
    binding_2007 = Counter()
    for rec in records:
        if not rec['accepted_modern']:
            binding_modern[rec['binding_modern']] += 1
        if not rec['accepted_2007']:
            binding_2007[rec['binding_2007']] += 1

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    # Left panel: acceptance rates
    ax = axes[0]
    eras = ['2007-era', 'Modern (2024+)']
    accepted = [n_accepted_2007, n_accepted_modern]
    rejected = [n_total - n_accepted_2007, n_total - n_accepted_modern]
    pct_acc = [100 * a / n_total for a in accepted]

    bars_acc = ax.bar(eras, accepted, color=['#2ca02c', '#2ca02c'], alpha=0.8,
                      label='Accepted')
    bars_rej = ax.bar(eras, rejected, bottom=accepted,
                      color=['#d62728', '#d62728'], alpha=0.6,
                      label='Rejected')

    for i, (a, p) in enumerate(zip(accepted, pct_acc)):
        ax.text(i, a / 2, f'{a:,}\n({p:.1f}%)', ha='center', va='center',
                fontsize=12, fontweight='bold', color='white')

    ax.set_ylabel('Number of scan points', fontsize=13)
    ax.set_title('Acceptance rate comparison', fontsize=13)
    ax.legend(fontsize=11, loc='upper right')
    ax.tick_params(labelsize=12)

    # Right panel: binding system breakdown for rejected points
    ax = axes[1]
    x_pos = np.arange(len(SYSTEM_IDS))
    width = 0.35

    counts_2007 = [binding_2007.get(s, 0) for s in SYSTEM_IDS]
    counts_modern = [binding_modern.get(s, 0) for s in SYSTEM_IDS]

    bars1 = ax.bar(x_pos - width / 2, counts_2007, width,
                   label='2007-era', color='#d62728', alpha=0.7)
    bars2 = ax.bar(x_pos + width / 2, counts_modern, width,
                   label='Modern (2024+)', color='#1f77b4', alpha=0.7)

    ax.set_xticks(x_pos)
    ax.set_xticklabels([SYSTEM_LABELS[s] for s in SYSTEM_IDS], fontsize=11)
    ax.set_ylabel('Number of rejected points (binding)', fontsize=12)
    ax.set_title('Binding constraint for rejected points', fontsize=13)
    ax.legend(fontsize=11)
    ax.tick_params(labelsize=12)

    fig.tight_layout()
    path = os.path.join(outdir, 'fig_2007_vs_modern_acceptance.pdf')
    fig.savefig(path, bbox_inches='tight', dpi=150)
    print(f'  Saved: {path}')
    plt.close(fig)


def print_summary(records, agg):
    """Print summary statistics to stdout."""
    n = len(records)
    n_mod = sum(1 for r in records if r['accepted_modern'])
    n_2007 = sum(1 for r in records if r['accepted_2007'])

    print('\n' + '=' * 70)
    print('  2007 vs Modern Comparison Summary')
    print('=' * 70)
    print(f'  Total scan points:       {n:>8,}')
    print(f'  Accepted (2007-era):     {n_2007:>8,}  ({100*n_2007/n:.2f}%)')
    print(f'  Accepted (modern):       {n_mod:>8,}  ({100*n_mod/n:.2f}%)')
    print(f'  Difference:              {n_2007 - n_mod:>+8,}  '
          f'({100*(n_2007-n_mod)/n:+.2f}%)')
    print()

    # Rescaling factors
    print('  Bound rescaling factors (modern / 2007):')
    for s in SYSTEM_IDS:
        print(f'    {s:>12s}: {RESCALE[s]:.4f}  '
              f'(modern {RESCALE[s]:.1%} of 2007)')
    print()

    # Binding system breakdown for modern-rejected but 2007-accepted
    lost = Counter()
    for rec in records:
        if rec['accepted_2007'] and not rec['accepted_modern']:
            lost[rec['binding_modern']] += 1

    print(f'  Points accepted in 2007 but rejected by modern ({sum(lost.values()):,}):')
    for s in SYSTEM_IDS:
        if lost[s] > 0:
            print(f'    {s:>12s}: {lost[s]:>6,}  '
                  f'({100*lost[s]/sum(lost.values()):.1f}%)')
    print()

    # Grid-level stats
    n_grid = len(agg['r'])
    n_grid_mod = agg['accepted_modern'].sum()
    n_grid_2007 = agg['accepted_2007'].sum()
    print(f'  Grid points (best overall_scale): {n_grid:,}')
    print(f'    Accepted 2007-era:  {n_grid_2007:>5}')
    print(f'    Accepted modern:    {n_grid_mod:>5}')
    print('=' * 70)
    print()


def main():
    parser = argparse.ArgumentParser(
        description='Compare 2007-era vs modern quark exclusion bounds',
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

    print('Generating figures ...')
    plot_exclusion_comparison(agg, args.output_dir)
    plot_mkk_bound(agg, args.output_dir)
    plot_tightening_driver(agg, args.output_dir)
    plot_acceptance_bar(records, args.output_dir)

    print('\nDone.')


if __name__ == '__main__':
    main()
