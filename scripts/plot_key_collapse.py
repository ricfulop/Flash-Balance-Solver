#!/usr/bin/env python3
"""
Create the key collapse figure for publication.
This is the "money plot" showing Flash as a universal state of matter.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

kB = 8.617e-5  # Boltzmann constant (eV/K)

FAMILY_COLORS = {
    'fluorite': '#1f77b4', 'perovskite': '#ff7f0e', 'spinel': '#2ca02c',
    'rutile': '#d62728', 'wurtzite': '#9467bd', 'corundum': '#8c564b',
    'carbide': '#e377c2', 'nitride': '#7f7f7f', 'garnet': '#bcbd22',
    'oxide': '#17becf', 'metal': '#000000', 'glass': '#aec7e8',
    'glass_ceramic': '#98df8a', 'ferrite': '#c49c94', 'composite': '#ffbb78',
    'high_entropy_oxide': '#ff9896', 'pyrochlore': '#c5b0d5',
    'actinide-oxide': '#f7b6d2', 'oxide_composite': '#c7c7c7',
}

def main():
    # Load data
    data_dir = Path(__file__).parent.parent / 'data'
    df = pd.read_csv(data_dir / 'complete_validation_table.csv')
    
    # Calculate conductivity
    df['sigma_onset'] = df['sigma_0(S/m)'] * np.exp(-df['Ea(eV)'] / (kB * df['T_onset(K)']))
    df = df.dropna(subset=['sigma_onset', 'E(V/cm)'])
    
    # Universal parameter
    df['E_sigma_half'] = df['E(V/cm)'] * np.sqrt(df['sigma_onset'])
    
    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # ========== PLOT 1: THE COLLAPSE ==========
    ax = axes[0]
    
    families = sorted(df['Family'].unique())
    for family in families:
        mask = df['Family'] == family
        color = FAMILY_COLORS.get(family, '#333333')
        ax.scatter(df.loc[mask, 'T_onset(K)'], df.loc[mask, 'E_sigma_half'],
                   c=color, label=family, alpha=0.7, s=80, edgecolors='white', linewidth=0.5)
    
    ax.set_yscale('log')
    
    # Universal band
    log_vals = np.log10(df['E_sigma_half'].dropna())
    mean_val = 10**log_vals.mean()
    std_decades = log_vals.std()
    
    ax.axhline(mean_val, color='red', linestyle='-', linewidth=3, alpha=0.8)
    ax.axhspan(10**(log_vals.mean() - std_decades), 10**(log_vals.mean() + std_decades),
               alpha=0.15, color='red')
    
    ax.set_xlabel('Flash Onset Temperature (K)', fontsize=14)
    ax.set_ylabel('E × σ$^{1/2}$ (V·S$^{1/2}$·cm$^{-1}$·m$^{-1/2}$)', fontsize=14)
    ax.set_title(f'Universal Flash Condition\n'
                 f'E × √σ = {mean_val:.0f} ± {std_decades:.2f} decades', 
                 fontsize=16, fontweight='bold')
    ax.legend(loc='upper left', fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(200, 2500)
    
    # ========== PLOT 2: HISTOGRAM ==========
    ax = axes[1]
    
    ax.hist(log_vals, bins=25, edgecolor='black', alpha=0.7, color='steelblue', density=True)
    
    # Fit normal distribution
    x_fit = np.linspace(log_vals.min(), log_vals.max(), 100)
    y_fit = stats.norm.pdf(x_fit, log_vals.mean(), log_vals.std())
    ax.plot(x_fit, y_fit, 'r-', linewidth=3, label=f'Normal fit\nσ = {std_decades:.2f} decades')
    
    ax.axvline(log_vals.mean(), color='red', linestyle='--', linewidth=2)
    
    ax.set_xlabel('log₁₀(E × σ$^{1/2}$)', fontsize=14)
    ax.set_ylabel('Probability Density', fontsize=14)
    ax.set_title('Distribution of Universal Parameter\n'
                 f'n = {len(log_vals)} materials, 19 families', fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    
    # ========== PLOT 3: E vs σ^(-1/2) ==========
    ax = axes[2]
    
    sigma_inv_sqrt = 1 / np.sqrt(df['sigma_onset'])
    
    for family in families:
        mask = df['Family'] == family
        color = FAMILY_COLORS.get(family, '#333333')
        ax.scatter(sigma_inv_sqrt[mask], df.loc[mask, 'E(V/cm)'],
                   c=color, label=family, alpha=0.7, s=60)
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    # Power law fit
    x = np.log10(sigma_inv_sqrt)
    y = np.log10(df['E(V/cm)'])
    valid = np.isfinite(x) & np.isfinite(y)
    slope, intercept, r_value, _, _ = stats.linregress(x[valid], y[valid])
    
    x_fit = np.logspace(np.log10(sigma_inv_sqrt.min()), np.log10(sigma_inv_sqrt.max()), 100)
    y_fit = 10**(slope * np.log10(x_fit) + intercept)
    ax.plot(x_fit, y_fit, 'k--', linewidth=3, 
            label=f'Fit: E ∝ σ$^{{{-slope/2:.2f}}}$\nR² = {r_value**2:.3f}')
    
    # Theoretical line (slope = 1, i.e., E ∝ σ^(-0.5))
    y_theory = mean_val / np.sqrt(1/x_fit**2)  # E = const × σ^(-0.5)
    ax.plot(x_fit, mean_val * x_fit, 'r-', linewidth=2, alpha=0.5,
            label='Theory: E ∝ σ$^{-0.5}$')
    
    ax.set_xlabel('σ$^{-1/2}$ [(S/m)$^{-1/2}$]', fontsize=14)
    ax.set_ylabel('E-field (V/cm)', fontsize=14)
    ax.set_title('Critical Field Scaling\n'
                 f'Theory: E ∝ σ$^{{-0.5}}$, Measured: σ$^{{{-slope/2:.2f}}}$', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')
    
    plt.tight_layout()
    
    # Save
    output_dir = data_dir / 'collapse_analysis'
    output_dir.mkdir(exist_ok=True)
    plt.savefig(output_dir / 'KEY_COLLAPSE_FIGURE.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'KEY_COLLAPSE_FIGURE.pdf', bbox_inches='tight')
    print(f"Saved to {output_dir / 'KEY_COLLAPSE_FIGURE.png'}")
    
    # Print summary
    print("\n" + "="*70)
    print("KEY RESULT: UNIVERSAL FLASH CONDITION")
    print("="*70)
    print(f"\n  E × √σ = {mean_val:.1f} at flash onset")
    print(f"  Spread = {std_decades:.2f} decades (log₁₀)")
    print(f"  Data: {len(df)} experiments, 19 material families")
    print(f"\n  E scales as σ^{-slope/2:.2f} (theory predicts -0.50)")
    print(f"  Scaling R² = {r_value**2:.3f}")
    print("\n" + "="*70)
    print("INTERPRETATION:")
    print("="*70)
    if std_decades < 1.0:
        print("\n  ✓ STRONG COLLAPSE: Data from 19 material families collapse")
        print("    onto a single universal band spanning less than 1 decade!")
        print("\n  This is compelling evidence that Flash represents a")
        print("  UNIVERSAL THERMODYNAMIC STATE accessible through")
        print("  electric field - phonon coupling.")
    elif std_decades < 1.5:
        print("\n  ○ MODERATE COLLAPSE: Good universality observed")
    else:
        print("\n  × Weak collapse - family-specific effects dominate")
    print("="*70)


if __name__ == '__main__':
    main()
