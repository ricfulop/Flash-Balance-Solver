#!/usr/bin/env python3
"""
Flash Ridge Collapse Analysis - Empirical Test

This script tests ChatGPT's ridge collapse methodology using REAL onset data
from published papers to see if the universal exponents ν₋ ≈ 1 and ν₊ ≈ 2
actually emerge from empirical measurements.

The transformation:
  x = E / E_R  (reduced field)
  y = J / J_R  (reduced current density)
  
Plot: log₁₀(|y - 1|) vs log₁₀(|x - 1|)

Expected slopes if universality holds:
  - Below ridge (x < 1): ν₋ ≈ 1 (Ohmic/defect-limited)
  - Above ridge (x > 1): ν₊ ≈ 2 (Runaway/Flash)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

def load_onset_data(data_dir):
    """Load onset data from validation table."""
    onset_file = data_dir / 'complete_validation_table.csv'
    df = pd.read_csv(onset_file)
    df.columns = df.columns.str.strip()
    
    # Rename first column
    if df.columns[0] == '#':
        df = df.rename(columns={'#': 'ID'})
    
    return df

def get_materials_with_multiple_points(df, min_points=3):
    """Find materials with enough data points for collapse analysis."""
    
    # Filter for valid E and J data
    df_valid = df[df['E(V/cm)'].notna() & (df['E(V/cm)'] > 0)]
    
    # Use J_calc if available
    if 'J_calc(mA/mm²)' in df_valid.columns:
        df_valid = df_valid[df_valid['J_calc(mA/mm²)'].notna() & (df_valid['J_calc(mA/mm²)'] > 0)]
    
    # Count points per material
    counts = df_valid.groupby('Material').size()
    valid_materials = counts[counts >= min_points].index.tolist()
    
    return valid_materials, df_valid

def calculate_ridge_collapse(df, material):
    """
    Calculate ridge-normalized coordinates for a single material.
    
    E_R = median field (or could be defined by physics)
    J_R = current density at E_R
    
    Returns x = E/E_R, y = J/J_R for all data points.
    """
    mat_data = df[df['Material'] == material].copy()
    
    E = mat_data['E(V/cm)'].astype(float).values
    J = mat_data['J_calc(mA/mm²)'].astype(float).values
    T = mat_data['T_onset(K)'].astype(float).values
    
    # Define ridge as median field point
    E_R = np.median(E)
    
    # Find J_R by interpolation at E_R
    # Sort by E for interpolation
    sort_idx = np.argsort(E)
    E_sorted = E[sort_idx]
    J_sorted = J[sort_idx]
    
    # Linear interpolation in log space for J_R
    if len(E) > 1:
        J_R = np.interp(E_R, E_sorted, J_sorted)
    else:
        J_R = J[0]
    
    # Calculate reduced coordinates
    x = E / E_R
    y = J / J_R
    
    return x, y, E_R, J_R, mat_data['Family'].iloc[0]

def create_ridge_collapse_plot(df, materials, output_dir):
    """
    Create the ridge collapse plot following ChatGPT's methodology.
    
    Plot log₁₀(y-1) vs log₁₀(|x-1|) with separate fits for x<1 and x>1.
    """
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Color scheme by family
    family_colors = {
        'fluorite': '#1f77b4',
        'wurtzite': '#ff7f0e',
        'perovskite': '#2ca02c',
        'spinel': '#d62728',
        'corundum': '#9467bd',
        'carbide': '#8c564b',
        'glass': '#e377c2',
        'oxide': '#7f7f7f',
        'rutile': '#bcbd22',
        'garnet': '#17becf',
    }
    
    # Collect all data for global fit
    all_x_below = []
    all_y_below = []
    all_x_above = []
    all_y_above = []
    
    plotted_materials = []
    material_stats = []
    
    # Panel 1: Raw E-J data with ridge markers
    ax1 = axes[0]
    
    for material in materials:
        try:
            x, y, E_R, J_R, family = calculate_ridge_collapse(df, material)
            
            if len(x) < 3:
                continue
            
            color = family_colors.get(family, '#333333')
            
            # Plot raw E vs J
            mat_data = df[df['Material'] == material]
            E = mat_data['E(V/cm)'].astype(float).values
            J = mat_data['J_calc(mA/mm²)'].astype(float).values
            
            ax1.scatter(E, J, c=color, s=40, alpha=0.7, label=f'{material}')
            ax1.axvline(E_R, color=color, linestyle='--', alpha=0.3)
            ax1.scatter([E_R], [J_R], c=color, s=100, marker='*', edgecolors='black')
            
            plotted_materials.append(material)
            material_stats.append({
                'Material': material,
                'Family': family,
                'E_R': E_R,
                'J_R': J_R,
                'n_points': len(x)
            })
            
        except Exception as e:
            print(f"Error with {material}: {e}")
            continue
    
    ax1.set_xlabel('E (V/cm)', fontsize=12)
    ax1.set_ylabel('J (mA/mm²)', fontsize=12)
    ax1.set_title('Raw Onset Data with Ridge Points (*)', fontsize=14, fontweight='bold')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.legend(loc='upper left', fontsize=7, ncol=2)
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Ridge collapse plot
    ax2 = axes[1]
    
    for material in plotted_materials:
        try:
            x, y, E_R, J_R, family = calculate_ridge_collapse(df, material)
            
            color = family_colors.get(family, '#333333')
            
            # Calculate collapse coordinates
            # Exclude points exactly at ridge (x=1)
            mask_valid = np.abs(x - 1) > 0.01
            x_valid = x[mask_valid]
            y_valid = y[mask_valid]
            
            if len(x_valid) < 2:
                continue
            
            # Separate below and above ridge
            mask_below = x_valid < 1
            mask_above = x_valid > 1
            
            # Transform to log coordinates
            for x_pt, y_pt in zip(x_valid, y_valid):
                log_x = np.log10(np.abs(x_pt - 1))
                
                # For y, we need y > 1 for log(y-1) to work
                # If y < 1, we use |y-1| and note it's below ridge response
                if y_pt > 1:
                    log_y = np.log10(y_pt - 1)
                else:
                    log_y = np.log10(np.abs(y_pt - 1))
                
                if x_pt < 1:
                    all_x_below.append(log_x)
                    all_y_below.append(log_y)
                else:
                    all_x_above.append(log_x)
                    all_y_above.append(log_y)
                
                marker = 'o' if x_pt > 1 else 's'
                ax2.scatter(log_x, log_y, c=color, s=40, alpha=0.7, marker=marker)
            
        except Exception as e:
            continue
    
    # Fit slopes for both branches
    results = {}
    
    if len(all_x_below) >= 3:
        all_x_below = np.array(all_x_below)
        all_y_below = np.array(all_y_below)
        slope_below, intercept_below, r_below, p_below, se_below = stats.linregress(all_x_below, all_y_below)
        
        x_line = np.linspace(min(all_x_below), max(all_x_below), 50)
        y_line = slope_below * x_line + intercept_below
        ax2.plot(x_line, y_line, 'b-', linewidth=2.5, 
                label=f'Below ridge: ν₋ = {slope_below:.2f} (R² = {r_below**2:.2f})')
        
        results['nu_minus'] = slope_below
        results['r2_minus'] = r_below**2
    
    if len(all_x_above) >= 3:
        all_x_above = np.array(all_x_above)
        all_y_above = np.array(all_y_above)
        slope_above, intercept_above, r_above, p_above, se_above = stats.linregress(all_x_above, all_y_above)
        
        x_line = np.linspace(min(all_x_above), max(all_x_above), 50)
        y_line = slope_above * x_line + intercept_above
        ax2.plot(x_line, y_line, 'r-', linewidth=2.5,
                label=f'Above ridge: ν₊ = {slope_above:.2f} (R² = {r_above**2:.2f})')
        
        results['nu_plus'] = slope_above
        results['r2_plus'] = r_above**2
    
    # Add reference lines for expected slopes
    x_ref = np.linspace(-2, 0.5, 50)
    ax2.plot(x_ref, 1.0 * x_ref, 'k--', alpha=0.4, linewidth=1.5, label='Reference: ν = 1')
    ax2.plot(x_ref, 2.0 * x_ref, 'k:', alpha=0.4, linewidth=1.5, label='Reference: ν = 2')
    
    ax2.axvline(x=0, color='gray', linestyle='-', alpha=0.3)
    ax2.axhline(y=0, color='gray', linestyle='-', alpha=0.3)
    
    ax2.set_xlabel(r'$\log_{10}(|E/E_R - 1|)$', fontsize=12)
    ax2.set_ylabel(r'$\log_{10}(|J/J_R - 1|)$', fontsize=12)
    ax2.set_title('Ridge Collapse: Empirical Test of Universal Exponents', fontsize=14, fontweight='bold')
    ax2.legend(loc='lower right', fontsize=9)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(-2.5, 1.0)
    ax2.set_ylim(-3.0, 2.0)
    
    # Add annotation
    annotation = f"""
    ◼ = Below ridge (E < E_R)
    ○ = Above ridge (E > E_R)
    
    ChatGPT prediction:
      ν₋ ≈ 1 (Ohmic)
      ν₊ ≈ 2 (Runaway)
    """
    ax2.annotate(annotation, xy=(0.02, 0.98), xycoords='axes fraction',
                fontsize=8, ha='left', va='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    fig.suptitle('Flash Ridge Universality Test: Real Onset Data from 5+ Papers',
                fontsize=16, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    
    # Save
    output_file = output_dir / 'RIDGE_COLLAPSE_EMPIRICAL.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"✓ Saved: {output_file}")
    
    pdf_file = output_dir / 'RIDGE_COLLAPSE_EMPIRICAL.pdf'
    plt.savefig(pdf_file, bbox_inches='tight', facecolor='white')
    print(f"✓ Saved: {pdf_file}")
    
    plt.close()
    
    return results, material_stats

def main():
    # Set up paths
    script_dir = Path(__file__).parent
    data_dir = script_dir.parent / 'data'
    output_dir = data_dir / 'collapse_analysis'
    output_dir.mkdir(exist_ok=True)
    
    print("=" * 60)
    print("RIDGE COLLAPSE ANALYSIS - EMPIRICAL TEST")
    print("Testing ChatGPT's Universal Exponents with Real Data")
    print("=" * 60)
    
    # Load data
    print("\nLoading onset data...")
    df = load_onset_data(data_dir)
    print(f"  Total experiments: {len(df)}")
    
    # Find materials with multiple points
    valid_materials, df_valid = get_materials_with_multiple_points(df, min_points=3)
    print(f"\nMaterials with ≥3 data points: {len(valid_materials)}")
    for mat in valid_materials[:10]:
        count = len(df_valid[df_valid['Material'] == mat])
        print(f"  - {mat}: {count} points")
    if len(valid_materials) > 10:
        print(f"  ... and {len(valid_materials) - 10} more")
    
    # Create ridge collapse plot
    print("\nCreating ridge collapse plot...")
    results, material_stats = create_ridge_collapse_plot(df_valid, valid_materials, output_dir)
    
    # Print summary
    print("\n" + "=" * 60)
    print("RESULTS: EMPIRICAL EXPONENTS")
    print("=" * 60)
    
    if 'nu_minus' in results:
        print(f"\nBelow-ridge exponent (ν₋):")
        print(f"  Measured: {results['nu_minus']:.2f}")
        print(f"  Expected: 1.0 (Ohmic)")
        print(f"  R²: {results['r2_minus']:.3f}")
        
    if 'nu_plus' in results:
        print(f"\nAbove-ridge exponent (ν₊):")
        print(f"  Measured: {results['nu_plus']:.2f}")
        print(f"  Expected: 2.0 (Runaway)")
        print(f"  R²: {results['r2_plus']:.3f}")
    
    print("\n" + "=" * 60)
    print("INTERPRETATION")
    print("=" * 60)
    
    if 'nu_minus' in results and 'nu_plus' in results:
        nu_m = results['nu_minus']
        nu_p = results['nu_plus']
        
        if 0.5 < nu_m < 1.5 and 1.5 < nu_p < 2.5:
            print("\n✓ UNIVERSAL EXPONENTS SUPPORTED!")
            print("  Empirical data collapses with ν₋ ≈ 1 and ν₊ ≈ 2")
            print("  This supports flash as a universal critical phenomenon.")
        elif abs(nu_m - nu_p) < 0.5:
            print("\n⚠ Single slope regime")
            print("  Both branches show similar exponents.")
            print("  The ridge may not be well-defined in this data.")
        else:
            print("\n⚠ Mixed results")
            print(f"  ν₋ = {nu_m:.2f} (expected ~1)")
            print(f"  ν₊ = {nu_p:.2f} (expected ~2)")
            print("  Universality partially supported.")
    
    print("\n" + "=" * 60)
    print("MATERIALS USED IN ANALYSIS")
    print("=" * 60)
    
    print("\n{:<25} {:<12} {:>10} {:>12} {:>6}".format(
        'Material', 'Family', 'E_R (V/cm)', 'J_R (mA/mm²)', 'N'))
    print("-" * 65)
    for stat in material_stats:
        print("{:<25} {:<12} {:>10.1f} {:>12.2f} {:>6}".format(
            stat['Material'][:24], stat['Family'][:11], 
            stat['E_R'], stat['J_R'], stat['n_points']))

if __name__ == '__main__':
    main()
