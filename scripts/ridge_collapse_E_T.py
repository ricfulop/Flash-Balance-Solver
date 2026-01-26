#!/usr/bin/env python3
"""
Flash Ridge Collapse Analysis - E-T Onset Boundary

This script tests ChatGPT's ridge collapse using the E-T onset relationship,
which is the primary experimentally measured quantity (higher E → lower T_onset).

The transformation:
  x = E / E_R  (reduced field)
  y = T / T_R  (reduced onset temperature)
  
Plot: log₁₀(|y - 1|) vs log₁₀(|x - 1|)

For onset curves, the expected behavior is:
  - Higher E → Lower T (inverse relationship)
  - So y < 1 when x > 1 and vice versa
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
    
    if df.columns[0] == '#':
        df = df.rename(columns={'#': 'ID'})
    
    return df

def get_materials_with_E_T_variation(df, min_points=4):
    """
    Find materials where E and T both vary (not constant T experiments).
    """
    valid_materials = []
    
    for material in df['Material'].unique():
        mat_data = df[df['Material'] == material]
        mat_data = mat_data[mat_data['E(V/cm)'].notna() & mat_data['T_onset(K)'].notna()]
        mat_data = mat_data[mat_data['E(V/cm)'] > 0]
        
        if len(mat_data) < min_points:
            continue
        
        # Check that both E and T actually vary
        E_range = mat_data['E(V/cm)'].max() / mat_data['E(V/cm)'].min()
        T_range = mat_data['T_onset(K)'].max() - mat_data['T_onset(K)'].min()
        
        # Need at least 50% variation in E and 50K variation in T
        if E_range > 1.5 and T_range > 50:
            valid_materials.append(material)
    
    return valid_materials

def calculate_E_T_collapse(df, material):
    """
    Calculate ridge-normalized coordinates for E-T onset data.
    """
    mat_data = df[df['Material'] == material].copy()
    mat_data = mat_data[mat_data['E(V/cm)'].notna() & mat_data['T_onset(K)'].notna()]
    mat_data = mat_data[mat_data['E(V/cm)'] > 0]
    
    E = mat_data['E(V/cm)'].astype(float).values
    T = mat_data['T_onset(K)'].astype(float).values
    
    # Define ridge as median E, interpolated T
    E_R = np.median(E)
    
    # Interpolate T_R at E_R
    sort_idx = np.argsort(E)
    E_sorted = E[sort_idx]
    T_sorted = T[sort_idx]
    T_R = np.interp(E_R, E_sorted, T_sorted)
    
    # Reduced coordinates
    x = E / E_R
    y = T / T_R
    
    return x, y, E_R, T_R, E, T, mat_data['Family'].iloc[0]

def create_E_T_collapse_plot(df, materials, output_dir):
    """
    Create the E-T ridge collapse plot.
    """
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Color scheme
    colors = plt.cm.tab20(np.linspace(0, 1, len(materials)))
    
    # Panel 1: Raw E-T onset curves
    ax1 = axes[0]
    
    all_log_x = []
    all_log_y = []
    material_data = []
    
    for i, material in enumerate(materials):
        try:
            x, y, E_R, T_R, E, T, family = calculate_E_T_collapse(df, material)
            
            if len(x) < 4:
                continue
            
            # Plot raw E vs T
            ax1.scatter(E, T, c=[colors[i]], s=50, alpha=0.8, label=f'{material[:15]}')
            ax1.scatter([E_R], [T_R], c=[colors[i]], s=150, marker='*', edgecolors='black')
            
            # Store for collapse plot
            material_data.append({
                'material': material,
                'family': family,
                'x': x,
                'y': y,
                'E_R': E_R,
                'T_R': T_R,
                'color': colors[i]
            })
            
        except Exception as e:
            print(f"Error with {material}: {e}")
            continue
    
    ax1.set_xlabel('E (V/cm)', fontsize=12)
    ax1.set_ylabel('T_onset (K)', fontsize=12)
    ax1.set_title('Raw E-T Onset Curves with Ridge Points (*)', fontsize=14, fontweight='bold')
    ax1.legend(loc='upper right', fontsize=7, ncol=2)
    ax1.grid(True, alpha=0.3)
    ax1.invert_yaxis()  # Higher E → lower T
    
    # Panel 2: Ridge collapse
    ax2 = axes[1]
    
    for mat in material_data:
        x = mat['x']
        y = mat['y']
        color = mat['color']
        
        # Transform to log coordinates
        for xi, yi in zip(x, y):
            if abs(xi - 1) < 0.02 or abs(yi - 1) < 0.001:
                continue  # Skip points too close to ridge
            
            log_x = np.log10(np.abs(xi - 1))
            log_y = np.log10(np.abs(yi - 1))
            
            all_log_x.append(log_x)
            all_log_y.append(log_y)
            
            # Use different markers for above/below ridge
            marker = 'o' if xi > 1 else 's'
            ax2.scatter(log_x, log_y, c=[color], s=50, alpha=0.8, marker=marker)
    
    # Fit global slope
    if len(all_log_x) >= 5:
        all_log_x = np.array(all_log_x)
        all_log_y = np.array(all_log_y)
        
        slope, intercept, r, p, se = stats.linregress(all_log_x, all_log_y)
        
        x_line = np.linspace(min(all_log_x) - 0.2, max(all_log_x) + 0.2, 50)
        y_line = slope * x_line + intercept
        ax2.plot(x_line, y_line, 'k-', linewidth=2.5, 
                label=f'Empirical fit: slope = {slope:.2f} (R² = {r**2:.2f})')
    
    # Reference lines
    x_ref = np.linspace(-2, 1, 50)
    ax2.plot(x_ref, 1.0 * x_ref - 1, 'b--', alpha=0.5, linewidth=1.5, label='ν = 1 (shifted)')
    ax2.plot(x_ref, 0.5 * x_ref - 0.5, 'g--', alpha=0.5, linewidth=1.5, label='ν = 0.5 (shifted)')
    
    ax2.axhline(y=0, color='gray', linestyle='-', alpha=0.3)
    ax2.axvline(x=0, color='gray', linestyle='-', alpha=0.3)
    
    ax2.set_xlabel(r'$\log_{10}(|E/E_R - 1|)$', fontsize=12)
    ax2.set_ylabel(r'$\log_{10}(|T/T_R - 1|)$', fontsize=12)
    ax2.set_title('E-T Ridge Collapse: Universal Scaling Test', fontsize=14, fontweight='bold')
    ax2.legend(loc='lower right', fontsize=9)
    ax2.grid(True, alpha=0.3)
    
    # Add annotation
    ax2.annotate(
        '○ = E > E_R (high field)\n◼ = E < E_R (low field)',
        xy=(0.02, 0.98), xycoords='axes fraction',
        fontsize=9, ha='left', va='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    fig.suptitle('Flash Ridge Collapse: E-T Onset Boundary (Empirical)',
                fontsize=16, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    
    # Save
    output_file = output_dir / 'RIDGE_COLLAPSE_E_T_EMPIRICAL.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"✓ Saved: {output_file}")
    
    pdf_file = output_dir / 'RIDGE_COLLAPSE_E_T_EMPIRICAL.pdf'
    plt.savefig(pdf_file, bbox_inches='tight', facecolor='white')
    print(f"✓ Saved: {pdf_file}")
    
    plt.close()
    
    return slope, r**2, material_data

def main():
    script_dir = Path(__file__).parent
    data_dir = script_dir.parent / 'data'
    output_dir = data_dir / 'collapse_analysis'
    output_dir.mkdir(exist_ok=True)
    
    print("=" * 60)
    print("RIDGE COLLAPSE: E-T ONSET BOUNDARY")
    print("=" * 60)
    
    df = load_onset_data(data_dir)
    print(f"\nTotal experiments: {len(df)}")
    
    # Find materials with E-T variation
    materials = get_materials_with_E_T_variation(df, min_points=4)
    print(f"\nMaterials with E-T variation (≥4 points): {len(materials)}")
    for mat in materials:
        mat_data = df[df['Material'] == mat]
        mat_data = mat_data[mat_data['E(V/cm)'].notna() & mat_data['T_onset(K)'].notna()]
        print(f"  - {mat}: {len(mat_data)} points")
    
    print("\nCreating E-T ridge collapse plot...")
    slope, r2, material_data = create_E_T_collapse_plot(df, materials, output_dir)
    
    print("\n" + "=" * 60)
    print("RESULTS")
    print("=" * 60)
    print(f"\nGlobal slope: {slope:.2f}")
    print(f"R²: {r2:.3f}")
    
    print("\n" + "=" * 60)
    print("RIDGE PARAMETERS")
    print("=" * 60)
    print("\n{:<20} {:>10} {:>10}".format('Material', 'E_R (V/cm)', 'T_R (K)'))
    print("-" * 45)
    for mat in material_data:
        print("{:<20} {:>10.1f} {:>10.0f}".format(
            mat['material'][:19], mat['E_R'], mat['T_R']))

if __name__ == '__main__':
    main()
