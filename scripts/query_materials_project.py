#!/usr/bin/env python3
"""
Query Materials Project for First-Principles Properties
=========================================================

This script queries the Materials Project database for Born effective charges,
dielectric constants, and other phonon-related properties that could be used
to derive β (ridge parameter) from first principles.

Goal: Test if β can be predicted from independently measured material properties,
making the k_soft collapse analysis non-circular.

Requirements:
    pip install mp-api

Usage:
    1. Get an API key from https://materialsproject.org/api
    2. Set environment variable: export MP_API_KEY="your_key_here"
    3. Run: python query_materials_project.py

Author: Flash Balance Solver Project
Date: January 2026
"""

import os
import json
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, Dict, List
import pandas as pd
import numpy as np

# Material mapping: validation table materials -> base compounds for MP query
# Many flash materials are doped versions - we query the base compound
MATERIAL_TO_MP_FORMULA = {
    # Fluorites
    '8YSZ': 'ZrO2',
    '3YSZ': 'ZrO2', 
    'GDC10': 'CeO2',
    'SDC': 'CeO2',
    '(ZrO₂)₀.₈(CeO₂)₀.₂': 'ZrO2',
    'Irradiated-YSZ': 'ZrO2',
    
    # Rutiles
    'TiO2': 'TiO2',
    'TiO₂': 'TiO2',
    'SnO2': 'SnO2',
    'SnO₂-AC': 'SnO2',
    
    # Perovskites
    'BaTiO3': 'BaTiO3',
    'SrTiO3': 'SrTiO3',
    'SrTiO₃': 'SrTiO3',
    'NaNbO3': 'NaNbO3',
    'KNN': 'KNbO3',  # K0.5Na0.5NbO3 approximated
    'PZT': 'PbTiO3',  # Lead zirconate titanate approximated
    
    # Spinels/Corundum
    'Al2O3': 'Al2O3',
    'Al₂O₃': 'Al2O3',
    'MgAl2O4': 'MgAl2O4',
    'MnCo2O4': 'MnCo2O4',
    
    # Other oxides
    'ZnO': 'ZnO',
    'MgO': 'MgO',
    'SiO2': 'SiO2',
    
    # Carbides/Nitrides
    'SiC': 'SiC',
    'TiN': 'TiN',
    'Si3N4': 'Si3N4',
    
    # Metals
    'W': 'W',
    'Ni': 'Ni',
    'Re': 'Re',
}

# Known β values from the material database (for correlation testing)
KNOWN_BETA_VALUES = {
    'ZrO2': 1.69,   # Fluorites
    'CeO2': 1.69,
    'TiO2': 1.43,   # Rutiles
    'SnO2': 1.40,
    'BaTiO3': 1.50, # Perovskites
    'SrTiO3': 1.50,
    'PbTiO3': 1.50,
    'Al2O3': 1.20,  # Corundum
    'MgAl2O4': 1.20, # Spinel
    'SiC': 0.23,    # Carbide
    'TiN': 0.40,    # Nitride
    'W': 0.50,      # Metal
    'Ni': 0.50,
}

# Known Debye temperatures (K) from CRC Handbook
KNOWN_T_DEBYE = {
    'ZrO2': 660,
    'CeO2': 370,
    'TiO2': 670,
    'SnO2': 650,
    'BaTiO3': 600,
    'SrTiO3': 420,
    'Al2O3': 1047,
    'MgAl2O4': 850,
    'MgO': 946,
    'ZnO': 416,
    'SiC': 1200,
    'TiN': 580,
    'Si3N4': 920,
    'W': 343,
    'Ni': 450,
}

# Nominal ionic charges
NOMINAL_CHARGES = {
    'Zr': 4, 'Ce': 4, 'Ti': 4, 'Sn': 4, 'Ba': 2, 'Sr': 2,
    'Pb': 2, 'Al': 3, 'Mg': 2, 'Zn': 2, 'Si': 4, 'N': -3,
    'O': -2, 'C': -4, 'W': 0, 'Ni': 0, 'Re': 0, 'K': 1, 'Na': 1, 'Nb': 5,
}


@dataclass
class MPMaterialData:
    """Data retrieved from Materials Project."""
    formula: str
    mp_id: str
    # Dielectric properties
    epsilon_electronic: Optional[float] = None  # ε∞
    epsilon_ionic: Optional[float] = None       # ε_ionic
    epsilon_total: Optional[float] = None       # ε_0 = ε∞ + ε_ionic
    # Born effective charges (average magnitude)
    born_charge_cation: Optional[float] = None
    born_charge_anion: Optional[float] = None
    born_charge_anomaly: Optional[float] = None  # |Z*|/|Z_nominal| - 1
    # Other properties
    band_gap: Optional[float] = None
    formation_energy: Optional[float] = None
    # Phonon properties (if available)
    # Note: PhononDB would have more detailed phonon info


def query_materials_project(api_key: str, formulas: List[str]) -> Dict[str, MPMaterialData]:
    """
    Query Materials Project for dielectric and Born charge data.
    
    Args:
        api_key: Materials Project API key
        formulas: List of chemical formulas to query
        
    Returns:
        Dictionary mapping formula to MPMaterialData
    """
    try:
        from mp_api.client import MPRester
    except ImportError:
        print("ERROR: mp-api not installed. Run: pip install mp-api")
        return {}
    
    results = {}
    
    with MPRester(api_key) as mpr:
        for formula in formulas:
            print(f"Querying {formula}...")
            try:
                # Search for the material
                docs = mpr.materials.summary.search(
                    formula=formula,
                    fields=["material_id", "formula_pretty", "band_gap", 
                           "formation_energy_per_atom"]
                )
                
                if not docs:
                    print(f"  No results for {formula}")
                    continue
                
                # Get the most stable phase (lowest formation energy)
                doc = min(docs, key=lambda x: x.formation_energy_per_atom or 0)
                mp_id = doc.material_id
                
                data = MPMaterialData(
                    formula=formula,
                    mp_id=str(mp_id),
                    band_gap=doc.band_gap,
                    formation_energy=doc.formation_energy_per_atom,
                )
                
                # Try to get dielectric data
                try:
                    dielectric_docs = mpr.materials.dielectric.search(
                        material_ids=[mp_id]
                    )
                    if dielectric_docs:
                        diel = dielectric_docs[0]
                        # Dielectric tensor is 3x3, take average of diagonal
                        if hasattr(diel, 'e_electronic') and diel.e_electronic is not None:
                            data.epsilon_electronic = np.mean(np.diag(diel.e_electronic))
                        if hasattr(diel, 'e_ionic') and diel.e_ionic is not None:
                            data.epsilon_ionic = np.mean(np.diag(diel.e_ionic))
                        if hasattr(diel, 'e_total') and diel.e_total is not None:
                            data.epsilon_total = np.mean(np.diag(diel.e_total))
                        
                        # Born effective charges
                        if hasattr(diel, 'born') and diel.born is not None:
                            # Born is a list of 3x3 tensors for each atom
                            born_charges = []
                            for bc in diel.born:
                                # Average of diagonal elements
                                avg_charge = np.mean(np.diag(bc))
                                born_charges.append(abs(avg_charge))
                            
                            if born_charges:
                                # Separate cations (positive) and anions (negative)
                                # For simplicity, take max as cation, min as anion proxy
                                data.born_charge_cation = max(born_charges)
                                data.born_charge_anion = min(born_charges)
                                
                except Exception as e:
                    print(f"  No dielectric data for {formula}: {e}")
                
                results[formula] = data
                print(f"  Found: {mp_id}, Eg={data.band_gap:.2f} eV" if data.band_gap else f"  Found: {mp_id}")
                
            except Exception as e:
                print(f"  Error querying {formula}: {e}")
    
    return results


def calculate_born_anomaly(data: MPMaterialData, formula: str) -> Optional[float]:
    """
    Calculate Born effective charge anomaly: |Z*|/|Z_nominal| - 1
    
    Positive values indicate enhanced dynamic charge transfer (covalent character).
    """
    if data.born_charge_cation is None:
        return None
    
    # Get nominal cation charge for this compound
    # Simple heuristic: first element is usually the cation
    cation = formula.replace('O2', '').replace('O3', '').replace('O4', '').replace('O', '')
    cation = ''.join([c for c in cation if c.isalpha()])[:2]  # First 2 letters
    
    nominal = NOMINAL_CHARGES.get(cation, 4)  # Default to +4 for oxides
    if nominal == 0:  # Metals
        return None
    
    anomaly = abs(data.born_charge_cation) / abs(nominal) - 1
    return anomaly


def test_beta_correlation(mp_data: Dict[str, MPMaterialData]) -> None:
    """
    Test if β correlates with first-principles properties.
    """
    print("\n" + "="*70)
    print("TESTING β CORRELATION WITH FIRST-PRINCIPLES PROPERTIES")
    print("="*70)
    
    # Collect data for correlation
    rows = []
    for formula, data in mp_data.items():
        if formula not in KNOWN_BETA_VALUES:
            continue
        
        beta = KNOWN_BETA_VALUES[formula]
        T_D = KNOWN_T_DEBYE.get(formula)
        
        # Calculate Born anomaly
        born_anomaly = calculate_born_anomaly(data, formula)
        
        rows.append({
            'formula': formula,
            'beta': beta,
            'T_debye': T_D,
            'epsilon_electronic': data.epsilon_electronic,
            'epsilon_total': data.epsilon_total,
            'born_cation': data.born_charge_cation,
            'born_anomaly': born_anomaly,
            'band_gap': data.band_gap,
        })
    
    df = pd.DataFrame(rows)
    
    if len(df) < 3:
        print("Not enough data for correlation analysis")
        return
    
    print(f"\nData collected for {len(df)} materials:")
    print(df.to_string(index=False))
    
    # Test correlations
    print("\n" + "-"*70)
    print("CORRELATION ANALYSIS")
    print("-"*70)
    
    from scipy import stats
    
    # β vs various properties
    correlations = []
    
    for col in ['T_debye', 'epsilon_electronic', 'epsilon_total', 'born_cation', 'born_anomaly', 'band_gap']:
        valid = df[['beta', col]].dropna()
        if len(valid) >= 3:
            r, p = stats.pearsonr(valid['beta'], valid[col])
            correlations.append({
                'property': col,
                'r': r,
                'r²': r**2,
                'p_value': p,
                'n': len(valid),
            })
    
    corr_df = pd.DataFrame(correlations)
    if len(corr_df) > 0:
        print("\nCorrelation of β with material properties:")
        print(corr_df.to_string(index=False))
        
        # Find best predictor
        best = corr_df.loc[corr_df['r²'].idxmax()]
        print(f"\nBest predictor: {best['property']} (R² = {best['r²']:.3f}, p = {best['p_value']:.4f})")
    
    # Test combined formula: β ∝ f(T_D, ε, Z*)
    print("\n" + "-"*70)
    print("TESTING SEMI-EMPIRICAL FORMULA")
    print("-"*70)
    
    # Try: β ∝ (1 + born_anomaly)² / (T_D / T_D_ref)
    T_D_ref = 600  # Reference Debye temperature
    
    valid = df[['beta', 'T_debye', 'born_anomaly']].dropna()
    if len(valid) >= 3:
        # Proposed formula
        valid['beta_pred'] = (1 + valid['born_anomaly'].clip(lower=0))**2 / (valid['T_debye'] / T_D_ref)
        
        # Scale to match
        scale = valid['beta'].mean() / valid['beta_pred'].mean()
        valid['beta_pred_scaled'] = valid['beta_pred'] * scale
        
        r, p = stats.pearsonr(valid['beta'], valid['beta_pred_scaled'])
        print(f"\nFormula: β ∝ (1 + Z*_anomaly)² / (T_D / {T_D_ref}K)")
        print(f"Correlation: R² = {r**2:.3f}, p = {p:.4f}")
        print(f"Scale factor: {scale:.2f}")
        
        print("\nPredicted vs Actual β:")
        for _, row in valid.iterrows():
            print(f"  {row.name}: actual={row['beta']:.2f}, predicted={row['beta_pred_scaled']:.2f}")


def save_results(mp_data: Dict[str, MPMaterialData], output_path: Path) -> None:
    """Save queried data to JSON for future use."""
    data_dict = {}
    for formula, data in mp_data.items():
        data_dict[formula] = {
            'mp_id': data.mp_id,
            'epsilon_electronic': data.epsilon_electronic,
            'epsilon_ionic': data.epsilon_ionic,
            'epsilon_total': data.epsilon_total,
            'born_charge_cation': data.born_charge_cation,
            'born_charge_anion': data.born_charge_anion,
            'band_gap': data.band_gap,
            'formation_energy': data.formation_energy,
        }
    
    with open(output_path, 'w') as f:
        json.dump(data_dict, f, indent=2)
    print(f"\n✓ Saved results to {output_path}")


def main():
    """Main function."""
    print("="*70)
    print("MATERIALS PROJECT QUERY FOR FIRST-PRINCIPLES β DERIVATION")
    print("="*70)
    
    # Check for API key
    api_key = os.environ.get('MP_API_KEY')
    
    # Always run literature analysis for comparison
    print("\nRunning analysis with LITERATURE VALUES (Born charges from DFT papers)...")
    
    # Literature values from DFT studies (PRB, JACS, etc.)
    # Born charges: Z* for cation from various DFT papers
    literature_data = {
            'ZrO2': MPMaterialData(
                formula='ZrO2', mp_id='mp-2858',
                epsilon_electronic=4.9, epsilon_total=22.0,
                born_charge_cation=5.8, born_charge_anion=2.9,
                band_gap=3.9
            ),
            'TiO2': MPMaterialData(
                formula='TiO2', mp_id='mp-2657',
                epsilon_electronic=7.4, epsilon_total=89.0,
                born_charge_cation=7.3, born_charge_anion=2.0,
                band_gap=2.0
            ),
            'BaTiO3': MPMaterialData(
                formula='BaTiO3', mp_id='mp-5020',
                epsilon_electronic=6.8, epsilon_total=130.0,
                born_charge_cation=7.5, born_charge_anion=2.1,
                band_gap=2.0
            ),
            'SrTiO3': MPMaterialData(
                formula='SrTiO3', mp_id='mp-5229',
                epsilon_electronic=6.0, epsilon_total=310.0,
                born_charge_cation=7.1, born_charge_anion=2.1,
                band_gap=2.2
            ),
            'Al2O3': MPMaterialData(
                formula='Al2O3', mp_id='mp-1143',
                epsilon_electronic=3.1, epsilon_total=10.5,
                born_charge_cation=2.9, born_charge_anion=1.9,
                band_gap=6.0
            ),
            'CeO2': MPMaterialData(
                formula='CeO2', mp_id='mp-20194',
                epsilon_electronic=5.3, epsilon_total=23.0,
                born_charge_cation=5.4, born_charge_anion=2.7,
                band_gap=2.5
            ),
            'MgO': MPMaterialData(
                formula='MgO', mp_id='mp-1265',
                epsilon_electronic=3.0, epsilon_total=9.8,
                born_charge_cation=1.98, born_charge_anion=1.98,
                band_gap=4.8
            ),
        }
        
    test_beta_correlation(literature_data)
    
    # Save literature data
    output_dir = Path(__file__).parent.parent / 'data'
    save_results(literature_data, output_dir / 'materials_project_literature.json')
    
    # If API key available, also query Materials Project
    if api_key:
        formulas = list(set(KNOWN_BETA_VALUES.keys()))
        print(f"\n\nQuerying {len(formulas)} materials from Materials Project API...")
        mp_data = query_materials_project(api_key, formulas)
        
        if mp_data:
            print("\n" + "="*70)
            print("MATERIALS PROJECT API RESULTS")
            print("="*70)
            test_beta_correlation(mp_data)
            save_results(mp_data, output_dir / 'materials_project_api.json')


if __name__ == '__main__':
    main()
