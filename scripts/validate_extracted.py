#!/usr/bin/env python3
"""
Validate Extracted Material Parameters

Takes automatically extracted Arrhenius data (Ea, σ₀) from the 
auto_extract_arrhenius.py script and validates them against 
experimental flash sintering data.

Workflow:
1. Load extracted Arrhenius parameters (Ea, σ₀)
2. Combine with thermodynamic data (ΔH, ΔS) from NIST
3. Calibrate r_eff to match experimental T_onset
4. Generate validation table with accuracy metrics

Usage:
    python scripts/validate_extracted.py
"""

import sys
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import numpy as np
from scipy.optimize import brentq

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from flash_balance_solver import (
    MaterialParameters, MaterialFamily, FlashBalanceSolver
)


# =============================================================================
# EXTRACTED DATA FROM AUTO_EXTRACT_ARRHENIUS.PY
# =============================================================================

@dataclass
class ExtractedData:
    """Data extracted from papers."""
    filename: str
    material: str
    doi: str
    Ea_eV: float
    sigma_0_Sm: float
    R2: float
    confidence: str
    # Experimental flash onset data
    T_onset_C: float
    E_field_Vcm: float


# Extracted Arrhenius parameters + experimental onset data from pdf_extraction_results.csv
EXTRACTED_DATA = [
    ExtractedData(
        filename="Al2O3:TZP Steil.pdf",
        material="Al2O3:TZP",
        doi="10.1016/j.jeurceramsoc.2015.02.033",
        Ea_eV=0.782,
        sigma_0_Sm=4.36e6,
        R2=0.818,
        confidence="D",
        T_onset_C=900,  # From pdf_extraction_results.csv
        E_field_Vcm=100,  # First field value
    ),
    ExtractedData(
        filename="LLZO Huang.pdf",
        material="LLZO",
        doi="10.1007/s11581-024-06046-7",
        Ea_eV=0.304,
        sigma_0_Sm=0.765,
        R2=0.952,
        confidence="M",
        T_onset_C=700,  # From pdf_extraction_results.csv
        E_field_Vcm=90,
    ),
    ExtractedData(
        filename="Sm2O3 doped CeO2 Jiang.pdf",
        material="SDC",
        doi="10.1016/j.scriptamat.2018.09.014",
        Ea_eV=0.570,
        sigma_0_Sm=212,
        R2=0.923,
        confidence="L",
        T_onset_C=578,  # Lower onset from pdf_extraction_results.csv
        E_field_Vcm=30,  # First field value
    ),
]


# =============================================================================
# THERMODYNAMIC DATA (NIST-JANAF)
# =============================================================================

THERMO_DATA = {
    "Al2O3:TZP": {
        "delta_H": -1676000,  # J/mol (Al2O3 dominant)
        "delta_S": -210,       # J/(mol·K)
        "n_electrons": 4,
        "family": MaterialFamily.SPINEL,
        "beta": 1.20,          # Spinel-like
        "alpha_res": 0.20,
        "gamma": 1.8,
    },
    "LLZO": {
        "delta_H": -1750000,  # J/mol (garnet approximation)
        "delta_S": -200,
        "n_electrons": 4,
        "family": MaterialFamily.FLUORITE,  # Treat as fluorite-like
        "beta": 1.50,
        "alpha_res": 0.20,
        "gamma": 1.6,
    },
    "SDC": {
        "delta_H": -1024000,  # J/mol (CeO2)
        "delta_S": -195,
        "n_electrons": 4,
        "family": MaterialFamily.FLUORITE,
        "beta": 1.69,
        "alpha_res": 0.16,
        "gamma": 2.0,
    },
}


# =============================================================================
# CALIBRATION AND VALIDATION
# =============================================================================

def calibrate_r_eff(params: MaterialParameters, 
                    T_target_K: float, 
                    E_field_Vcm: float) -> float:
    """
    Calibrate r_eff to match experimental onset temperature.
    
    Args:
        params: MaterialParameters with all values except r_eff
        T_target_K: Target onset temperature (K)
        E_field_Vcm: Electric field (V/cm)
    
    Returns:
        Calibrated r_eff (m)
    """
    E_field = E_field_Vcm * 100  # Convert to V/m
    
    def objective(r_eff: float) -> float:
        test_params = MaterialParameters(
            name=params.name,
            family=params.family,
            Ea=params.Ea,
            sigma_0=params.sigma_0,
            beta=params.beta,
            alpha_res=params.alpha_res,
            gamma=params.gamma,
            delta_H=params.delta_H,
            delta_S=params.delta_S,
            n_electrons=params.n_electrons,
            r_eff=r_eff,
        )
        solver = FlashBalanceSolver(test_params)
        T_pred = solver.solve_onset_temperature(E_field)
        
        if T_pred is None:
            return 1e6  # Large error if no solution
        
        return T_pred - T_target_K
    
    # Search for r_eff in reasonable range (1 nm to 1 mm)
    try:
        r_eff = brentq(objective, 1e-9, 1e-3, rtol=1e-4)
        return r_eff
    except ValueError:
        # If brentq fails, try a grid search
        best_r_eff = 20e-6
        best_error = float('inf')
        
        for r_eff in np.logspace(-8, -3, 100):
            error = abs(objective(r_eff))
            if error < best_error:
                best_error = error
                best_r_eff = r_eff
        
        return best_r_eff


def create_material_params(data: ExtractedData) -> Tuple[MaterialParameters, float]:
    """
    Create MaterialParameters from extracted data.
    
    Returns:
        (MaterialParameters with calibrated r_eff, calibrated r_eff value)
    """
    thermo = THERMO_DATA[data.material]
    
    # Initial parameters (r_eff will be calibrated)
    initial_params = MaterialParameters(
        name=data.material,
        family=thermo["family"],
        Ea=data.Ea_eV,
        sigma_0=data.sigma_0_Sm,
        beta=thermo["beta"],
        alpha_res=thermo["alpha_res"],
        gamma=thermo["gamma"],
        delta_H=thermo["delta_H"],
        delta_S=thermo["delta_S"],
        n_electrons=thermo["n_electrons"],
        r_eff=20e-6,  # Initial guess
    )
    
    # Calibrate r_eff
    T_target_K = data.T_onset_C + 273.15
    r_eff = calibrate_r_eff(initial_params, T_target_K, data.E_field_Vcm)
    
    # Create final parameters
    final_params = MaterialParameters(
        name=data.material,
        family=thermo["family"],
        Ea=data.Ea_eV,
        sigma_0=data.sigma_0_Sm,
        beta=thermo["beta"],
        alpha_res=thermo["alpha_res"],
        gamma=thermo["gamma"],
        delta_H=thermo["delta_H"],
        delta_S=thermo["delta_S"],
        n_electrons=thermo["n_electrons"],
        r_eff=r_eff,
    )
    
    return final_params, r_eff


def validate_material(params: MaterialParameters, 
                      T_exp_K: float, 
                      E_field_Vcm: float) -> Dict:
    """
    Validate material parameters against experimental data.
    
    Returns:
        Dictionary with validation results
    """
    E_field = E_field_Vcm * 100  # V/m
    solver = FlashBalanceSolver(params)
    
    T_pred = solver.solve_onset_temperature(E_field)
    
    if T_pred is not None:
        error_K = T_pred - T_exp_K
        error_pct = 100 * error_K / T_exp_K
    else:
        error_K = None
        error_pct = None
    
    return {
        "T_exp_K": T_exp_K,
        "T_pred_K": T_pred,
        "error_K": error_K,
        "error_pct": error_pct,
        "k_soft": params.get_ksoft(),
    }


def run_validation():
    """Run full validation on extracted data."""
    results = []
    
    print("\n" + "=" * 130)
    print("FLASH BALANCE SOLVER - VALIDATION OF AUTO-EXTRACTED PARAMETERS")
    print("=" * 130)
    
    # Process each extracted material
    for data in EXTRACTED_DATA:
        params, r_eff = create_material_params(data)
        T_exp_K = data.T_onset_C + 273.15
        
        validation = validate_material(params, T_exp_K, data.E_field_Vcm)
        
        results.append({
            "material": data.material,
            "filename": data.filename,
            "doi": data.doi,
            "Ea": data.Ea_eV,
            "sigma_0": data.sigma_0_Sm,
            "beta": params.beta,
            "alpha_res": params.alpha_res,
            "r_eff": r_eff,
            "k_soft": validation["k_soft"],
            "E_field": data.E_field_Vcm,
            "T_exp_K": T_exp_K,
            "T_pred_K": validation["T_pred_K"],
            "error_pct": validation["error_pct"],
            "R2": data.R2,
            "confidence": data.confidence,
        })
    
    return results


def print_validation_table(results: List[Dict]):
    """Print formatted validation table."""
    
    print(f"\n{'Material':<12} | {'Ea(eV)':<7} | {'β':<5} | {'r_eff(μm)':<10} | {'k_soft':<6} | "
          f"{'E(V/cm)':<8} | {'T_exp(K)':<9} | {'T_pred(K)':<10} | {'Error%':<8} | DOI")
    print("-" * 130)
    
    errors = []
    for r in results:
        T_pred_str = f"{r['T_pred_K']:.0f}" if r['T_pred_K'] else "N/A"
        error_str = f"{r['error_pct']:+.1f}%" if r['error_pct'] is not None else "N/A"
        doi_short = r['doi'][:25] + "..." if len(r['doi']) > 25 else r['doi']
        
        print(f"{r['material']:<12} | {r['Ea']:<7.3f} | {r['beta']:<5.2f} | "
              f"{r['r_eff']*1e6:<10.1f} | {r['k_soft']:<6.3f} | "
              f"{r['E_field']:<8.0f} | {r['T_exp_K']:<9.0f} | {T_pred_str:<10} | "
              f"{error_str:<8} | {doi_short}")
        
        if r['error_pct'] is not None:
            errors.append(abs(r['error_pct']))
    
    print("-" * 130)
    
    # Summary statistics
    if errors:
        mae = np.mean(errors)
        within_15 = sum(1 for e in errors if e < 15) / len(errors) * 100
        max_error = max(errors)
        
        print(f"\nSUMMARY STATISTICS:")
        print(f"  Materials validated: {len(errors)}")
        print(f"  Mean Absolute Error: {mae:.1f}%")
        print(f"  Max Error: {max_error:.1f}%")
        print(f"  Within ±15% threshold: {within_15:.0f}%")
    
    print("=" * 130)


def print_detailed_table(results: List[Dict]):
    """Print detailed table with all parameters."""
    
    print("\n" + "=" * 150)
    print("DETAILED PARAMETERS TABLE")
    print("=" * 150)
    
    print(f"\n{'Material':<12} | {'Ea':<6} | {'σ₀(S/m)':<12} | {'β':<5} | {'α_res':<6} | "
          f"{'r_eff(μm)':<10} | {'ΔH(kJ/mol)':<12} | {'k_soft':<6} | {'R²':<6} | Conf")
    print("-" * 150)
    
    for r in results:
        thermo = THERMO_DATA[r['material']]
        delta_H_kJ = thermo['delta_H'] / 1000
        
        print(f"{r['material']:<12} | {r['Ea']:<6.3f} | {r['sigma_0']:<12.2e} | "
              f"{r['beta']:<5.2f} | {r['alpha_res']:<6.2f} | {r['r_eff']*1e6:<10.1f} | "
              f"{delta_H_kJ:<12.0f} | {r['k_soft']:<6.3f} | {r['R2']:<6.3f} | {r['confidence']}")
    
    print("=" * 150)


def print_comparison_table(results: List[Dict]):
    """Print comparison showing calibration vs prediction."""
    
    print("\n" + "=" * 100)
    print("CALIBRATION VALIDATION (r_eff calibrated to match T_exp)")
    print("=" * 100)
    print("\nNote: Since r_eff was calibrated to match T_exp at the given E_field,")
    print("      Error% should be ~0%. Non-zero errors indicate calibration issues.\n")
    
    print(f"{'Material':<12} | {'E_cal(V/cm)':<12} | {'T_exp(K)':<10} | {'T_pred(K)':<10} | "
          f"{'Error':<8} | {'Status':<10}")
    print("-" * 100)
    
    for r in results:
        T_pred_str = f"{r['T_pred_K']:.0f}" if r['T_pred_K'] else "N/A"
        
        if r['error_pct'] is not None:
            error_str = f"{r['error_pct']:+.2f}%"
            status = "✓ OK" if abs(r['error_pct']) < 1 else "⚠ Check"
        else:
            error_str = "N/A"
            status = "✗ Failed"
        
        print(f"{r['material']:<12} | {r['E_field']:<12.0f} | {r['T_exp_K']:<10.0f} | "
              f"{T_pred_str:<10} | {error_str:<8} | {status:<10}")
    
    print("=" * 100)


def save_results_csv(results: List[Dict], output_path: str):
    """Save results to CSV."""
    import csv
    
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'material', 'doi', 'Ea_eV', 'sigma_0_Sm', 'beta', 'alpha_res', 
            'r_eff_m', 'k_soft', 'E_field_Vcm', 'T_exp_K', 'T_pred_K', 
            'error_pct', 'R2_arrhenius', 'confidence'
        ])
        
        for r in results:
            writer.writerow([
                r['material'], r['doi'], r['Ea'], r['sigma_0'], r['beta'],
                r['alpha_res'], r['r_eff'], r['k_soft'], r['E_field'],
                r['T_exp_K'], r['T_pred_K'], r['error_pct'], r['R2'], r['confidence']
            ])
    
    print(f"\nResults saved to: {output_path}")


def main():
    """Main entry point."""
    print("\n" + "=" * 130)
    print("FLASH BALANCE SOLVER - EXTRACTED DATA VALIDATION")
    print("=" * 130)
    print("\nInput: Automatically extracted Arrhenius parameters from PDFs")
    print("Process: Calibrate r_eff, predict T_onset, calculate accuracy")
    
    # Run validation
    results = run_validation()
    
    # Print tables
    print_validation_table(results)
    print_detailed_table(results)
    print_comparison_table(results)
    
    # Save to CSV
    output_path = Path(__file__).parent.parent / "data" / "extracted_validation.csv"
    save_results_csv(results, str(output_path))
    
    return results


if __name__ == "__main__":
    main()
