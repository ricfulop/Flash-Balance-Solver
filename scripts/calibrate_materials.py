#!/usr/bin/env python3
"""
Calibrate r_eff for Extracted Materials

Automatically adjusts r_eff (effective localization length) for each material
to match experimental Flash onset data.

Usage:
    python scripts/calibrate_materials.py
"""

import sys
from pathlib import Path
from scipy.optimize import brentq
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from flash_balance_solver import (
    MaterialParameters, MaterialFamily, FlashBalanceSolver
)


# =============================================================================
# MATERIALS TO CALIBRATE (with initial r_eff guesses)
# =============================================================================

MATERIALS_TO_CALIBRATE = {
    # =========================================================================
    # FIRST PRINCIPLES PARAMETERS - for r_eff calibration
    # =========================================================================
    
    # BiFeO3: Multiferroic perovskite
    "BiFeO3": {
        "params": MaterialParameters(
            name="Bismuth Ferrite",
            family=MaterialFamily.PEROVSKITE,
            Ea=0.70,            # Higher: mixed conduction
            sigma_0=1.0e4,      # Higher
            beta=1.50,
            alpha_res=0.35,     # Strong multiferroic coupling
            gamma=1.4,
            delta_H=-945000,    # NIST
            delta_S=-175,
            n_electrons=4,
            r_eff=40e-6,        # Initial guess (perovskite scale)
        ),
        "calibration_point": (923, 100),
    },
    
    # LLZO: Superionic Li-ion conductor
    "LLZO": {
        "params": MaterialParameters(
            name="Li7La3Zr2O12 (Garnet)",
            family=MaterialFamily.FLUORITE,
            Ea=0.30,            # Low for Li-ion
            sigma_0=1.0e5,      # High for superionic
            beta=1.50,
            alpha_res=0.20,
            gamma=1.6,
            delta_H=-1750000,
            delta_S=-200,
            n_electrons=4,
            r_eff=25e-6,
        ),
        "calibration_point": (673, 50),
    },
    
    # BaTiO3: Standard perovskite
    "BaTiO3": {
        "params": MaterialParameters(
            name="Barium Titanate",
            family=MaterialFamily.PEROVSKITE,
            Ea=0.48,
            sigma_0=8.0e2,
            beta=1.50,
            alpha_res=0.28,
            gamma=1.4,
            delta_H=-1660000,
            delta_S=-192,
            n_electrons=4,
            r_eff=40e-6,        # Perovskite scale
        ),
        "calibration_point": (923, 100),
    },
    
    # WC: Metallic carbide - use ceramic parameters but flag for special handling
    # Note: WC may not follow standard Flash Balance model due to metallic behavior
    "WC": {
        "params": MaterialParameters(
            name="Tungsten Carbide",
            family=MaterialFamily.CARBIDE,
            Ea=0.30,            # Lower than typical ceramics
            sigma_0=1.0e5,      # High conductivity
            beta=0.80,          # Lower than ceramics
            alpha_res=0.10,     # Lower coupling
            gamma=1.8,
            delta_H=-40000,     # NIST-JANAF: correct
            delta_S=-10,
            n_electrons=4,
            r_eff=10e-6,        # Try ceramic scale first
        ),
        "calibration_point": (1473, 80),
    },
    
    # MgAl2O4: Spinel
    "MgAl2O4": {
        "params": MaterialParameters(
            name="Magnesium Aluminate Spinel",
            family=MaterialFamily.SPINEL,
            Ea=1.10,
            sigma_0=1.0e5,      # Increased
            beta=1.20,
            alpha_res=0.20,
            gamma=1.8,
            delta_H=-2300000,
            delta_S=-145,
            n_electrons=4,
            r_eff=25e-6,        # Similar to Al2O3
        ),
        "calibration_point": (1323, 150),
    },
    
    # ZnO: Wurtzite n-type
    "ZnO": {
        "params": MaterialParameters(
            name="Zinc Oxide",
            family=MaterialFamily.WURTZITE,
            Ea=0.40,
            sigma_0=1.0e4,      # Higher
            beta=1.10,          # Lower for wurtzite
            alpha_res=0.22,
            gamma=1.6,
            delta_H=-350000,
            delta_S=-43,
            n_electrons=4,
            r_eff=15e-6,
        ),
        "calibration_point": (873, 80),
    },
    
    # TiO2: Standard rutile
    "TiO2": {
        "params": MaterialParameters(
            name="Titanium Dioxide",
            family=MaterialFamily.RUTILE,
            Ea=0.52,
            sigma_0=1.5e3,
            beta=1.43,
            alpha_res=0.22,
            gamma=1.8,
            delta_H=-944000,
            delta_S=-186,
            n_electrons=4,
            r_eff=18e-6,
        ),
        "calibration_point": (973, 120),
    },
    
    # 3YSZ: Standard fluorite
    "3YSZ": {
        "params": MaterialParameters(
            name="3YSZ",
            family=MaterialFamily.FLUORITE,
            Ea=0.98,
            sigma_0=2.5e4,
            beta=1.69,
            alpha_res=0.15,
            gamma=2.0,
            delta_H=-1085000,
            delta_S=-178,
            n_electrons=4,
            r_eff=14e-6,
        ),
        "calibration_point": (1173, 100),
    },
}


def calibrate_r_eff(material_params: MaterialParameters, 
                    T_target: float, E_Vcm: float,
                    r_eff_min: float = 1e-7,
                    r_eff_max: float = 1e-2) -> float:
    """
    Find r_eff that makes T_onset_predicted = T_target.
    
    Args:
        material_params: Material parameters (r_eff will be modified)
        T_target: Target onset temperature (K)
        E_Vcm: Electric field (V/cm)
        r_eff_min/max: Search bounds for r_eff
        
    Returns:
        Calibrated r_eff value (m)
    """
    E_field = E_Vcm * 100  # Convert to V/m
    
    def objective(r_eff):
        # Create new material with this r_eff
        params = MaterialParameters(
            name=material_params.name,
            family=material_params.family,
            Ea=material_params.Ea,
            sigma_0=material_params.sigma_0,
            beta=material_params.beta,
            alpha_res=material_params.alpha_res,
            gamma=material_params.gamma,
            delta_H=material_params.delta_H,
            delta_S=material_params.delta_S,
            n_electrons=material_params.n_electrons,
            r_eff=r_eff,
        )
        
        solver = FlashBalanceSolver(params)
        T_pred = solver.solve_onset_temperature(E_field)
        
        if T_pred is None:
            # If no solution, return large error
            return 1000
        
        return T_pred - T_target
    
    try:
        # Check if solution exists in range
        f_min = objective(r_eff_min)
        f_max = objective(r_eff_max)
        
        # Need opposite signs for brentq
        if f_min * f_max > 0:
            # No root in range, return best guess
            if abs(f_min) < abs(f_max):
                return r_eff_min
            else:
                return r_eff_max
        
        r_eff_calibrated = brentq(objective, r_eff_min, r_eff_max, rtol=1e-4)
        return r_eff_calibrated
        
    except Exception as e:
        print(f"  Warning: Calibration failed - {e}")
        return material_params.r_eff


def run_calibration():
    """Run calibration for all materials."""
    print("\n" + "="*80)
    print("CALIBRATING r_eff FOR EXTRACTED MATERIALS")
    print("="*80)
    print(f"\n{'Material':<15} {'T_target (K)':<12} {'E (V/cm)':<10} {'r_eff_old':<12} {'r_eff_new':<12}")
    print("-"*80)
    
    calibrated = {}
    
    for name, data in MATERIALS_TO_CALIBRATE.items():
        params = data["params"]
        T_target, E_field = data["calibration_point"]
        
        old_r_eff = params.r_eff
        new_r_eff = calibrate_r_eff(params, T_target, E_field)
        
        print(f"{name:<15} {T_target:<12.0f} {E_field:<10.0f} {old_r_eff:<12.2e} {new_r_eff:<12.2e}")
        
        # Store calibrated material
        calibrated[name] = MaterialParameters(
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
            r_eff=new_r_eff,
        )
    
    print("-"*80)
    return calibrated


def validate_calibrated(calibrated_materials: dict):
    """Validate calibrated materials against all data points."""
    print("\n" + "="*100)
    print("VALIDATION WITH CALIBRATED r_eff")
    print("="*100)
    
    # Test data points (including multiple fields per material)
    test_data = [
        ("BiFeO3", 923, 100),
        ("BiFeO3", 873, 150),
        ("LLZO", 673, 50),
        ("BaTiO3", 923, 100),
        ("BaTiO3", 873, 150),
        ("3YSZ", 1173, 100),
        ("3YSZ", 1123, 150),
        ("WC", 1473, 80),
        ("MgAl2O4", 1323, 150),
        ("ZnO", 873, 80),
        ("TiO2", 973, 120),
    ]
    
    print(f"\n{'Material':<15} {'E (V/cm)':<10} {'T_exp (K)':<10} {'T_pred (K)':<11} {'Error (%)':<10} {'Status':<8}")
    print("-"*100)
    
    total_error = 0
    count = 0
    pass_count = 0
    
    for name, T_exp, E_Vcm in test_data:
        if name not in calibrated_materials:
            continue
        
        params = calibrated_materials[name]
        solver = FlashBalanceSolver(params)
        
        E_field = E_Vcm * 100
        T_pred = solver.solve_onset_temperature(E_field)
        
        if T_pred is not None:
            error_pct = 100 * (T_pred - T_exp) / T_exp
            passed = abs(error_pct) < 15
            symbol = "✓" if passed else "✗"
            status = "PASS" if passed else "FAIL"
            
            print(f"{name:<15} {E_Vcm:<10.0f} {T_exp:<10.0f} {T_pred:<11.0f} {error_pct:<+10.1f} {symbol} {status}")
            
            total_error += abs(error_pct)
            count += 1
            if passed:
                pass_count += 1
        else:
            print(f"{name:<15} {E_Vcm:<10.0f} {T_exp:<10.0f} {'N/A':<11} {'N/A':<10} ✗ FAIL")
    
    print("-"*100)
    
    if count > 0:
        print(f"\nAverage error: {total_error/count:.1f}%")
        print(f"Passed: {pass_count}/{count} ({100*pass_count/count:.0f}%)")
    
    print("="*100)
    return calibrated_materials


def generate_python_code(calibrated_materials: dict):
    """Generate Python code for calibrated materials."""
    print("\n" + "="*80)
    print("CALIBRATED MATERIAL ENTRIES (copy to solver)")
    print("="*80 + "\n")
    
    for name, params in calibrated_materials.items():
        print(f'''    "{name}": MaterialParameters(
        name="{params.name}",
        family=MaterialFamily.{params.family.name},
        Ea={params.Ea:.3f},
        sigma_0={params.sigma_0:.2e},
        beta={params.beta:.2f},
        alpha_res={params.alpha_res:.2f},
        gamma={params.gamma:.1f},
        delta_H={int(params.delta_H)},
        delta_S={int(params.delta_S)},
        n_electrons={params.n_electrons},
        r_eff={params.r_eff:.2e},  # Calibrated
    ),
''')


def main():
    print("\n" + "#"*80)
    print("#  CALIBRATING EXTRACTED MATERIALS")
    print("#"*80)
    
    # Run calibration
    calibrated = run_calibration()
    
    # Validate
    validate_calibrated(calibrated)
    
    # Generate code
    generate_python_code(calibrated)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
