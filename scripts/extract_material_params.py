#!/usr/bin/env python3
"""
Flash Balance First Principles Parameter Extraction Tool

This interactive script guides users through extracting Flash Balance parameters
for new materials using first principles data from thermodynamic databases
and electrical property measurements.

The workflow:
1. Material identification and classification
2. Thermodynamic data lookup (ΔH, ΔS)
3. Electrical properties extraction (Ea, σ₀)
4. Family-based parameter estimation (β, α_res, γ)
5. r_eff calibration using experimental onset data
6. Validation and output

Usage:
    python scripts/extract_material_params.py
    python scripts/extract_material_params.py --non-interactive --config material.json

For supplementary material documentation, see:
    docs/FIRST_PRINCIPLES_WORKFLOW.md

Author: Flash Balance Solver Team
"""

import sys
import json
import argparse
from pathlib import Path
from dataclasses import dataclass, asdict
from typing import Optional, Dict, Tuple
from scipy.optimize import brentq
import numpy as np

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from flash_balance_solver import (
    MaterialParameters, MaterialFamily, FlashBalanceSolver
)


# =============================================================================
# PARAMETER GUIDELINES BY MATERIAL FAMILY
# =============================================================================

FAMILY_GUIDELINES = {
    MaterialFamily.FLUORITE: {
        "description": "Fluorite oxides (e.g., YSZ, GDC, CeO2)",
        "beta": {"typical": 1.69, "range": (1.5, 1.9), "uncertainty": 0.10},
        "alpha_res": {"typical": 0.15, "range": (0.12, 0.18), "uncertainty": 0.03},
        "gamma": {"typical": 2.0, "range": (1.8, 2.2), "uncertainty": 0.2},
        "ksoft": {"typical": 0.10, "note": "Calculated from β"},
        "Ea_range": (0.7, 1.1),  # eV
        "sigma_0_range": (1e4, 1e5),  # S/m
        "r_eff_range": (10e-6, 30e-6),  # m
        "n_electrons": 4,  # Per O2
        "examples": ["8YSZ", "3YSZ", "GDC10", "CeO2"],
    },
    MaterialFamily.RUTILE: {
        "description": "Rutile oxides (e.g., TiO2, SnO2)",
        "beta": {"typical": 1.43, "range": (1.3, 1.5), "uncertainty": 0.08},
        "alpha_res": {"typical": 0.22, "range": (0.18, 0.26), "uncertainty": 0.05},
        "gamma": {"typical": 1.8, "range": (1.6, 2.0), "uncertainty": 0.2},
        "ksoft": {"typical": 0.24, "note": "Calculated from β"},
        "Ea_range": (0.4, 0.6),  # eV
        "sigma_0_range": (1e3, 5e3),  # S/m
        "r_eff_range": (12e-6, 25e-6),  # m
        "n_electrons": 4,
        "examples": ["TiO2", "SnO2"],
    },
    MaterialFamily.PEROVSKITE: {
        "description": "Perovskite oxides (e.g., SrTiO3, BaTiO3, BiFeO3)",
        "beta": {"typical": 1.50, "range": (1.3, 1.7), "uncertainty": 0.10},
        "alpha_res": {"typical": 0.28, "range": (0.20, 0.35), "uncertainty": 0.05},
        "gamma": {"typical": 1.4, "range": (1.2, 1.6), "uncertainty": 0.2},
        "ksoft": {"typical": 0.20, "note": "Calculated from β"},
        "Ea_range": (0.3, 0.8),  # eV
        "sigma_0_range": (5e2, 1e4),  # S/m
        "r_eff_range": (15e-6, 55e-6),  # m (larger due to polar domains)
        "n_electrons": 4,
        "examples": ["SrTiO3", "BaTiO3", "BiFeO3", "PZT"],
        "note": "Multiferroics may have higher α_res (0.30-0.40)",
    },
    MaterialFamily.SPINEL: {
        "description": "Spinel/Corundum oxides (e.g., Al2O3, MgAl2O4)",
        "beta": {"typical": 1.20, "range": (1.0, 1.4), "uncertainty": 0.15},
        "alpha_res": {"typical": 0.20, "range": (0.15, 0.25), "uncertainty": 0.04},
        "gamma": {"typical": 1.8, "range": (1.6, 2.0), "uncertainty": 0.2},
        "ksoft": {"typical": 0.36, "note": "Calculated from β"},
        "Ea_range": (1.0, 1.4),  # eV (high activation)
        "sigma_0_range": (5e4, 2e5),  # S/m
        "r_eff_range": (20e-6, 120e-6),  # m
        "n_electrons": 4,
        "examples": ["Al2O3", "MgAl2O4", "MgO"],
    },
    MaterialFamily.WURTZITE: {
        "description": "Wurtzite oxides (e.g., ZnO)",
        "beta": {"typical": 1.10, "range": (0.9, 1.3), "uncertainty": 0.10},
        "alpha_res": {"typical": 0.22, "range": (0.18, 0.26), "uncertainty": 0.04},
        "gamma": {"typical": 1.6, "range": (1.4, 1.8), "uncertainty": 0.2},
        "ksoft": {"typical": 0.41, "note": "Calculated from β"},
        "Ea_range": (0.3, 0.5),  # eV
        "sigma_0_range": (5e3, 2e4),  # S/m
        "r_eff_range": (10e-6, 25e-6),  # m
        "n_electrons": 4,
        "examples": ["ZnO"],
    },
    MaterialFamily.CARBIDE: {
        "description": "Carbides (e.g., SiC, WC, TiC)",
        "beta": {"typical": 1.20, "range": (0.8, 1.4), "uncertainty": 0.20},
        "alpha_res": {"typical": 0.08, "range": (0.05, 0.12), "uncertainty": 0.03},
        "gamma": {"typical": 1.8, "range": (1.6, 2.0), "uncertainty": 0.2},
        "ksoft": {"typical": 0.88, "note": "High for covalent bonds"},
        "Ea_range": (0.3, 1.0),  # eV (varies widely)
        "sigma_0_range": (1e4, 1e6),  # S/m
        "r_eff_range": (4e-6, 15e-6),  # m
        "n_electrons": 4,
        "examples": ["SiC", "WC", "TiC"],
        "note": "WC is metallic - may not follow standard model",
    },
    MaterialFamily.NITRIDE: {
        "description": "Nitrides (e.g., Si3N4, TiN, AlN)",
        "beta": {"typical": 0.40, "range": (0.3, 0.6), "uncertainty": 0.20},
        "alpha_res": {"typical": 0.12, "range": (0.08, 0.16), "uncertainty": 0.03},
        "gamma": {"typical": 1.6, "range": (1.4, 1.8), "uncertainty": 0.2},
        "ksoft": {"typical": 0.79, "note": "Specified override"},
        "Ea_range": (0.05, 0.9),  # eV
        "sigma_0_range": (1e3, 5e6),  # S/m
        "r_eff_range": (3e-3, 20e-6),  # m (TiN is metallic)
        "n_electrons": 6,  # Per N2
        "examples": ["Si3N4", "TiN", "AlN"],
    },
    MaterialFamily.METAL: {
        "description": "Metals (electroplasticity regime)",
        "beta": {"typical": 0.05, "range": (0.02, 0.08), "uncertainty": 0.05},
        "alpha_res": {"typical": 0.05, "range": (0.02, 0.08), "uncertainty": 0.02},
        "gamma": {"typical": 1.2, "range": (1.0, 1.4), "uncertainty": 0.2},
        "ksoft": {"typical": 0.98, "note": "Nearly unity for metals"},
        "Ea_range": (0.0, 0.05),  # eV
        "sigma_0_range": (1e6, 6e7),  # S/m
        "r_eff_range": (1e-3, 10e-3),  # m (gauge length scale)
        "n_electrons": 2,  # Per metal atom
        "examples": ["Cu", "Ni", "W", "Al"],
        "note": "Flash in metals is Joule heating dominated",
    },
    MaterialFamily.HYDRIDE: {
        "description": "Hydrides (e.g., ZrH2, TiH2)",
        "beta": {"typical": 0.55, "range": (0.4, 0.7), "uncertainty": 0.15},
        "alpha_res": {"typical": 0.30, "range": (0.25, 0.35), "uncertainty": 0.05},
        "gamma": {"typical": 1.5, "range": (1.3, 1.7), "uncertainty": 0.2},
        "ksoft": {"typical": 0.71, "note": "Calculated from β"},
        "Ea_range": (0.2, 0.4),  # eV
        "sigma_0_range": (1e4, 1e5),  # S/m
        "r_eff_range": (8e-6, 15e-6),  # m
        "n_electrons": 2,  # Per H2
        "examples": ["ZrH2", "TiH2"],
    },
}


# =============================================================================
# THERMODYNAMIC DATA SOURCES
# =============================================================================

THERMO_SOURCES = {
    "NIST-JANAF": {
        "url": "https://janaf.nist.gov/",
        "description": "NIST-JANAF Thermochemical Tables",
        "reliability": "High",
        "coverage": "Common oxides, carbides, nitrides",
    },
    "HSC Chemistry": {
        "url": "https://www.hsc-chemistry.com/",
        "description": "HSC Chemistry Database",
        "reliability": "High",
        "coverage": "Comprehensive inorganic compounds",
    },
    "FactSage": {
        "url": "https://www.factsage.com/",
        "description": "FactSage Thermodynamic Database",
        "reliability": "High",
        "coverage": "Slags, solutions, complex systems",
    },
    "SGTE": {
        "url": "https://www.sgte.net/",
        "description": "SGTE Pure Substance Database",
        "reliability": "High",
        "coverage": "Pure substances, alloys",
    },
    "Literature": {
        "url": "N/A",
        "description": "Published experimental data",
        "reliability": "Variable",
        "coverage": "Specific to paper",
    },
    "Estimated": {
        "url": "N/A",
        "description": "Estimated from similar compounds",
        "reliability": "Low",
        "coverage": "When no data available",
    },
}


# =============================================================================
# COMMON FORMATION REACTIONS
# =============================================================================

FORMATION_REACTIONS = {
    "oxide_from_elements": "{metal}(s) + O2(g) → {formula}(s)",
    "perovskite_from_oxides": "{A}O(s) + {B}2O3(s) → 2{A}{B}O3(s)",
    "spinel_from_oxides": "{A}O(s) + {B}2O3(s) → {A}{B}2O4(s)",
    "carbide_from_elements": "{metal}(s) + C(s) → {formula}(s)",
    "nitride_from_elements": "{metal}(s) + N2(g) → {formula}(s)",
}


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def calculate_ksoft(beta: float) -> float:
    """Calculate softening factor from beta parameter."""
    q_ratio = 0.73
    p = 2
    return max(0.01, 1 - beta * (q_ratio ** p))


def estimate_delta_S(delta_H: float, family: MaterialFamily) -> float:
    """
    Estimate formation entropy from formation enthalpy.
    
    Uses typical ΔS/ΔH ratios for different material families.
    This is a rough estimate - use actual data when available.
    """
    # Typical ratios in J/(mol·K) per kJ/mol
    ratios = {
        MaterialFamily.FLUORITE: 0.16,
        MaterialFamily.RUTILE: 0.20,
        MaterialFamily.PEROVSKITE: 0.12,
        MaterialFamily.SPINEL: 0.09,
        MaterialFamily.WURTZITE: 0.12,
        MaterialFamily.CARBIDE: 0.25,
        MaterialFamily.NITRIDE: 0.45,
        MaterialFamily.METAL: 0.60,
        MaterialFamily.HYDRIDE: 0.80,
    }
    ratio = ratios.get(family, 0.15)
    # delta_H in J/mol, convert to kJ for ratio
    return -abs(delta_H / 1000) * ratio


def calibrate_r_eff(params: MaterialParameters, 
                    T_target: float, 
                    E_Vcm: float,
                    r_eff_min: float = 1e-7,
                    r_eff_max: float = 1e-2) -> Tuple[float, bool]:
    """
    Calibrate r_eff to match experimental onset temperature.
    
    Returns:
        (calibrated_r_eff, hit_bound): The calibrated value and whether it hit a bound
    """
    E_field = E_Vcm * 100  # Convert to V/m
    
    def objective(r_eff):
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
            ksoft=params.ksoft,
        )
        solver = FlashBalanceSolver(test_params)
        T_pred = solver.solve_onset_temperature(E_field)
        
        if T_pred is None:
            return 1000
        return T_pred - T_target
    
    try:
        f_min = objective(r_eff_min)
        f_max = objective(r_eff_max)
        
        if f_min * f_max > 0:
            # No root in range
            if abs(f_min) < abs(f_max):
                return r_eff_min, True
            else:
                return r_eff_max, True
        
        r_eff_cal = brentq(objective, r_eff_min, r_eff_max, rtol=1e-4)
        return r_eff_cal, False
        
    except Exception:
        return params.r_eff, True


def validate_material(params: MaterialParameters, 
                      T_exp: float, 
                      E_Vcm: float) -> Dict:
    """Validate material parameters against experimental data."""
    solver = FlashBalanceSolver(params)
    E_field = E_Vcm * 100
    T_pred = solver.solve_onset_temperature(E_field)
    
    if T_pred is not None:
        error = T_pred - T_exp
        error_pct = 100 * error / T_exp
        passed = abs(error_pct) < 15
    else:
        error = None
        error_pct = None
        passed = False
    
    return {
        "T_experimental": T_exp,
        "T_predicted": T_pred,
        "error_K": error,
        "error_percent": error_pct,
        "passed": passed,
    }


# =============================================================================
# INTERACTIVE EXTRACTION WORKFLOW
# =============================================================================

def print_header():
    """Print the tool header."""
    print("\n" + "=" * 70)
    print("  FLASH BALANCE FIRST PRINCIPLES PARAMETER EXTRACTION TOOL")
    print("=" * 70)
    print("\nThis tool guides you through extracting Flash Balance parameters")
    print("for new materials using first principles thermodynamic and")
    print("electrical property data.\n")


def print_family_options():
    """Print available material families."""
    print("\nAvailable material families:")
    print("-" * 50)
    for i, (family, info) in enumerate(FAMILY_GUIDELINES.items(), 1):
        examples = ", ".join(info["examples"][:3])
        print(f"  {i:2}. {family.value:<15} - {info['description'][:30]}")
        print(f"      Examples: {examples}")
    print()


def get_user_input(prompt: str, default=None, input_type=str, 
                   valid_range=None) -> any:
    """Get validated user input."""
    while True:
        if default is not None:
            full_prompt = f"{prompt} [{default}]: "
        else:
            full_prompt = f"{prompt}: "
        
        response = input(full_prompt).strip()
        
        if response == "" and default is not None:
            return default
        
        try:
            value = input_type(response)
            
            if valid_range is not None:
                if value < valid_range[0] or value > valid_range[1]:
                    print(f"  Warning: Value outside typical range {valid_range}")
                    confirm = input("  Use anyway? (y/n): ").lower()
                    if confirm != 'y':
                        continue
            
            return value
            
        except ValueError:
            print(f"  Invalid input. Please enter a {input_type.__name__}.")


def select_family() -> MaterialFamily:
    """Interactive family selection."""
    print_family_options()
    
    families = list(FAMILY_GUIDELINES.keys())
    while True:
        try:
            choice = int(input("Select family number: "))
            if 1 <= choice <= len(families):
                return families[choice - 1]
            print(f"  Please enter a number between 1 and {len(families)}")
        except ValueError:
            print("  Please enter a valid number")


def extract_parameters_interactive() -> Tuple[MaterialParameters, Dict]:
    """
    Interactive parameter extraction workflow.
    
    Returns:
        (MaterialParameters, metadata_dict)
    """
    print_header()
    
    # =========================================================================
    # STEP 1: Material Identification
    # =========================================================================
    print("\n" + "=" * 50)
    print("STEP 1: Material Identification")
    print("=" * 50)
    
    name = get_user_input("Material name (e.g., 'Bismuth Ferrite')")
    formula = get_user_input("Chemical formula (e.g., 'BiFeO3')")
    family = select_family()
    
    guidelines = FAMILY_GUIDELINES[family]
    print(f"\n  Selected: {family.value}")
    print(f"  Description: {guidelines['description']}")
    if "note" in guidelines:
        print(f"  Note: {guidelines['note']}")
    
    # =========================================================================
    # STEP 2: Thermodynamic Data
    # =========================================================================
    print("\n" + "=" * 50)
    print("STEP 2: Thermodynamic Data")
    print("=" * 50)
    
    print("\nRecommended data sources:")
    for source, info in list(THERMO_SOURCES.items())[:4]:
        print(f"  - {source}: {info['description']}")
    
    print(f"\nEnter formation enthalpy for: {formula}")
    print("  (Negative for exothermic formation)")
    
    delta_H_kJ = get_user_input("ΔH_f (kJ/mol)", input_type=float)
    delta_H = delta_H_kJ * 1000  # Convert to J/mol
    
    # Estimate or input ΔS
    estimated_delta_S = estimate_delta_S(delta_H, family)
    print(f"\n  Estimated ΔS from correlation: {estimated_delta_S:.0f} J/(mol·K)")
    
    use_estimate = input("  Use estimate? (y/n) [y]: ").lower()
    if use_estimate == 'n':
        delta_S = get_user_input("ΔS_f (J/(mol·K))", input_type=float)
    else:
        delta_S = estimated_delta_S
    
    thermo_source = get_user_input(
        "Data source (NIST-JANAF/HSC/FactSage/Literature/Estimated)",
        default="Literature"
    )
    
    # =========================================================================
    # STEP 3: Electrical Properties
    # =========================================================================
    print("\n" + "=" * 50)
    print("STEP 3: Electrical Properties (from Arrhenius plot)")
    print("=" * 50)
    
    print(f"\n  Typical Ea range for {family.value}: "
          f"{guidelines['Ea_range'][0]:.2f} - {guidelines['Ea_range'][1]:.2f} eV")
    print(f"  Typical σ₀ range: {guidelines['sigma_0_range'][0]:.0e} - "
          f"{guidelines['sigma_0_range'][1]:.0e} S/m")
    
    Ea = get_user_input(
        "Activation energy Ea (eV)",
        default=np.mean(guidelines['Ea_range']),
        input_type=float,
        valid_range=guidelines['Ea_range']
    )
    
    sigma_0 = get_user_input(
        "Pre-exponential σ₀ (S/m)",
        default=np.sqrt(guidelines['sigma_0_range'][0] * guidelines['sigma_0_range'][1]),
        input_type=float
    )
    
    print("\n  Conduction types: ionic, electronic, mixed")
    cond_type = get_user_input("Conduction type", default="ionic")
    
    # =========================================================================
    # STEP 4: Family-Based Parameter Estimates
    # =========================================================================
    print("\n" + "=" * 50)
    print("STEP 4: Family-Based Parameter Estimates")
    print("=" * 50)
    
    beta_info = guidelines["beta"]
    alpha_info = guidelines["alpha_res"]
    gamma_info = guidelines["gamma"]
    
    print(f"\n  Based on {family.value} family:")
    print(f"  β (ridge parameter): {beta_info['typical']:.2f} ± {beta_info['uncertainty']:.2f}")
    print(f"  α_res (coupling):    {alpha_info['typical']:.2f} ± {alpha_info['uncertainty']:.2f}")
    print(f"  γ (damping):         {gamma_info['typical']:.1f} ± {gamma_info['uncertainty']:.1f}")
    
    use_defaults = input("\n  Use family defaults? (y/n) [y]: ").lower()
    
    if use_defaults == 'n':
        beta = get_user_input("β", default=beta_info['typical'], input_type=float)
        alpha_res = get_user_input("α_res", default=alpha_info['typical'], input_type=float)
        gamma = get_user_input("γ", default=gamma_info['typical'], input_type=float)
    else:
        beta = beta_info['typical']
        alpha_res = alpha_info['typical']
        gamma = gamma_info['typical']
    
    # Calculate ksoft
    ksoft = calculate_ksoft(beta)
    print(f"\n  Calculated k_soft = {ksoft:.3f}")
    
    # Check for ksoft override
    if "ksoft" in guidelines and guidelines["ksoft"].get("note") == "Specified override":
        print(f"  Note: {family.value} typically uses override k_soft = {guidelines['ksoft']['typical']}")
        use_override = input("  Use override? (y/n) [n]: ").lower()
        if use_override == 'y':
            ksoft = guidelines['ksoft']['typical']
    
    # n_electrons
    n_electrons = guidelines.get("n_electrons", 4)
    print(f"\n  Using n_electrons = {n_electrons} (per reaction unit)")
    
    # =========================================================================
    # STEP 5: Calibration Point
    # =========================================================================
    print("\n" + "=" * 50)
    print("STEP 5: Experimental Calibration Point")
    print("=" * 50)
    
    print("\n  Enter ONE experimental Flash onset data point for calibration.")
    print("  This will be used to optimize r_eff.")
    
    T_onset = get_user_input("T_onset (K)", input_type=float)
    E_field = get_user_input("E_field (V/cm)", input_type=float)
    
    # Initial r_eff estimate
    r_eff_range = guidelines['r_eff_range']
    r_eff_init = np.sqrt(r_eff_range[0] * r_eff_range[1])
    
    # =========================================================================
    # STEP 6: Optimize r_eff
    # =========================================================================
    print("\n" + "=" * 50)
    print("STEP 6: Calibrating r_eff")
    print("=" * 50)
    
    # Create initial parameters
    initial_params = MaterialParameters(
        name=name,
        family=family,
        Ea=Ea,
        sigma_0=sigma_0,
        beta=beta,
        alpha_res=alpha_res,
        gamma=gamma,
        delta_H=delta_H,
        delta_S=delta_S,
        n_electrons=n_electrons,
        r_eff=r_eff_init,
        ksoft=ksoft if family in [MaterialFamily.METAL, MaterialFamily.NITRIDE, 
                                   MaterialFamily.CARBIDE] else None,
    )
    
    print(f"\n  Calibrating r_eff to match T_onset = {T_onset} K at E = {E_field} V/cm...")
    
    r_eff_cal, hit_bound = calibrate_r_eff(initial_params, T_onset, E_field)
    
    if hit_bound:
        print(f"\n  ⚠️  WARNING: r_eff hit calibration bound ({r_eff_cal:.2e} m)")
        print("     This may indicate other parameters need adjustment.")
        print(f"     Typical range for {family.value}: "
              f"{r_eff_range[0]*1e6:.1f} - {r_eff_range[1]*1e6:.1f} μm")
    else:
        print(f"\n  ✓ Calibrated r_eff = {r_eff_cal:.2e} m ({r_eff_cal*1e6:.1f} μm)")
    
    # Create final parameters
    final_params = MaterialParameters(
        name=name,
        family=family,
        Ea=Ea,
        sigma_0=sigma_0,
        beta=beta,
        alpha_res=alpha_res,
        gamma=gamma,
        delta_H=delta_H,
        delta_S=delta_S,
        n_electrons=n_electrons,
        r_eff=r_eff_cal,
        ksoft=ksoft if family in [MaterialFamily.METAL, MaterialFamily.NITRIDE,
                                   MaterialFamily.CARBIDE] else None,
    )
    
    # Validate
    validation = validate_material(final_params, T_onset, E_field)
    
    print(f"\n  Validation at calibration point:")
    print(f"    T_exp = {T_onset:.0f} K, T_pred = {validation['T_predicted']:.0f} K")
    print(f"    Error = {validation['error_percent']:.1f}%")
    
    # =========================================================================
    # OUTPUT
    # =========================================================================
    print("\n" + "=" * 50)
    print("OUTPUT: MaterialParameters Entry")
    print("=" * 50)
    
    ksoft_str = f"ksoft={ksoft:.2f}," if final_params.ksoft else ""
    
    print(f'''
    "{formula}": MaterialParameters(
        name="{name}",
        family=MaterialFamily.{family.name},
        Ea={Ea:.3f},           # eV
        sigma_0={sigma_0:.2e},  # S/m
        beta={beta:.2f},
        alpha_res={alpha_res:.2f},
        gamma={gamma:.1f},
        delta_H={int(delta_H)},   # J/mol
        delta_S={int(delta_S)},   # J/(mol·K)
        n_electrons={n_electrons},
        r_eff={r_eff_cal:.2e},    # m (calibrated)
        {ksoft_str}
    ),
''')
    
    # Create metadata
    metadata = {
        "formula": formula,
        "thermo_source": thermo_source,
        "conduction_type": cond_type,
        "calibration_point": {
            "T_onset_K": T_onset,
            "E_field_Vcm": E_field,
        },
        "r_eff_hit_bound": hit_bound,
        "validation": validation,
    }
    
    return final_params, metadata


def save_output(params: MaterialParameters, metadata: Dict, output_dir: Path):
    """Save extracted parameters to files."""
    output_dir.mkdir(exist_ok=True)
    
    # Save as JSON
    output_data = {
        "parameters": {
            "name": params.name,
            "family": params.family.value,
            "Ea": params.Ea,
            "sigma_0": params.sigma_0,
            "beta": params.beta,
            "alpha_res": params.alpha_res,
            "gamma": params.gamma,
            "delta_H": params.delta_H,
            "delta_S": params.delta_S,
            "n_electrons": params.n_electrons,
            "r_eff": params.r_eff,
            "ksoft": params.ksoft,
        },
        "metadata": metadata,
    }
    
    formula = metadata["formula"]
    json_path = output_dir / f"{formula}_params.json"
    with open(json_path, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\n  Saved to: {json_path}")
    return json_path


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Flash Balance First Principles Parameter Extraction Tool"
    )
    parser.add_argument(
        "--output-dir", "-o",
        default="data/extracted_params",
        help="Output directory for saved parameters"
    )
    parser.add_argument(
        "--non-interactive",
        action="store_true",
        help="Run in non-interactive mode with config file"
    )
    parser.add_argument(
        "--config",
        help="JSON config file for non-interactive mode"
    )
    
    args = parser.parse_args()
    
    if args.non_interactive:
        if not args.config:
            print("Error: --config required for non-interactive mode")
            return 1
        print("Non-interactive mode not yet implemented")
        return 1
    
    # Interactive mode
    try:
        params, metadata = extract_parameters_interactive()
        
        # Offer to save
        save = input("\nSave parameters to file? (y/n) [y]: ").lower()
        if save != 'n':
            output_dir = Path(args.output_dir)
            save_output(params, metadata, output_dir)
        
        print("\n" + "=" * 70)
        print("  Parameter extraction complete!")
        print("=" * 70)
        print("\nNext steps:")
        print("  1. Add the MaterialParameters entry to your solver")
        print("  2. Test with additional experimental data points")
        print("  3. Refine parameters if needed")
        print("\nFor documentation, see: docs/FIRST_PRINCIPLES_WORKFLOW.md")
        
        return 0
        
    except KeyboardInterrupt:
        print("\n\nExtraction cancelled.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
