#!/usr/bin/env python3
"""
Process Full Data Papers for Flash Balance Solver

This script reads the filled-in templates for the 5 FULL_DATA papers
and generates calibrated MaterialParameters entries.

Input files:
    data/full_data_papers.csv       - Paper metadata
    data/full_data_extraction.csv   - T_onset, E_field values
    data/full_data_arrhenius.csv    - Ea, sigma_0 from digitization
    data/full_data_thermodynamics.csv - ΔH, ΔS values

Output:
    data/calibrated_materials.json  - Ready for solver
    data/calibrated_materials.py    - Python code to paste

Usage:
    python scripts/process_full_data.py
"""

import csv
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass, asdict

sys.path.insert(0, str(Path(__file__).parent.parent))

from flash_balance_solver import MaterialParameters, MaterialFamily, FlashBalanceSolver
from scipy.optimize import brentq


@dataclass
class PaperData:
    paper_id: int
    filename: str
    doi: str
    material_formula: str
    material_name: str
    material_family: str


@dataclass
class ExperimentalPoint:
    paper_id: int
    condition_num: int
    T_onset_C: float
    T_onset_K: float
    E_field_Vcm: float
    J_onset: Optional[float] = None


@dataclass
class ArrheniusData:
    paper_id: int
    Ea_eV: float
    sigma_0_S_m: float


@dataclass
class ThermoData:
    paper_id: int
    delta_H: float
    delta_S: float
    n_electrons: int


def load_papers(filepath: str) -> Dict[int, PaperData]:
    """Load paper metadata."""
    papers = {}
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            pid = int(row['paper_id'])
            papers[pid] = PaperData(
                paper_id=pid,
                filename=row['filename'],
                doi=row['doi'],
                material_formula=row['material_formula'],
                material_name=row['material_name'],
                material_family=row['material_family'],
            )
    return papers


def load_experimental(filepath: str) -> Dict[int, List[ExperimentalPoint]]:
    """Load experimental data points."""
    data = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(',')
            if len(parts) < 5 or not parts[0].isdigit():
                continue
            
            pid = int(parts[0])
            if not parts[2]:  # Skip empty T_onset
                continue
            
            point = ExperimentalPoint(
                paper_id=pid,
                condition_num=int(parts[1]),
                T_onset_C=float(parts[2]) if parts[2] else 0,
                T_onset_K=float(parts[3]) if parts[3] else (float(parts[2]) + 273.15 if parts[2] else 0),
                E_field_Vcm=float(parts[4]) if parts[4] else 0,
                J_onset=float(parts[5]) if len(parts) > 5 and parts[5] else None,
            )
            
            if pid not in data:
                data[pid] = []
            data[pid].append(point)
    
    return data


def load_arrhenius(filepath: str) -> Dict[int, ArrheniusData]:
    """Load Arrhenius parameters."""
    data = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(',')
            if len(parts) < 3 or not parts[0].isdigit():
                continue
            
            pid = int(parts[0])
            if not parts[1] or not parts[2]:  # Skip if Ea or sigma_0 missing
                continue
            
            data[pid] = ArrheniusData(
                paper_id=pid,
                Ea_eV=float(parts[1]),
                sigma_0_S_m=float(parts[2]),
            )
    
    return data


def load_thermo(filepath: str) -> Dict[int, ThermoData]:
    """Load thermodynamic data."""
    data = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(',')
            if len(parts) < 4 or not parts[0].isdigit():
                continue
            
            pid = int(parts[0])
            if not parts[1]:  # Skip if delta_H missing
                continue
            
            data[pid] = ThermoData(
                paper_id=pid,
                delta_H=float(parts[1]),
                delta_S=float(parts[2]) if parts[2] else -178,  # Default
                n_electrons=int(parts[3]) if parts[3] else 4,
            )
    
    return data


def get_family_params(family: str) -> Dict:
    """Get default family parameters."""
    families = {
        'fluorite': {'beta': 1.69, 'alpha_res': 0.15, 'gamma': 2.0},
        'perovskite': {'beta': 1.50, 'alpha_res': 0.28, 'gamma': 1.4},
        'spinel': {'beta': 1.20, 'alpha_res': 0.20, 'gamma': 1.8},
        'garnet': {'beta': 1.50, 'alpha_res': 0.20, 'gamma': 1.6},
        'rutile': {'beta': 1.43, 'alpha_res': 0.22, 'gamma': 1.8},
    }
    return families.get(family.lower(), families['fluorite'])


def get_material_family(family_str: str) -> MaterialFamily:
    """Convert string to MaterialFamily enum."""
    mapping = {
        'fluorite': MaterialFamily.FLUORITE,
        'perovskite': MaterialFamily.PEROVSKITE,
        'spinel': MaterialFamily.SPINEL,
        'garnet': MaterialFamily.FLUORITE,  # Treat as fluorite-like
        'rutile': MaterialFamily.RUTILE,
    }
    return mapping.get(family_str.lower(), MaterialFamily.FLUORITE)


def calibrate_r_eff(params: MaterialParameters, T_target: float, E_Vcm: float) -> float:
    """Calibrate r_eff to match experimental onset."""
    E_field = E_Vcm * 100  # V/m
    
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
        )
        solver = FlashBalanceSolver(test_params)
        T_pred = solver.solve_onset_temperature(E_field)
        return (T_pred - T_target) if T_pred else 1000
    
    try:
        r_eff = brentq(objective, 1e-7, 1e-1, rtol=1e-4)
        return r_eff
    except:
        return 20e-6  # Default


def process_papers():
    """Process all filled templates and generate calibrated parameters."""
    data_dir = Path(__file__).parent.parent / "data"
    
    # Load all data
    print("Loading data files...")
    papers = load_papers(data_dir / "full_data_papers.csv")
    experimental = load_experimental(data_dir / "full_data_extraction.csv")
    arrhenius = load_arrhenius(data_dir / "full_data_arrhenius.csv")
    thermo = load_thermo(data_dir / "full_data_thermodynamics.csv")
    
    print(f"  Papers: {len(papers)}")
    print(f"  With experimental data: {len(experimental)}")
    print(f"  With Arrhenius data: {len(arrhenius)}")
    print(f"  With thermo data: {len(thermo)}")
    
    # Process each paper
    results = []
    
    print("\n" + "=" * 70)
    print("PROCESSING PAPERS")
    print("=" * 70)
    
    for pid, paper in papers.items():
        print(f"\n[{pid}] {paper.material_name}")
        print(f"    DOI: {paper.doi}")
        
        # Check data availability
        has_exp = pid in experimental and len(experimental[pid]) > 0
        has_arr = pid in arrhenius
        has_thermo = pid in thermo
        
        if not has_exp:
            print("    ⚠ Missing experimental data - skipping")
            continue
        
        if not has_arr:
            print("    ⚠ Missing Arrhenius data - using family defaults")
        
        if not has_thermo:
            print("    ⚠ Missing thermodynamic data - using defaults")
        
        # Get parameters
        family_params = get_family_params(paper.material_family)
        
        if has_arr:
            Ea = arrhenius[pid].Ea_eV
            sigma_0 = arrhenius[pid].sigma_0_S_m
        else:
            # Family-based defaults
            Ea = 0.9 if paper.material_family == 'fluorite' else 0.5
            sigma_0 = 1e4
        
        if has_thermo:
            delta_H = thermo[pid].delta_H
            delta_S = thermo[pid].delta_S
            n_electrons = thermo[pid].n_electrons
        else:
            delta_H = -1000000
            delta_S = -180
            n_electrons = 4
        
        # Use first experimental point for calibration
        exp_point = experimental[pid][0]
        T_cal = exp_point.T_onset_K
        E_cal = exp_point.E_field_Vcm
        
        print(f"    Calibration point: T={T_cal:.0f}K, E={E_cal:.0f}V/cm")
        
        # Create initial parameters
        initial_params = MaterialParameters(
            name=paper.material_name,
            family=get_material_family(paper.material_family),
            Ea=Ea,
            sigma_0=sigma_0,
            beta=family_params['beta'],
            alpha_res=family_params['alpha_res'],
            gamma=family_params['gamma'],
            delta_H=delta_H,
            delta_S=delta_S,
            n_electrons=n_electrons,
            r_eff=20e-6,  # Initial guess
        )
        
        # Calibrate r_eff
        r_eff = calibrate_r_eff(initial_params, T_cal, E_cal)
        print(f"    Calibrated r_eff: {r_eff*1e6:.1f} μm")
        
        # Create final parameters
        final_params = MaterialParameters(
            name=paper.material_name,
            family=get_material_family(paper.material_family),
            Ea=Ea,
            sigma_0=sigma_0,
            beta=family_params['beta'],
            alpha_res=family_params['alpha_res'],
            gamma=family_params['gamma'],
            delta_H=delta_H,
            delta_S=delta_S,
            n_electrons=n_electrons,
            r_eff=r_eff,
        )
        
        # Validate against all experimental points
        print("    Validation:")
        solver = FlashBalanceSolver(final_params)
        for exp in experimental[pid]:
            T_pred = solver.solve_onset_temperature(exp.E_field_Vcm * 100)
            if T_pred:
                error = 100 * (T_pred - exp.T_onset_K) / exp.T_onset_K
                status = "✓" if abs(error) < 15 else "✗"
                print(f"      E={exp.E_field_Vcm}V/cm: T_exp={exp.T_onset_K:.0f}K, "
                      f"T_pred={T_pred:.0f}K, error={error:+.1f}% {status}")
        
        results.append({
            'paper_id': pid,
            'doi': paper.doi,
            'material': paper.material_formula,
            'name': paper.material_name,
            'family': paper.material_family,
            'Ea': Ea,
            'sigma_0': sigma_0,
            'beta': family_params['beta'],
            'alpha_res': family_params['alpha_res'],
            'gamma': family_params['gamma'],
            'delta_H': delta_H,
            'delta_S': delta_S,
            'n_electrons': n_electrons,
            'r_eff': r_eff,
            'calibration_T': T_cal,
            'calibration_E': E_cal,
        })
    
    # Save results
    if results:
        # JSON output
        json_path = data_dir / "calibrated_materials.json"
        with open(json_path, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\n\nJSON saved to: {json_path}")
        
        # Python code output
        py_path = data_dir / "calibrated_materials.py"
        with open(py_path, 'w') as f:
            f.write('"""Calibrated MaterialParameters from Full Data Papers."""\n\n')
            f.write('from flash_balance_solver import MaterialParameters, MaterialFamily\n\n')
            f.write('CALIBRATED_MATERIALS = {\n')
            
            for r in results:
                family_enum = f"MaterialFamily.{r['family'].upper()}"
                f.write(f'''    "{r['material']}": MaterialParameters(
        name="{r['name']}",
        family={family_enum},
        Ea={r['Ea']:.3f},
        sigma_0={r['sigma_0']:.2e},
        beta={r['beta']:.2f},
        alpha_res={r['alpha_res']:.2f},
        gamma={r['gamma']:.1f},
        delta_H={int(r['delta_H'])},
        delta_S={int(r['delta_S'])},
        n_electrons={r['n_electrons']},
        r_eff={r['r_eff']:.2e},  # Calibrated at T={r['calibration_T']:.0f}K, E={r['calibration_E']:.0f}V/cm
    ),  # DOI: {r['doi']}
''')
            
            f.write('}\n')
        print(f"Python code saved to: {py_path}")
    
    return results


def main():
    print("=" * 70)
    print("FLASH BALANCE - PROCESS FULL DATA PAPERS")
    print("=" * 70)
    
    results = process_papers()
    
    if not results:
        print("\n⚠ No complete data found!")
        print("Please fill in the template files:")
        print("  - data/full_data_extraction.csv (T_onset, E_field)")
        print("  - data/full_data_arrhenius.csv (Ea, sigma_0)")
        print("\nSee data/EXTRACTION_GUIDE.md for instructions.")
        return 1
    
    print("\n" + "=" * 70)
    print(f"Successfully processed {len(results)} materials!")
    print("=" * 70)
    
    return 0


if __name__ == "__main__":
    exit(main())
