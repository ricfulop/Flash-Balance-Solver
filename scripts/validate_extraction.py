#!/usr/bin/env python3
"""
Cross-Validation and Outlier Detection for Extracted Data

Compares extracted values across papers for the same material,
flags outliers and inconsistencies, and generates a quality report.

Usage:
    python validate_extraction.py <extraction_csv> [options]

Checks:
    1. Cross-paper consistency: Same material should have similar values
    2. Family consistency: Values should match typical ranges for material family
    3. Arrhenius consistency: Ea and sigma_0 should be physically reasonable
    4. Onset consistency: T_onset should decrease with increasing E_field
"""

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import defaultdict
import numpy as np


# =============================================================================
# FAMILY EXPECTATIONS
# =============================================================================

FAMILY_EXPECTATIONS = {
    'fluorite': {
        'Ea_eV': (0.7, 1.2),
        'T_onset_typical': (1000, 1300),  # K at ~100 V/cm
        'description': 'Ionic conductors (YSZ, GDC, CeO2)'
    },
    'rutile': {
        'Ea_eV': (0.3, 0.7),
        'T_onset_typical': (800, 1100),
        'description': 'N-type semiconductors (TiO2, SnO2)'
    },
    'perovskite': {
        'Ea_eV': (0.3, 0.7),
        'T_onset_typical': (700, 1000),
        'description': 'Mixed conductors (SrTiO3, BaTiO3, BiFeO3)'
    },
    'spinel': {
        'Ea_eV': (1.0, 1.8),
        'T_onset_typical': (1200, 1600),
        'description': 'Ionic insulators (Al2O3, MgAl2O4)'
    },
    'carbide': {
        'Ea_eV': (0.5, 1.5),
        'T_onset_typical': (1300, 1800),
        'description': 'Covalent semiconductors (SiC, WC, B4C)'
    },
    'nitride': {
        'Ea_eV': (0.05, 1.2),
        'T_onset_typical': (1200, 1600),
        'description': 'Metallic to semiconductor (TiN, ZrN, Si3N4)'
    },
    'garnet': {
        'Ea_eV': (0.2, 0.5),
        'T_onset_typical': (600, 900),
        'description': 'Fast ion conductors (LLZO)'
    },
    'wurtzite': {
        'Ea_eV': (0.2, 0.5),
        'T_onset_typical': (700, 1000),
        'description': 'Wide bandgap semiconductors (ZnO)'
    },
    'metal': {
        'Ea_eV': (0.0, 0.1),
        'T_onset_typical': (300, 800),
        'description': 'Metallic conductors (W, Ni, Re)'
    },
    'composite_oxide': {
        'Ea_eV': (0.3, 1.5),
        'T_onset_typical': (800, 1400),
        'description': 'Multi-phase composites'
    }
}


# =============================================================================
# DATA LOADING
# =============================================================================

def load_extraction_data(filepath: str) -> List[Dict]:
    """Load extraction CSV into list of dictionaries."""
    data = []
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Convert numeric fields
            for field in ['T_onset_K', 'T_onset_C', 'E_field_Vcm', 'J_crit_Acm2',
                          'Ea_eV', 'sigma_0_Sm', 'delta_H_Jmol', 'delta_S_JmolK',
                          'grain_size_um', 'completeness_score']:
                if row.get(field):
                    try:
                        row[field] = float(row[field])
                    except ValueError:
                        row[field] = None
            data.append(row)
    return data


# =============================================================================
# CROSS-PAPER CONSISTENCY
# =============================================================================

def check_cross_paper_consistency(data: List[Dict], 
                                   threshold: float = 0.3) -> List[Dict]:
    """
    Check if values for the same material are consistent across papers.
    
    Args:
        data: List of extraction rows
        threshold: Maximum allowed relative deviation (default 30%)
        
    Returns:
        List of inconsistency flags
    """
    flags = []
    
    # Group by material name
    material_groups = defaultdict(list)
    for row in data:
        material = row.get('material_name', '').strip()
        if material:
            material_groups[material].append(row)
    
    # Check each material with multiple entries
    for material, rows in material_groups.items():
        if len(rows) < 2:
            continue
        
        # Check key parameters
        for param in ['Ea_eV', 'T_onset_K']:
            values = [r[param] for r in rows if r.get(param) is not None]
            if len(values) < 2:
                continue
            
            mean_val = np.mean(values)
            std_val = np.std(values)
            cv = std_val / mean_val if mean_val > 0 else 0  # Coefficient of variation
            
            if cv > threshold:
                flags.append({
                    'type': 'cross_paper_mismatch',
                    'material': material,
                    'parameter': param,
                    'values': values,
                    'mean': mean_val,
                    'std': std_val,
                    'cv': cv,
                    'sources': [r.get('filename', '') for r in rows if r.get(param)],
                    'message': f"{material}: {param} varies by {cv*100:.1f}% across papers"
                })
    
    return flags


# =============================================================================
# FAMILY CONSISTENCY
# =============================================================================

def check_family_consistency(data: List[Dict]) -> List[Dict]:
    """
    Check if values are consistent with expected ranges for material family.
    
    Returns:
        List of outlier flags
    """
    flags = []
    
    for row in data:
        material = row.get('material_name', '')
        family = row.get('material_family', '').lower()
        
        if family not in FAMILY_EXPECTATIONS:
            continue
        
        expectations = FAMILY_EXPECTATIONS[family]
        
        # Check Ea
        Ea = row.get('Ea_eV')
        if Ea is not None:
            Ea_min, Ea_max = expectations['Ea_eV']
            if Ea < Ea_min * 0.5 or Ea > Ea_max * 2:
                flags.append({
                    'type': 'Ea_family_outlier',
                    'material': material,
                    'family': family,
                    'value': Ea,
                    'expected_range': expectations['Ea_eV'],
                    'source': row.get('filename', ''),
                    'message': f"{material}: Ea={Ea:.3f} eV outside typical range for {family}"
                })
        
        # Check T_onset (approximate, assuming ~100 V/cm)
        T_onset = row.get('T_onset_K')
        E_field = row.get('E_field_Vcm', 100)
        if T_onset is not None:
            # Adjust expected range based on field (rough approximation)
            field_factor = (100 / E_field) ** 0.3 if E_field > 0 else 1
            T_min = expectations['T_onset_typical'][0] * field_factor * 0.7
            T_max = expectations['T_onset_typical'][1] * field_factor * 1.3
            
            if T_onset < T_min or T_onset > T_max:
                flags.append({
                    'type': 'T_onset_family_outlier',
                    'material': material,
                    'family': family,
                    'value': T_onset,
                    'E_field': E_field,
                    'expected_range': (T_min, T_max),
                    'source': row.get('filename', ''),
                    'message': f"{material}: T_onset={T_onset:.0f}K at {E_field}V/cm unusual for {family}"
                })
    
    return flags


# =============================================================================
# ARRHENIUS CONSISTENCY
# =============================================================================

def check_arrhenius_consistency(data: List[Dict]) -> List[Dict]:
    """
    Check if Ea and sigma_0 are physically reasonable.
    
    Returns:
        List of flags
    """
    flags = []
    
    for row in data:
        material = row.get('material_name', '')
        Ea = row.get('Ea_eV')
        sigma_0 = row.get('sigma_0_Sm')
        
        # Check Ea
        if Ea is not None:
            if Ea < 0:
                flags.append({
                    'type': 'negative_Ea',
                    'material': material,
                    'value': Ea,
                    'source': row.get('filename', ''),
                    'message': f"{material}: Negative Ea={Ea:.3f} eV - check sign or data"
                })
            elif Ea > 3.0:
                flags.append({
                    'type': 'high_Ea',
                    'material': material,
                    'value': Ea,
                    'source': row.get('filename', ''),
                    'message': f"{material}: Very high Ea={Ea:.3f} eV - unusual for ceramics"
                })
        
        # Check sigma_0
        if sigma_0 is not None:
            if sigma_0 < 1:
                flags.append({
                    'type': 'low_sigma_0',
                    'material': material,
                    'value': sigma_0,
                    'source': row.get('filename', ''),
                    'message': f"{material}: Very low sigma_0={sigma_0:.2e} S/m - check units"
                })
            elif sigma_0 > 1e10:
                flags.append({
                    'type': 'high_sigma_0',
                    'material': material,
                    'value': sigma_0,
                    'source': row.get('filename', ''),
                    'message': f"{material}: Very high sigma_0={sigma_0:.2e} S/m - verify"
                })
    
    return flags


# =============================================================================
# ONSET CONSISTENCY
# =============================================================================

def check_onset_field_relationship(data: List[Dict]) -> List[Dict]:
    """
    Check if T_onset decreases with increasing E_field for same material.
    
    Returns:
        List of flags
    """
    flags = []
    
    # Group by material
    material_groups = defaultdict(list)
    for row in data:
        material = row.get('material_name', '').strip()
        T_onset = row.get('T_onset_K')
        E_field = row.get('E_field_Vcm')
        
        if material and T_onset and E_field:
            material_groups[material].append({
                'T_onset': T_onset,
                'E_field': E_field,
                'source': row.get('filename', '')
            })
    
    # Check relationship for materials with multiple conditions
    for material, points in material_groups.items():
        if len(points) < 2:
            continue
        
        # Sort by E_field
        sorted_points = sorted(points, key=lambda x: x['E_field'])
        
        # Check if T_onset generally decreases with E_field
        for i in range(len(sorted_points) - 1):
            E1, T1 = sorted_points[i]['E_field'], sorted_points[i]['T_onset']
            E2, T2 = sorted_points[i+1]['E_field'], sorted_points[i+1]['T_onset']
            
            # If E increases significantly but T also increases, flag it
            if E2 > E1 * 1.2 and T2 > T1 * 1.05:
                flags.append({
                    'type': 'inverse_onset_trend',
                    'material': material,
                    'point1': {'E': E1, 'T': T1},
                    'point2': {'E': E2, 'T': T2},
                    'message': f"{material}: T_onset increases with E_field - check data"
                })
    
    return flags


# =============================================================================
# COMPLETENESS CHECK
# =============================================================================

def check_data_completeness(data: List[Dict]) -> Dict:
    """
    Calculate completeness statistics for the dataset.
    
    Returns:
        Dictionary with completeness metrics
    """
    total = len(data)
    if total == 0:
        return {'total': 0}
    
    # Count non-empty fields
    field_counts = defaultdict(int)
    required_fields = ['T_onset_K', 'E_field_Vcm']
    important_fields = ['Ea_eV', 'sigma_0_Sm', 'J_crit_Acm2']
    supplementary_fields = ['delta_H_Jmol', 'delta_S_JmolK', 'grain_size_um']
    
    for row in data:
        for field in required_fields + important_fields + supplementary_fields:
            if row.get(field):
                field_counts[field] += 1
    
    # Calculate scores
    ready_for_solver = sum(1 for row in data 
                          if row.get('T_onset_K') and row.get('E_field_Vcm'))
    
    full_data = sum(1 for row in data
                   if row.get('T_onset_K') and row.get('E_field_Vcm')
                   and row.get('Ea_eV') and row.get('sigma_0_Sm'))
    
    return {
        'total_rows': total,
        'ready_for_solver': ready_for_solver,
        'full_data': full_data,
        'field_coverage': {k: v/total*100 for k, v in field_counts.items()},
        'required_complete': all(field_counts[f] == total for f in required_fields)
    }


# =============================================================================
# REPORT GENERATION
# =============================================================================

def generate_validation_report(data: List[Dict]) -> Dict:
    """
    Generate comprehensive validation report.
    
    Returns:
        Dictionary with all validation results
    """
    report = {
        'summary': {},
        'completeness': check_data_completeness(data),
        'cross_paper_flags': check_cross_paper_consistency(data),
        'family_flags': check_family_consistency(data),
        'arrhenius_flags': check_arrhenius_consistency(data),
        'onset_flags': check_onset_field_relationship(data),
    }
    
    # Calculate summary
    all_flags = (report['cross_paper_flags'] + report['family_flags'] +
                 report['arrhenius_flags'] + report['onset_flags'])
    
    report['summary'] = {
        'total_rows': len(data),
        'total_flags': len(all_flags),
        'cross_paper_issues': len(report['cross_paper_flags']),
        'family_outliers': len(report['family_flags']),
        'arrhenius_issues': len(report['arrhenius_flags']),
        'onset_issues': len(report['onset_flags']),
        'ready_for_solver': report['completeness']['ready_for_solver'],
        'quality_score': max(0, 100 - len(all_flags) * 5)
    }
    
    return report


def print_validation_report(report: Dict):
    """Print formatted validation report."""
    print("\n" + "="*70)
    print("EXTRACTION VALIDATION REPORT")
    print("="*70)
    
    s = report['summary']
    print(f"\nTotal rows: {s['total_rows']}")
    print(f"Ready for solver: {s['ready_for_solver']}")
    print(f"Quality score: {s['quality_score']}/100")
    print(f"\nTotal flags: {s['total_flags']}")
    print(f"  - Cross-paper inconsistencies: {s['cross_paper_issues']}")
    print(f"  - Family outliers: {s['family_outliers']}")
    print(f"  - Arrhenius issues: {s['arrhenius_issues']}")
    print(f"  - Onset trend issues: {s['onset_issues']}")
    
    # Completeness
    print("\n" + "-"*70)
    print("DATA COMPLETENESS")
    print("-"*70)
    c = report['completeness']
    if c.get('field_coverage'):
        for field, pct in sorted(c['field_coverage'].items(), key=lambda x: -x[1]):
            bar = "█" * int(pct/5) + "░" * (20 - int(pct/5))
            print(f"  {field:<20} {bar} {pct:5.1f}%")
    
    # Show flags
    all_flags = (report['cross_paper_flags'] + report['family_flags'] +
                 report['arrhenius_flags'] + report['onset_flags'])
    
    if all_flags:
        print("\n" + "-"*70)
        print("FLAGS AND WARNINGS")
        print("-"*70)
        for flag in all_flags[:20]:
            print(f"  ⚠️  [{flag['type']}] {flag['message']}")
        if len(all_flags) > 20:
            print(f"  ... and {len(all_flags)-20} more flags")
    else:
        print("\n✓ No validation issues found")
    
    print("\n" + "="*70)


def save_validation_report(report: Dict, output_path: str):
    """Save validation report to JSON."""
    # Convert numpy types
    def convert(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (np.integer, np.floating)):
            return float(obj)
        elif isinstance(obj, dict):
            return {k: convert(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert(v) for v in obj]
        return obj
    
    clean_report = convert(report)
    
    with open(output_path, 'w') as f:
        json.dump(clean_report, f, indent=2)
    print(f"\nValidation report saved to: {output_path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Validate extracted Flash data and flag inconsistencies',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python validate_extraction.py data/material_extraction.csv
  python validate_extraction.py data/material_extraction.csv --output report.json

Validation Checks:
  - Cross-paper consistency (same material, different papers)
  - Family consistency (values match expected ranges)
  - Arrhenius consistency (Ea, sigma_0 physically reasonable)
  - Onset trend (T_onset decreases with E_field)
        """
    )
    
    parser.add_argument('csv_file', help='Path to extraction CSV file')
    parser.add_argument('--output', '-o',
                        help='Output JSON file for validation report')
    parser.add_argument('--threshold', '-t', type=float, default=0.3,
                        help='Cross-paper consistency threshold (default: 0.3 = 30%%)')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress console output')
    
    args = parser.parse_args()
    
    # Check input file exists
    if not Path(args.csv_file).exists():
        print(f"Error: File not found: {args.csv_file}")
        sys.exit(1)
    
    try:
        # Load data
        data = load_extraction_data(args.csv_file)
        
        if len(data) == 0:
            print("Warning: No data rows found in CSV")
            sys.exit(0)
        
        # Generate report
        report = generate_validation_report(data)
        
        # Print report
        if not args.quiet:
            print_validation_report(report)
        
        # Save report if output specified
        if args.output:
            save_validation_report(report, args.output)
        
        # Return code based on quality
        if report['summary']['quality_score'] < 50:
            sys.exit(2)  # Low quality warning
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
