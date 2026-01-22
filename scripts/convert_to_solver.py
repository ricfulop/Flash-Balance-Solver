#!/usr/bin/env python3
"""
Convert Extracted Data to Flash Balance Solver Format

Generates:
1. JSON database (extracted_materials.json)
2. Python code for MATERIAL_DATABASE entries
3. Python code for ExperimentalData entries
4. Summary statistics

Usage:
    python convert_to_solver.py <extraction_csv> [options]
"""

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime
from collections import defaultdict


# =============================================================================
# FAMILY DEFAULTS
# =============================================================================

FAMILY_DEFAULTS = {
    'fluorite': {
        'beta': 1.69,
        'alpha_res': 0.15,
        'gamma': 2.0,
        'n_electrons': 4,
        'r_eff_default': 15e-6,
    },
    'rutile': {
        'beta': 1.43,
        'alpha_res': 0.22,
        'gamma': 1.8,
        'n_electrons': 4,
        'r_eff_default': 18e-6,
    },
    'perovskite': {
        'beta': 1.50,
        'alpha_res': 0.28,
        'gamma': 1.4,
        'n_electrons': 4,
        'r_eff_default': 40e-6,
    },
    'spinel': {
        'beta': 1.20,
        'alpha_res': 0.20,
        'gamma': 1.8,
        'n_electrons': 4,
        'r_eff_default': 22e-6,
    },
    'carbide': {
        'beta': 1.20,
        'alpha_res': 0.08,
        'gamma': 2.0,
        'n_electrons': 4,
        'r_eff_default': 6e-6,
    },
    'nitride': {
        'beta': 0.40,
        'alpha_res': 0.12,
        'gamma': 1.5,
        'n_electrons': 3,
        'r_eff_default': 10e-6,
    },
    'garnet': {
        'beta': 1.60,
        'alpha_res': 0.18,
        'gamma': 1.6,
        'n_electrons': 4,
        'r_eff_default': 12e-6,
    },
    'wurtzite': {
        'beta': 1.35,
        'alpha_res': 0.20,
        'gamma': 1.6,
        'n_electrons': 4,
        'r_eff_default': 15e-6,
    },
    'metal': {
        'beta': 0.05,
        'alpha_res': 0.05,
        'gamma': 1.2,
        'n_electrons': 2,
        'r_eff_default': 5e-3,
    },
    'composite_oxide': {
        'beta': 1.40,
        'alpha_res': 0.18,
        'gamma': 1.8,
        'n_electrons': 4,
        'r_eff_default': 20e-6,
    },
}

FAMILY_ENUM_MAP = {
    'fluorite': 'FLUORITE',
    'rutile': 'RUTILE',
    'perovskite': 'PEROVSKITE',
    'spinel': 'SPINEL',
    'carbide': 'CARBIDE',
    'nitride': 'NITRIDE',
    'garnet': 'FLUORITE',  # Use fluorite as closest
    'wurtzite': 'WURTZITE',
    'metal': 'METAL',
    'composite_oxide': 'SPINEL',
}


# =============================================================================
# DATA LOADING
# =============================================================================

def load_extraction_data(filepath: str) -> List[Dict]:
    """Load extraction CSV."""
    data = []
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Convert numeric fields
            numeric_fields = [
                'T_onset_K', 'T_onset_C', 'E_field_Vcm', 'J_crit_Acm2',
                'Ea_eV', 'sigma_0_Sm', 'delta_H_Jmol', 'delta_S_JmolK',
                'grain_size_um', 'completeness_score', 'Ea_R2_fit',
                'green_density_percent', 'final_density_percent',
                'heating_rate_Cmin', 'power_density_Wcm3', 'sigma_at_onset_Sm',
                'n_electrons'
            ]
            for field in numeric_fields:
                if row.get(field):
                    try:
                        row[field] = float(row[field])
                    except ValueError:
                        row[field] = None
            
            # Convert boolean fields
            bool_fields = ['is_composite', 'ready_for_solver']
            for field in bool_fields:
                val = row.get(field, '').lower()
                row[field] = val in ('true', 'yes', '1', 't', 'y')
            
            data.append(row)
    return data


def load_references(filepath: str) -> Dict[str, Dict]:
    """Load references CSV."""
    refs = {}
    if not Path(filepath).exists():
        return refs
    
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get('doi'):
                refs[row['doi']] = row
            elif row.get('filename'):
                refs[row['filename']] = row
    return refs


# =============================================================================
# JSON GENERATION
# =============================================================================

def generate_json_database(data: List[Dict], references: Dict) -> Dict:
    """Generate JSON database from extracted data."""
    materials = []
    experimental_data = []
    refs_list = []
    
    # Track unique materials for MATERIAL_DATABASE
    unique_materials = {}
    
    for row in data:
        material_name = row.get('material_name', '').strip()
        if not material_name:
            continue
        
        # Create material entry
        entry = {
            'id': row.get('extraction_id', f"{material_name}_{len(materials)+1:03d}"),
            'name': material_name,
            'formula': row.get('material_formula', ''),
            'family': row.get('material_family', '').lower(),
            'is_composite': row.get('is_composite', False),
            'source': {
                'doi': row.get('doi', ''),
                'filename': row.get('filename', ''),
                'figure_table_ref': row.get('figure_table_ref', ''),
                'page_num': row.get('page_num', ''),
            },
            'conditions': {
                'grain_size_um': row.get('grain_size_um'),
                'green_density_percent': row.get('green_density_percent'),
                'final_density_percent': row.get('final_density_percent'),
                'atmosphere': row.get('atmosphere', ''),
                'heating_rate_Cmin': row.get('heating_rate_Cmin'),
                'sample_geometry': row.get('sample_geometry', ''),
            },
            'onset': {
                'T_K': row.get('T_onset_K'),
                'T_C': row.get('T_onset_C'),
                'E_Vcm': row.get('E_field_Vcm'),
                'J_Acm2': row.get('J_crit_Acm2'),
                'power_density_Wcm3': row.get('power_density_Wcm3'),
                'condition_type': row.get('condition_type', 'at_onset'),
                'T_confidence': row.get('T_onset_confidence', ''),
                'E_confidence': row.get('E_field_confidence', ''),
            },
            'conductivity': {
                'Ea_eV': row.get('Ea_eV'),
                'sigma_0_Sm': row.get('sigma_0_Sm'),
                'sigma_at_onset_Sm': row.get('sigma_at_onset_Sm'),
                'Ea_R2': row.get('Ea_R2_fit'),
                'Ea_confidence': row.get('Ea_confidence', ''),
                'Ea_source': row.get('Ea_source', ''),
            },
            'thermodynamics': {
                'delta_H_Jmol': row.get('delta_H_Jmol'),
                'delta_S_JmolK': row.get('delta_S_JmolK'),
                'n_electrons': row.get('n_electrons'),
                'delta_H_source': row.get('delta_H_source', ''),
                'delta_S_source': row.get('delta_S_source', ''),
            },
            'quality': {
                'completeness_score': row.get('completeness_score'),
                'ready_for_solver': row.get('ready_for_solver', False),
                'notes': row.get('data_quality_notes', ''),
                'graphs_digitized': row.get('graphs_digitized', ''),
            }
        }
        
        # Remove None values
        entry = clean_none_values(entry)
        materials.append(entry)
        
        # Add to experimental data if has onset
        if row.get('T_onset_K') and row.get('E_field_Vcm'):
            exp_entry = {
                'material': material_name,
                'T_onset_K': row.get('T_onset_K'),
                'E_field_Vcm': row.get('E_field_Vcm'),
                'J_crit_Acm2': row.get('J_crit_Acm2'),
                'source_doi': row.get('doi', ''),
                'source_file': row.get('filename', ''),
            }
            experimental_data.append(exp_entry)
        
        # Track best entry for each material (for MATERIAL_DATABASE)
        if material_name not in unique_materials:
            unique_materials[material_name] = row
        else:
            current_score = row.get('completeness_score') or 0
            existing_score = unique_materials[material_name].get('completeness_score') or 0
            if current_score > existing_score:
                unique_materials[material_name] = row
        
        # Add reference if not already added
        doi = row.get('doi', '')
        filename = row.get('filename', '')
        ref_key = doi or filename
        if ref_key and ref_key not in [r.get('doi') or r.get('filename') for r in refs_list]:
            refs_list.append({
                'doi': doi,
                'filename': filename,
            })
    
    return {
        'version': '1.0',
        'extraction_date': datetime.now().isoformat(),
        'total_entries': len(materials),
        'unique_materials': len(unique_materials),
        'materials': materials,
        'experimental_data': experimental_data,
        'references': refs_list,
        'unique_material_keys': list(unique_materials.keys()),
    }


def clean_none_values(d: Dict) -> Dict:
    """Recursively remove None values from dict."""
    if not isinstance(d, dict):
        return d
    return {k: clean_none_values(v) for k, v in d.items() 
            if v is not None and v != '' and v != {}}


# =============================================================================
# PYTHON CODE GENERATION
# =============================================================================

def generate_material_database_entry(row: Dict) -> str:
    """Generate Python code for a single MaterialParameters entry."""
    material_name = row.get('material_name', 'Unknown')
    family = row.get('material_family', 'spinel').lower()
    
    # Get family defaults
    defaults = FAMILY_DEFAULTS.get(family, FAMILY_DEFAULTS['spinel'])
    family_enum = FAMILY_ENUM_MAP.get(family, 'SPINEL')
    
    # Extract values or use defaults
    Ea = row.get('Ea_eV') or 0.8
    sigma_0 = row.get('sigma_0_Sm') or 1e4
    delta_H = row.get('delta_H_Jmol') or -500000
    delta_S = row.get('delta_S_JmolK') or -150
    n_electrons = row.get('n_electrons') or defaults['n_electrons']
    
    # Generate comment with source info
    doi = row.get('doi', '')
    filename = row.get('filename', '')
    completeness = row.get('completeness_score', 0)
    Ea_source = row.get('Ea_source', 'estimated')
    
    source_comment = f"# Source: {filename}"
    if doi:
        source_comment += f", DOI: {doi}"
    source_comment += f"\n    # Completeness: {completeness}%"
    
    # Format the entry
    code = f'''    "{material_name}": MaterialParameters(
        {source_comment}
        name="{material_name}",
        family=MaterialFamily.{family_enum},
        Ea={Ea:.3f},           # {Ea_source}
        sigma_0={sigma_0:.2e},
        beta={defaults['beta']:.2f},         # Family default
        alpha_res={defaults['alpha_res']:.2f},    # Family default
        gamma={defaults['gamma']:.1f},
        delta_H={int(delta_H)},
        delta_S={int(delta_S)},
        n_electrons={int(n_electrons)},
        r_eff={defaults['r_eff_default']:.0e},      # Family default, may need calibration
    ),'''
    
    return code


def generate_experimental_data_entry(row: Dict) -> str:
    """Generate Python code for ExperimentalData entry."""
    material = row.get('material_name', 'Unknown')
    T_onset = row.get('T_onset_K', 0)
    E_field = row.get('E_field_Vcm', 0)
    filename = row.get('filename', '')
    
    # Generate ref number (placeholder)
    return f'    ExperimentalData("{material}", {int(T_onset)}, {int(E_field)}, REF_NUM),  # {filename}'


def generate_python_code(data: List[Dict]) -> str:
    """Generate complete Python code for solver integration."""
    # Group by material and get best entry for each
    unique_materials = {}
    for row in data:
        material = row.get('material_name', '').strip()
        if not material:
            continue
        if material not in unique_materials:
            unique_materials[material] = row
        else:
            current_score = row.get('completeness_score') or 0
            existing_score = unique_materials[material].get('completeness_score') or 0
            if current_score > existing_score:
                unique_materials[material] = row
    
    # Generate material database entries
    material_entries = []
    for material, row in sorted(unique_materials.items()):
        if row.get('Ea_eV') or row.get('T_onset_K'):  # Has some useful data
            material_entries.append(generate_material_database_entry(row))
    
    # Generate experimental data entries
    exp_entries = []
    for row in data:
        if row.get('T_onset_K') and row.get('E_field_Vcm'):
            exp_entries.append(generate_experimental_data_entry(row))
    
    # Assemble code
    code = f'''"""
Auto-generated Flash Balance Solver entries
Generated: {datetime.now().isoformat()}
Source: material_extraction.csv

To use:
1. Add MaterialParameters entries to MATERIAL_DATABASE in flash_balance_solver.py
2. Add ExperimentalData entries to EXPERIMENTAL_DATA
3. Update REFERENCES dict with proper ref numbers
4. Run validation to check accuracy
"""

# =============================================================================
# MATERIAL DATABASE ENTRIES
# Add these to MATERIAL_DATABASE in flash_balance_solver.py
# =============================================================================

EXTRACTED_MATERIALS = {{
{chr(10).join(material_entries)}
}}


# =============================================================================
# EXPERIMENTAL DATA ENTRIES  
# Add these to EXPERIMENTAL_DATA in flash_balance_solver.py
# Replace REF_NUM with actual reference numbers
# =============================================================================

EXTRACTED_EXPERIMENTAL_DATA = [
{chr(10).join(exp_entries)}
]


# =============================================================================
# SUMMARY
# =============================================================================
# Total unique materials: {len(unique_materials)}
# Total experimental data points: {len(exp_entries)}
'''
    
    return code


# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

def generate_summary(data: List[Dict], json_db: Dict) -> str:
    """Generate markdown summary of extraction."""
    # Count by family
    family_counts = defaultdict(int)
    for row in data:
        family = row.get('material_family', 'unknown').lower()
        family_counts[family] += 1
    
    # Count ready for solver
    ready = sum(1 for row in data if row.get('T_onset_K') and row.get('E_field_Vcm'))
    full_data = sum(1 for row in data if row.get('T_onset_K') and row.get('Ea_eV'))
    
    # Unique materials
    unique = len(set(row.get('material_name', '') for row in data if row.get('material_name')))
    
    summary = f"""# Extraction Summary

Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}

## Overview

| Metric | Count |
|--------|-------|
| Total data points | {len(data)} |
| Unique materials | {unique} |
| Ready for solver (T_onset + E_field) | {ready} |
| Full data (+ Ea + σ₀) | {full_data} |

## By Material Family

| Family | Count |
|--------|-------|
"""
    for family, count in sorted(family_counts.items(), key=lambda x: -x[1]):
        summary += f"| {family} | {count} |\n"
    
    summary += f"""
## Data Quality

- Entries with onset data: {ready} ({100*ready/max(1,len(data)):.1f}%)
- Entries with conductivity: {full_data} ({100*full_data/max(1,len(data)):.1f}%)

## Files Generated

- `extracted_materials.json` - Full JSON database
- `generated_solver_entries.py` - Python code for solver
- `extraction_summary.md` - This summary

## Next Steps

1. Review generated Python entries in `generated_solver_entries.py`
2. Add entries to `flash_balance_solver.py`
3. Run validation: `python flash_balance_solver.py`
4. Calibrate `r_eff` values if predictions are off
"""
    
    return summary


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Convert extracted data to Flash Balance Solver format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python convert_to_solver.py data/material_extraction.csv
  python convert_to_solver.py data/material_extraction.csv --output-dir output/

Outputs:
  - extracted_materials.json: Full JSON database
  - generated_solver_entries.py: Python code for MATERIAL_DATABASE
  - extraction_summary.md: Summary statistics
        """
    )
    
    parser.add_argument('csv_file', help='Path to extraction CSV file')
    parser.add_argument('--output-dir', '-o', default='data',
                        help='Output directory (default: data/)')
    parser.add_argument('--references', '-r',
                        help='Path to references CSV file')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress console output')
    
    args = parser.parse_args()
    
    # Check input file exists
    if not Path(args.csv_file).exists():
        print(f"Error: File not found: {args.csv_file}")
        sys.exit(1)
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Load data
        data = load_extraction_data(args.csv_file)
        
        if len(data) == 0:
            print("Warning: No data rows found in CSV")
            sys.exit(0)
        
        # Load references if provided
        refs = {}
        if args.references:
            refs = load_references(args.references)
        
        # Generate JSON database
        json_db = generate_json_database(data, refs)
        json_path = output_dir / 'extracted_materials.json'
        with open(json_path, 'w') as f:
            json.dump(json_db, f, indent=2)
        
        # Generate Python code
        python_code = generate_python_code(data)
        python_path = output_dir / 'generated_solver_entries.py'
        with open(python_path, 'w') as f:
            f.write(python_code)
        
        # Generate summary
        summary = generate_summary(data, json_db)
        summary_path = output_dir / 'extraction_summary.md'
        with open(summary_path, 'w') as f:
            f.write(summary)
        
        if not args.quiet:
            print("\n" + "="*60)
            print("CONVERSION COMPLETE")
            print("="*60)
            print(f"\nProcessed {len(data)} data points")
            print(f"Unique materials: {json_db['unique_materials']}")
            print(f"\nFiles generated:")
            print(f"  - {json_path}")
            print(f"  - {python_path}")
            print(f"  - {summary_path}")
            print("\n" + "="*60)
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
