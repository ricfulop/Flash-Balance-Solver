#!/usr/bin/env python3
"""
Unit Converter and Validator for Flash Balance Data

Validates and converts extracted data to standard units.
Flags values outside expected ranges for review.

Usage:
    python unit_converter.py <extraction_csv> [options]

Standard Units:
    Temperature: K
    Electric field: V/cm
    Conductivity: S/m
    Current density: A/cm²
    Activation energy: eV
    Enthalpy: J/mol
    Entropy: J/(mol·K)
    Grain size: μm
"""

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import re


# =============================================================================
# UNIT DEFINITIONS AND CONVERSIONS
# =============================================================================

@dataclass
class UnitConversion:
    """Definition of a unit conversion."""
    from_unit: str
    to_unit: str
    factor: float
    offset: float = 0.0
    
    def convert(self, value: float) -> float:
        return value * self.factor + self.offset


# Temperature conversions
TEMP_CONVERSIONS = {
    'C_to_K': UnitConversion('°C', 'K', 1.0, 273.15),
    'F_to_K': UnitConversion('°F', 'K', 5/9, 273.15 - 32*5/9),
}

# Electric field conversions (target: V/cm)
FIELD_CONVERSIONS = {
    'V/m_to_V/cm': UnitConversion('V/m', 'V/cm', 0.01),
    'kV/cm_to_V/cm': UnitConversion('kV/cm', 'V/cm', 1000),
    'kV/m_to_V/cm': UnitConversion('kV/m', 'V/cm', 10),
    'V/mm_to_V/cm': UnitConversion('V/mm', 'V/cm', 10),
}

# Conductivity conversions (target: S/m)
CONDUCTIVITY_CONVERSIONS = {
    'S/cm_to_S/m': UnitConversion('S/cm', 'S/m', 100),
    '1/Ohm.cm_to_S/m': UnitConversion('(Ω·cm)⁻¹', 'S/m', 100),
    '1/Ohm.m_to_S/m': UnitConversion('(Ω·m)⁻¹', 'S/m', 1),
    'mS/cm_to_S/m': UnitConversion('mS/cm', 'S/m', 0.1),
}

# Current density conversions (target: A/cm²)
CURRENT_CONVERSIONS = {
    'mA/mm2_to_A/cm2': UnitConversion('mA/mm²', 'A/cm²', 0.01),
    'A/m2_to_A/cm2': UnitConversion('A/m²', 'A/cm²', 0.0001),
    'mA/cm2_to_A/cm2': UnitConversion('mA/cm²', 'A/cm²', 0.001),
    'kA/m2_to_A/cm2': UnitConversion('kA/m²', 'A/cm²', 0.1),
}

# Activation energy conversions (target: eV)
ENERGY_CONVERSIONS = {
    'kJ/mol_to_eV': UnitConversion('kJ/mol', 'eV', 1/96.485),
    'J/mol_to_eV': UnitConversion('J/mol', 'eV', 1/96485),
    'cal/mol_to_eV': UnitConversion('cal/mol', 'eV', 4.184/96485),
    'kcal/mol_to_eV': UnitConversion('kcal/mol', 'eV', 4.184/96.485),
}

# Enthalpy conversions (target: J/mol)
ENTHALPY_CONVERSIONS = {
    'kJ/mol_to_J/mol': UnitConversion('kJ/mol', 'J/mol', 1000),
    'cal/mol_to_J/mol': UnitConversion('cal/mol', 'J/mol', 4.184),
    'kcal/mol_to_J/mol': UnitConversion('kcal/mol', 'J/mol', 4184),
}

# Grain size conversions (target: μm)
SIZE_CONVERSIONS = {
    'nm_to_um': UnitConversion('nm', 'μm', 0.001),
    'mm_to_um': UnitConversion('mm', 'μm', 1000),
    'm_to_um': UnitConversion('m', 'μm', 1e6),
}


# =============================================================================
# VALIDATION RANGES
# =============================================================================

VALIDATION_RANGES = {
    'T_onset_K': {
        'min': 300, 'max': 2500,
        'typical_min': 500, 'typical_max': 1800,
        'description': 'Flash onset temperature'
    },
    'E_field_Vcm': {
        'min': 1, 'max': 10000,
        'typical_min': 10, 'typical_max': 500,
        'description': 'Electric field'
    },
    'Ea_eV': {
        'min': 0.01, 'max': 5.0,
        'typical_min': 0.1, 'typical_max': 2.0,
        'description': 'Activation energy'
    },
    'sigma_0_Sm': {
        'min': 1e-2, 'max': 1e10,
        'typical_min': 1e2, 'typical_max': 1e7,
        'description': 'Pre-exponential conductivity'
    },
    'J_crit_Acm2': {
        'min': 0.001, 'max': 1000,
        'typical_min': 0.1, 'typical_max': 100,
        'description': 'Critical current density'
    },
    'grain_size_um': {
        'min': 0.001, 'max': 1000,
        'typical_min': 0.01, 'typical_max': 100,
        'description': 'Grain size'
    },
    'delta_H_Jmol': {
        'min': -5e6, 'max': 0,
        'typical_min': -2e6, 'typical_max': -1e5,
        'description': 'Formation enthalpy (should be negative)'
    },
}


# =============================================================================
# UNIT DETECTION
# =============================================================================

def detect_unit_from_header(header: str) -> Tuple[str, str]:
    """
    Detect parameter and unit from column header.
    
    Returns:
        (parameter_name, detected_unit)
    """
    header_lower = header.lower().strip()
    
    # Temperature patterns
    if 't_onset' in header_lower or 'temperature' in header_lower:
        if '°c' in header_lower or 'celsius' in header_lower or '_c' in header_lower:
            return ('temperature', 'C')
        elif 'k' in header_lower or 'kelvin' in header_lower:
            return ('temperature', 'K')
        return ('temperature', 'unknown')
    
    # Electric field patterns
    if 'e_field' in header_lower or 'field' in header_lower or 'e_' in header_lower:
        if 'v/m' in header_lower and 'kv/m' not in header_lower:
            return ('field', 'V/m')
        elif 'kv/cm' in header_lower:
            return ('field', 'kV/cm')
        elif 'v/cm' in header_lower:
            return ('field', 'V/cm')
        return ('field', 'V/cm')  # Default
    
    # Conductivity patterns
    if 'sigma' in header_lower or 'conductivity' in header_lower:
        if 's/cm' in header_lower:
            return ('conductivity', 'S/cm')
        elif 's/m' in header_lower:
            return ('conductivity', 'S/m')
        return ('conductivity', 'S/m')  # Default
    
    # Activation energy patterns
    if 'ea' in header_lower or 'activation' in header_lower:
        if 'kj/mol' in header_lower:
            return ('energy', 'kJ/mol')
        elif 'j/mol' in header_lower and 'kj' not in header_lower:
            return ('energy', 'J/mol')
        elif 'ev' in header_lower:
            return ('energy', 'eV')
        return ('energy', 'eV')  # Default
    
    return ('unknown', 'unknown')


# =============================================================================
# CONVERSION FUNCTIONS
# =============================================================================

def convert_temperature(value: float, from_unit: str) -> Tuple[float, str]:
    """Convert temperature to Kelvin."""
    if from_unit == 'K':
        return value, 'K'
    elif from_unit == 'C':
        return value + 273.15, 'K'
    elif from_unit == 'F':
        return (value - 32) * 5/9 + 273.15, 'K'
    else:
        # Heuristic: if value < 500, assume Celsius
        if value < 500:
            return value + 273.15, 'K (assumed from °C)'
        return value, 'K (assumed)'


def convert_field(value: float, from_unit: str) -> Tuple[float, str]:
    """Convert electric field to V/cm."""
    conversions = {
        'V/cm': 1.0,
        'V/m': 0.01,
        'kV/cm': 1000.0,
        'kV/m': 10.0,
        'V/mm': 10.0,
    }
    factor = conversions.get(from_unit, 1.0)
    return value * factor, 'V/cm'


def convert_conductivity(value: float, from_unit: str) -> Tuple[float, str]:
    """Convert conductivity to S/m."""
    conversions = {
        'S/m': 1.0,
        'S/cm': 100.0,
        '(Ω·cm)⁻¹': 100.0,
        '(Ω·m)⁻¹': 1.0,
        'mS/cm': 0.1,
    }
    factor = conversions.get(from_unit, 1.0)
    return value * factor, 'S/m'


def convert_energy(value: float, from_unit: str) -> Tuple[float, str]:
    """Convert activation energy to eV."""
    conversions = {
        'eV': 1.0,
        'kJ/mol': 1/96.485,
        'J/mol': 1/96485,
        'kcal/mol': 4.184/96.485,
    }
    factor = conversions.get(from_unit, 1.0)
    return value * factor, 'eV'


def convert_current_density(value: float, from_unit: str) -> Tuple[float, str]:
    """Convert current density to A/cm²."""
    conversions = {
        'A/cm²': 1.0,
        'A/cm2': 1.0,
        'mA/mm²': 0.01,
        'mA/mm2': 0.01,
        'A/m²': 0.0001,
        'A/m2': 0.0001,
        'mA/cm²': 0.001,
        'mA/cm2': 0.001,
    }
    factor = conversions.get(from_unit, 1.0)
    return value * factor, 'A/cm²'


# =============================================================================
# VALIDATION
# =============================================================================

@dataclass
class ValidationResult:
    """Result of validating a single value."""
    field: str
    original_value: float
    original_unit: str
    converted_value: float
    converted_unit: str
    is_valid: bool
    is_typical: bool
    flags: List[str]


def validate_value(field: str, value: float) -> ValidationResult:
    """
    Validate a single value against expected ranges.
    
    Args:
        field: Field name (e.g., 'T_onset_K', 'Ea_eV')
        value: Value to validate
        
    Returns:
        ValidationResult with flags
    """
    flags = []
    is_valid = True
    is_typical = True
    
    if field not in VALIDATION_RANGES:
        return ValidationResult(
            field=field,
            original_value=value,
            original_unit='unknown',
            converted_value=value,
            converted_unit='unknown',
            is_valid=True,
            is_typical=True,
            flags=['no_validation_range']
        )
    
    ranges = VALIDATION_RANGES[field]
    
    # Check absolute validity
    if value < ranges['min']:
        is_valid = False
        flags.append(f"below_min ({ranges['min']})")
    elif value > ranges['max']:
        is_valid = False
        flags.append(f"above_max ({ranges['max']})")
    
    # Check if in typical range
    if value < ranges['typical_min']:
        is_typical = False
        flags.append(f"below_typical ({ranges['typical_min']})")
    elif value > ranges['typical_max']:
        is_typical = False
        flags.append(f"above_typical ({ranges['typical_max']})")
    
    return ValidationResult(
        field=field,
        original_value=value,
        original_unit='standard',
        converted_value=value,
        converted_unit='standard',
        is_valid=is_valid,
        is_typical=is_typical,
        flags=flags
    )


def validate_extraction_row(row: Dict) -> List[ValidationResult]:
    """Validate all values in an extraction row."""
    results = []
    
    # Map CSV columns to validation fields
    field_mapping = {
        'T_onset_K': 'T_onset_K',
        'E_field_Vcm': 'E_field_Vcm',
        'Ea_eV': 'Ea_eV',
        'sigma_0_Sm': 'sigma_0_Sm',
        'J_crit_Acm2': 'J_crit_Acm2',
        'grain_size_um': 'grain_size_um',
        'delta_H_Jmol': 'delta_H_Jmol',
    }
    
    for csv_col, val_field in field_mapping.items():
        if csv_col in row and row[csv_col]:
            try:
                value = float(row[csv_col])
                result = validate_value(val_field, value)
                results.append(result)
            except (ValueError, TypeError):
                pass
    
    return results


# =============================================================================
# FILE PROCESSING
# =============================================================================

def process_extraction_csv(filepath: str) -> Tuple[List[Dict], List[Dict]]:
    """
    Process extraction CSV file, validating all values.
    
    Returns:
        (validated_rows, validation_report)
    """
    validated_rows = []
    validation_report = []
    
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        
        for row_num, row in enumerate(reader, start=2):  # Start at 2 (header is 1)
            # Validate row
            validations = validate_extraction_row(row)
            
            # Build report entry
            row_flags = []
            for v in validations:
                if v.flags:
                    row_flags.extend([f"{v.field}: {flag}" for flag in v.flags])
            
            report_entry = {
                'row': row_num,
                'extraction_id': row.get('extraction_id', ''),
                'material': row.get('material_name', ''),
                'is_valid': all(v.is_valid for v in validations),
                'is_typical': all(v.is_typical for v in validations),
                'flags': row_flags,
                'validations': [
                    {
                        'field': v.field,
                        'value': v.original_value,
                        'is_valid': v.is_valid,
                        'is_typical': v.is_typical,
                        'flags': v.flags
                    }
                    for v in validations
                ]
            }
            
            validation_report.append(report_entry)
            validated_rows.append(row)
    
    return validated_rows, validation_report


def print_validation_report(report: List[Dict]):
    """Print validation report to console."""
    print("\n" + "="*70)
    print("VALIDATION REPORT")
    print("="*70)
    
    total = len(report)
    valid = sum(1 for r in report if r['is_valid'])
    typical = sum(1 for r in report if r['is_typical'])
    
    print(f"Total rows: {total}")
    print(f"Valid: {valid}/{total} ({100*valid/total:.1f}%)")
    print(f"Typical: {typical}/{total} ({100*typical/total:.1f}%)")
    print("-"*70)
    
    # Show flagged rows
    flagged = [r for r in report if r['flags']]
    if flagged:
        print(f"\nFlagged rows ({len(flagged)}):")
        print("-"*70)
        for r in flagged[:20]:  # Show first 20
            print(f"Row {r['row']}: {r['material']}")
            for flag in r['flags']:
                print(f"  ⚠️  {flag}")
        if len(flagged) > 20:
            print(f"... and {len(flagged)-20} more")
    else:
        print("\n✓ All values within expected ranges")
    
    print("="*70)


def save_validation_report(report: List[Dict], output_path: str):
    """Save validation report to JSON."""
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)
    print(f"\nValidation report saved to: {output_path}")


# =============================================================================
# COMMAND LINE INTERFACE
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Validate and convert units in extraction CSV',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python unit_converter.py data/material_extraction.csv
  python unit_converter.py data/material_extraction.csv --output validation_log.json
  python unit_converter.py data/material_extraction.csv --quiet

Standard Units:
  Temperature:      K
  Electric field:   V/cm
  Conductivity:     S/m
  Current density:  A/cm²
  Activation energy: eV
  Enthalpy:         J/mol
        """
    )
    
    parser.add_argument('csv_file', help='Path to extraction CSV file')
    parser.add_argument('--output', '-o',
                        help='Output JSON file for validation report')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress console output')
    
    args = parser.parse_args()
    
    # Check input file exists
    if not Path(args.csv_file).exists():
        print(f"Error: File not found: {args.csv_file}")
        sys.exit(1)
    
    try:
        # Process file
        validated_rows, validation_report = process_extraction_csv(args.csv_file)
        
        # Print report
        if not args.quiet:
            print_validation_report(validation_report)
        
        # Save report if output specified
        if args.output:
            save_validation_report(validation_report, args.output)
        
        # Return code based on validation
        if not all(r['is_valid'] for r in validation_report):
            sys.exit(2)  # Some invalid values
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
