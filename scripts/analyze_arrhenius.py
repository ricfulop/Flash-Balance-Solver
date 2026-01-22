#!/usr/bin/env python3
"""
Analyze Arrhenius Plot Data

Extract activation energy (Ea) and pre-exponential factor (sigma_0)
from digitized Arrhenius plot data.

Usage:
    python analyze_arrhenius.py <csv_file> [options]

Input CSV format (choose one):
    Option 1: T_K, sigma_Sm (Temperature in K, conductivity in S/m)
    Option 2: inv_T, ln_sigma (1/T in 1/K, ln(sigma) in ln(S/m))
    Option 3: inv_T_1000, log_sigma (1000/T in 1000/K, log10(sigma))

Output:
    - Ea (eV): Activation energy
    - sigma_0 (S/m): Pre-exponential factor
    - R2: Coefficient of determination
    - Confidence level based on R2
"""

import argparse
import numpy as np
import csv
import json
import sys
from pathlib import Path
from typing import Tuple, Dict, Optional

# Physical constants
kB_eV = 8.617e-5  # Boltzmann constant in eV/K


def read_csv_data(filepath: str) -> Tuple[np.ndarray, np.ndarray, str]:
    """
    Read CSV data and detect format.
    
    Returns:
        x_data: Array of x values
        y_data: Array of y values
        format_type: Detected format ('T_sigma', 'inv_T_ln', 'inv_T_1000_log')
    """
    x_data = []
    y_data = []
    
    with open(filepath, 'r') as f:
        reader = csv.reader(f)
        header = next(reader, None)
        
        # Detect format from header if available
        if header:
            header_lower = [h.lower().strip() for h in header]
            if 't_k' in header_lower or 'temperature' in header_lower:
                format_type = 'T_sigma'
            elif '1000/t' in header_lower or '1000_t' in header_lower:
                format_type = 'inv_T_1000_log'
            elif '1/t' in header_lower or 'inv_t' in header_lower:
                format_type = 'inv_T_ln'
            else:
                # Default: assume raw T and sigma
                format_type = 'T_sigma'
        else:
            format_type = 'T_sigma'
        
        for row in reader:
            if len(row) >= 2:
                try:
                    x_data.append(float(row[0]))
                    y_data.append(float(row[1]))
                except ValueError:
                    continue
    
    return np.array(x_data), np.array(y_data), format_type


def convert_to_arrhenius_coords(x: np.ndarray, y: np.ndarray, 
                                 format_type: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert data to Arrhenius coordinates (1/T vs ln(sigma)).
    
    Args:
        x: x-axis data
        y: y-axis data (conductivity or log/ln conductivity)
        format_type: Input format type
        
    Returns:
        inv_T: 1/T in K^-1
        ln_sigma: ln(sigma) where sigma is in S/m
    """
    if format_type == 'T_sigma':
        # Input: T(K), sigma(S/m)
        inv_T = 1.0 / x
        ln_sigma = np.log(y)
    elif format_type == 'inv_T_1000_log':
        # Input: 1000/T, log10(sigma)
        inv_T = x / 1000.0
        ln_sigma = y * np.log(10)  # Convert log10 to ln
    elif format_type == 'inv_T_ln':
        # Input: 1/T, ln(sigma)
        inv_T = x
        ln_sigma = y
    else:
        raise ValueError(f"Unknown format type: {format_type}")
    
    return inv_T, ln_sigma


def fit_arrhenius(inv_T: np.ndarray, ln_sigma: np.ndarray) -> Dict:
    """
    Perform linear fit to Arrhenius data.
    
    ln(sigma) = ln(sigma_0) - Ea/(kB*T)
    
    Args:
        inv_T: 1/T values in K^-1
        ln_sigma: ln(sigma) values
        
    Returns:
        Dictionary with Ea, sigma_0, R2, and fit quality metrics
    """
    # Linear fit: ln(sigma) = intercept + slope * (1/T)
    # slope = -Ea/kB
    # intercept = ln(sigma_0)
    
    coeffs = np.polyfit(inv_T, ln_sigma, 1)
    slope = coeffs[0]
    intercept = coeffs[1]
    
    # Calculate R^2
    y_pred = slope * inv_T + intercept
    ss_res = np.sum((ln_sigma - y_pred) ** 2)
    ss_tot = np.sum((ln_sigma - np.mean(ln_sigma)) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    # Extract parameters
    Ea_eV = -slope * kB_eV  # Ea in eV
    sigma_0 = np.exp(intercept)  # sigma_0 in S/m
    
    # Determine confidence based on R^2
    if r_squared >= 0.99:
        confidence = 'H'  # High
        confidence_desc = 'Excellent fit (R² ≥ 0.99)'
    elif r_squared >= 0.95:
        confidence = 'M'  # Medium
        confidence_desc = 'Good fit (R² ≥ 0.95)'
    elif r_squared >= 0.90:
        confidence = 'L'  # Low
        confidence_desc = 'Acceptable fit (R² ≥ 0.90)'
    else:
        confidence = 'L'  # Low
        confidence_desc = f'Poor fit (R² = {r_squared:.3f}), check data quality'
    
    return {
        'Ea_eV': Ea_eV,
        'sigma_0_Sm': sigma_0,
        'R2': r_squared,
        'slope': slope,
        'intercept': intercept,
        'confidence': confidence,
        'confidence_description': confidence_desc,
        'n_points': len(inv_T),
        'T_range_K': (1/np.max(inv_T), 1/np.min(inv_T))
    }


def analyze_arrhenius_data(filepath: str, 
                           format_override: Optional[str] = None,
                           sigma_unit: str = 'S/m') -> Dict:
    """
    Complete analysis of Arrhenius data from CSV file.
    
    Args:
        filepath: Path to CSV file
        format_override: Override automatic format detection
        sigma_unit: Unit of conductivity ('S/m' or 'S/cm')
        
    Returns:
        Dictionary with analysis results
    """
    # Read data
    x_data, y_data, detected_format = read_csv_data(filepath)
    
    if len(x_data) < 3:
        raise ValueError(f"Insufficient data points ({len(x_data)}), need at least 3")
    
    format_type = format_override if format_override else detected_format
    
    # Convert to Arrhenius coordinates
    inv_T, ln_sigma = convert_to_arrhenius_coords(x_data, y_data, format_type)
    
    # Apply unit conversion if needed
    if sigma_unit == 'S/cm':
        ln_sigma = ln_sigma + np.log(100)  # Convert S/cm to S/m
    
    # Fit Arrhenius equation
    results = fit_arrhenius(inv_T, ln_sigma)
    
    # Add metadata
    results['input_file'] = str(filepath)
    results['detected_format'] = detected_format
    results['used_format'] = format_type
    results['sigma_unit'] = sigma_unit
    
    return results


def print_results(results: Dict):
    """Print formatted results to console."""
    print("\n" + "="*60)
    print("ARRHENIUS ANALYSIS RESULTS")
    print("="*60)
    print(f"Input file: {results['input_file']}")
    print(f"Data points: {results['n_points']}")
    print(f"Temperature range: {results['T_range_K'][0]:.0f} - {results['T_range_K'][1]:.0f} K")
    print("-"*60)
    print(f"Activation Energy (Ea): {results['Ea_eV']:.4f} eV")
    print(f"                        {results['Ea_eV'] * 96485 / 1000:.2f} kJ/mol")
    print(f"Pre-exponential (σ₀):   {results['sigma_0_Sm']:.2e} S/m")
    print(f"                        {results['sigma_0_Sm'] / 100:.2e} S/cm")
    print("-"*60)
    print(f"R² (goodness of fit):   {results['R2']:.6f}")
    print(f"Confidence level:       {results['confidence']} - {results['confidence_description']}")
    print("="*60)
    
    # Warnings
    if results['Ea_eV'] < 0:
        print("⚠️  WARNING: Negative Ea suggests wrong slope sign or data order")
    if results['Ea_eV'] > 3.0:
        print("⚠️  WARNING: Ea > 3 eV is unusually high, verify data")
    if results['Ea_eV'] < 0.05 and results['Ea_eV'] > 0:
        print("⚠️  WARNING: Ea < 0.05 eV suggests metallic conduction")
    if results['R2'] < 0.90:
        print("⚠️  WARNING: Poor fit, data may not follow Arrhenius behavior")


def save_results(results: Dict, output_path: str):
    """Save results to JSON file."""
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Extract Ea and sigma_0 from digitized Arrhenius plot data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python analyze_arrhenius.py data.csv
  python analyze_arrhenius.py data.csv --format inv_T_1000_log --unit S/cm
  python analyze_arrhenius.py data.csv --output results.json

CSV Format Options:
  T_sigma:        Columns are T(K), sigma(S/m)
  inv_T_ln:       Columns are 1/T(K^-1), ln(sigma)
  inv_T_1000_log: Columns are 1000/T(K^-1), log10(sigma)
        """
    )
    
    parser.add_argument('csv_file', help='Path to CSV file with Arrhenius data')
    parser.add_argument('--format', '-f', 
                        choices=['T_sigma', 'inv_T_ln', 'inv_T_1000_log'],
                        help='Override automatic format detection')
    parser.add_argument('--unit', '-u', 
                        choices=['S/m', 'S/cm'], default='S/m',
                        help='Unit of conductivity in input data (default: S/m)')
    parser.add_argument('--output', '-o', 
                        help='Output JSON file for results')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress console output')
    
    args = parser.parse_args()
    
    # Check input file exists
    if not Path(args.csv_file).exists():
        print(f"Error: File not found: {args.csv_file}")
        sys.exit(1)
    
    try:
        # Analyze data
        results = analyze_arrhenius_data(
            args.csv_file,
            format_override=args.format,
            sigma_unit=args.unit
        )
        
        # Print results
        if not args.quiet:
            print_results(results)
        
        # Save results if output specified
        if args.output:
            save_results(results, args.output)
        
        # Return code based on confidence
        if results['R2'] < 0.90:
            sys.exit(2)  # Poor fit warning
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
