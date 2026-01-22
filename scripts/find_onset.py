#!/usr/bin/env python3
"""
Find Flash Onset Temperature

Detect flash onset temperature from current-temperature or power-temperature data.
Uses derivative analysis to find the inflection point indicating flash ignition.

Usage:
    python find_onset.py <csv_file> [options]

Input CSV format:
    Columns: temperature, current/power
    Temperature can be in K or °C (auto-detected or specified)
    Current can be in A, mA, A/cm², mA/mm², etc.

Output:
    - T_onset (K and °C): Flash onset temperature
    - J_crit or I_crit: Current/current density at onset
    - Confidence level based on sharpness of transition
"""

import argparse
import numpy as np
import csv
import json
import sys
from pathlib import Path
from typing import Tuple, Dict, Optional, List
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d


def read_csv_data(filepath: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Read CSV data.
    
    Returns:
        x_data: Temperature values
        y_data: Current/power values
    """
    x_data = []
    y_data = []
    
    with open(filepath, 'r') as f:
        reader = csv.reader(f)
        header = next(reader, None)
        
        for row in reader:
            if len(row) >= 2:
                try:
                    x_data.append(float(row[0]))
                    y_data.append(float(row[1]))
                except ValueError:
                    continue
    
    return np.array(x_data), np.array(y_data)


def detect_temperature_unit(T_data: np.ndarray) -> str:
    """
    Auto-detect if temperature is in K or °C.
    
    Heuristic: If max < 500, likely °C; if max > 500, likely K
    """
    T_max = np.max(T_data)
    if T_max > 500:
        return 'K'
    else:
        return 'C'


def convert_to_kelvin(T_data: np.ndarray, unit: str) -> np.ndarray:
    """Convert temperature to Kelvin."""
    if unit == 'C':
        return T_data + 273.15
    return T_data


def smooth_data(y_data: np.ndarray, method: str = 'savgol', 
                window: int = 5) -> np.ndarray:
    """
    Smooth noisy data before differentiation.
    
    Args:
        y_data: Raw y values
        method: 'savgol' or 'gaussian'
        window: Window size for smoothing
        
    Returns:
        Smoothed y values
    """
    if len(y_data) < window:
        return y_data
    
    if method == 'savgol':
        # Ensure window is odd
        if window % 2 == 0:
            window += 1
        return savgol_filter(y_data, window, 2)
    elif method == 'gaussian':
        sigma = window / 3
        return gaussian_filter1d(y_data, sigma)
    else:
        return y_data


def find_onset_derivative(T: np.ndarray, y: np.ndarray, 
                          smooth: bool = True) -> Dict:
    """
    Find onset using derivative analysis.
    
    The onset is identified as the point of maximum first derivative
    (steepest rise) in the current/power vs temperature curve.
    
    Args:
        T: Temperature array (K)
        y: Current/power array
        smooth: Whether to smooth data before differentiation
        
    Returns:
        Dictionary with onset temperature and related metrics
    """
    # Sort by temperature
    sort_idx = np.argsort(T)
    T = T[sort_idx]
    y = y[sort_idx]
    
    # Smooth if requested
    if smooth and len(y) > 5:
        y_smooth = smooth_data(y, method='savgol', window=min(7, len(y)//2*2+1))
    else:
        y_smooth = y
    
    # Calculate first derivative (dy/dT)
    dT = np.diff(T)
    dy = np.diff(y_smooth)
    
    # Avoid division by zero
    dT[dT == 0] = 1e-10
    
    derivative = dy / dT
    T_mid = (T[:-1] + T[1:]) / 2
    
    # Find maximum derivative (steepest rise)
    max_deriv_idx = np.argmax(derivative)
    T_onset = T_mid[max_deriv_idx]
    
    # Get y value at onset (interpolate)
    y_onset = np.interp(T_onset, T, y)
    
    # Calculate confidence metrics
    max_derivative = derivative[max_deriv_idx]
    mean_derivative = np.mean(np.abs(derivative))
    derivative_ratio = max_derivative / mean_derivative if mean_derivative > 0 else 0
    
    # Determine confidence based on how sharp the transition is
    if derivative_ratio > 5:
        confidence = 'H'
        confidence_desc = 'Sharp transition (derivative ratio > 5)'
    elif derivative_ratio > 2:
        confidence = 'M'
        confidence_desc = 'Clear transition (derivative ratio > 2)'
    else:
        confidence = 'L'
        confidence_desc = 'Gradual transition, onset may be uncertain'
    
    return {
        'T_onset_K': T_onset,
        'T_onset_C': T_onset - 273.15,
        'y_at_onset': y_onset,
        'max_derivative': max_derivative,
        'derivative_ratio': derivative_ratio,
        'confidence': confidence,
        'confidence_description': confidence_desc
    }


def find_onset_threshold(T: np.ndarray, y: np.ndarray, 
                         threshold_fraction: float = 0.1) -> Dict:
    """
    Find onset using threshold method.
    
    The onset is defined as the point where y exceeds a threshold
    above the baseline (e.g., 10% of the range).
    
    Args:
        T: Temperature array (K)
        y: Current/power array
        threshold_fraction: Fraction above baseline to define onset
        
    Returns:
        Dictionary with onset temperature
    """
    # Sort by temperature
    sort_idx = np.argsort(T)
    T = T[sort_idx]
    y = y[sort_idx]
    
    # Find baseline (first few points or minimum)
    n_baseline = max(3, len(y) // 10)
    baseline = np.mean(y[:n_baseline])
    
    # Define threshold
    y_range = np.max(y) - baseline
    threshold = baseline + threshold_fraction * y_range
    
    # Find first point above threshold
    above_threshold = np.where(y > threshold)[0]
    
    if len(above_threshold) == 0:
        return {'T_onset_K': None, 'error': 'No points above threshold'}
    
    onset_idx = above_threshold[0]
    
    # Interpolate for more precise onset
    if onset_idx > 0:
        T_below = T[onset_idx - 1]
        T_above = T[onset_idx]
        y_below = y[onset_idx - 1]
        y_above = y[onset_idx]
        
        # Linear interpolation
        if y_above != y_below:
            fraction = (threshold - y_below) / (y_above - y_below)
            T_onset = T_below + fraction * (T_above - T_below)
        else:
            T_onset = T_above
    else:
        T_onset = T[onset_idx]
    
    return {
        'T_onset_K': T_onset,
        'T_onset_C': T_onset - 273.15,
        'y_at_onset': threshold,
        'baseline': baseline,
        'threshold': threshold,
        'confidence': 'M',
        'confidence_description': 'Threshold-based detection'
    }


def analyze_onset_data(filepath: str,
                       temp_unit: Optional[str] = None,
                       y_label: str = 'current',
                       method: str = 'derivative') -> Dict:
    """
    Complete analysis of onset data from CSV file.
    
    Args:
        filepath: Path to CSV file
        temp_unit: Temperature unit ('K' or 'C'), auto-detected if None
        y_label: Label for y-axis data ('current', 'power', 'J')
        method: Detection method ('derivative' or 'threshold')
        
    Returns:
        Dictionary with analysis results
    """
    # Read data
    T_raw, y_data = read_csv_data(filepath)
    
    if len(T_raw) < 5:
        raise ValueError(f"Insufficient data points ({len(T_raw)}), need at least 5")
    
    # Detect/convert temperature unit
    if temp_unit is None:
        temp_unit = detect_temperature_unit(T_raw)
    
    T_K = convert_to_kelvin(T_raw, temp_unit)
    
    # Find onset using selected method
    if method == 'derivative':
        results = find_onset_derivative(T_K, y_data)
    elif method == 'threshold':
        results = find_onset_threshold(T_K, y_data)
    else:
        raise ValueError(f"Unknown method: {method}")
    
    # Add metadata
    results['input_file'] = str(filepath)
    results['detected_temp_unit'] = temp_unit
    results['y_label'] = y_label
    results['method'] = method
    results['n_points'] = len(T_raw)
    results['T_range_K'] = (float(np.min(T_K)), float(np.max(T_K)))
    results['y_range'] = (float(np.min(y_data)), float(np.max(y_data)))
    
    return results


def print_results(results: Dict):
    """Print formatted results to console."""
    print("\n" + "="*60)
    print("FLASH ONSET ANALYSIS RESULTS")
    print("="*60)
    print(f"Input file: {results['input_file']}")
    print(f"Data points: {results['n_points']}")
    print(f"Temperature range: {results['T_range_K'][0]:.0f} - {results['T_range_K'][1]:.0f} K")
    print(f"Detection method: {results['method']}")
    print("-"*60)
    
    if results.get('T_onset_K') is not None:
        print(f"Onset Temperature:  {results['T_onset_K']:.1f} K")
        print(f"                    {results['T_onset_C']:.1f} °C")
        print(f"{results['y_label']} at onset: {results['y_at_onset']:.4g}")
        print("-"*60)
        print(f"Confidence level:   {results['confidence']} - {results['confidence_description']}")
        
        if 'derivative_ratio' in results:
            print(f"Derivative ratio:   {results['derivative_ratio']:.2f}")
    else:
        print("ERROR: Could not determine onset temperature")
        if 'error' in results:
            print(f"Reason: {results['error']}")
    
    print("="*60)


def save_results(results: Dict, output_path: str):
    """Save results to JSON file."""
    # Convert numpy types to Python types for JSON serialization
    clean_results = {}
    for k, v in results.items():
        if isinstance(v, np.ndarray):
            clean_results[k] = v.tolist()
        elif isinstance(v, (np.integer, np.floating)):
            clean_results[k] = float(v)
        else:
            clean_results[k] = v
    
    with open(output_path, 'w') as f:
        json.dump(clean_results, f, indent=2)
    print(f"\nResults saved to: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Find flash onset temperature from current/power vs temperature data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python find_onset.py onset_data.csv
  python find_onset.py data.csv --temp-unit C --label "J (A/cm²)"
  python find_onset.py data.csv --method threshold --output results.json

Methods:
  derivative: Find point of maximum slope (steepest rise)
  threshold:  Find point where signal exceeds threshold above baseline
        """
    )
    
    parser.add_argument('csv_file', help='Path to CSV file with T vs current/power data')
    parser.add_argument('--temp-unit', '-t',
                        choices=['K', 'C'],
                        help='Temperature unit (auto-detected if not specified)')
    parser.add_argument('--label', '-l',
                        default='current',
                        help='Label for y-axis data (default: current)')
    parser.add_argument('--method', '-m',
                        choices=['derivative', 'threshold'],
                        default='derivative',
                        help='Detection method (default: derivative)')
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
        results = analyze_onset_data(
            args.csv_file,
            temp_unit=args.temp_unit,
            y_label=args.label,
            method=args.method
        )
        
        # Print results
        if not args.quiet:
            print_results(results)
        
        # Save results if output specified
        if args.output:
            save_results(results, args.output)
        
        # Return code based on success
        if results.get('T_onset_K') is None:
            sys.exit(2)
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
