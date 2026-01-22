#!/usr/bin/env python3
"""
Generate validation table excluding metals.
"""

from flash_balance_solver import validate_solver
import numpy as np

# Get validation results and filter out metals
results = validate_solver()
results_no_metals = [r for r in results if r['family'] != 'metal']

# Print validation table
print('\n' + '='*105)
print('FLASH BALANCE SOLVER - ACCURACY VALIDATION TABLE (NO METALS)')
print('='*105)
print(f"{'Material':<10} {'Family':<12} {'E (V/cm)':<10} {'T_exp (K)':<10} "
      f"{'T_pred (K)':<11} {'Error (K)':<10} {'Error (%)':<10} {'Status':<8}")
print('-'*105)

total_error = 0
count = 0
accurate_count = 0
family_errors = {}

for r in results_no_metals:
    if r['T_predicted'] is not None:
        status = 'PASS' if r['accuracy_met'] else 'FAIL'
        symbol = '✓' if r['accuracy_met'] else '✗'
        print(f"{r['material']:<10} {r['family']:<12} {r['E_field_Vcm']:<10.0f} "
              f"{r['T_experimental']:<10.0f} {r['T_predicted']:<11.0f} "
              f"{r['error_K']:<+10.0f} {r['error_percent']:<+10.1f} {symbol} {status:<6}")

        total_error += abs(r['error_percent'])
        count += 1
        if r['accuracy_met']:
            accurate_count += 1

        fam = r['family']
        if fam not in family_errors:
            family_errors[fam] = []
        family_errors[fam].append(abs(r['error_percent']))
    else:
        print(f"{r['material']:<10} {r['family']:<12} {r['E_field_Vcm']:<10.0f} "
              f"{r['T_experimental']:<10.0f} {'N/A':<11} {'N/A':<10} {'N/A':<10} ✗ FAIL")

print('-'*105)

# Summary
print(f"\n{'='*60}")
print('SUMMARY')
print(f"{'='*60}")
if count > 0:
    avg_error = total_error / count
    accuracy_rate = 100 * accurate_count / count
    print(f'Total samples validated:     {count}')
    print(f'Average absolute error:      {avg_error:.1f}%')
    print(f'Samples within ±15%:         {accurate_count}/{count} ({accuracy_rate:.0f}%)')

    print(f"\n{'Error by Material Family:':<30}")
    print('-'*40)
    for fam, errors in sorted(family_errors.items()):
        avg = np.mean(errors)
        print(f'  {fam:<20} {avg:>6.1f}% (n={len(errors)})')

print('='*105)

# Paper table
print('\n' + '='*100)
print('TABLE: Flash Balance Model Validation (No Metals)')
print('='*100)

print(f"\n{'Material':<10} {'Family':<12} {'E (V/cm)':<10} {'T_pred (°C)':<12} {'T_exp (°C)':<12} {'Error (%)':<10} {'Ref':<6}")
print('-'*100)

families_order = ['fluorite', 'rutile', 'perovskite', 'spinel', 'nitride', 'carbide']

for family in families_order:
    family_results = [r for r in results_no_metals if r['family'] == family]
    if not family_results:
        continue

    for r in family_results:
        T_exp_C = r['T_experimental'] - 273
        if r['T_predicted'] is not None:
            T_pred_C = r['T_predicted'] - 273
            error_str = f"{r['error_percent']:+.1f}"
        else:
            T_pred_C = 'N/A'
            error_str = 'N/A'

        print(f"{r['material']:<10} {r['family']:<12} {r['E_field_Vcm']:<10.0f} "
              f"{T_pred_C if isinstance(T_pred_C, str) else f'{T_pred_C:.0f}':<12} "
              f"{T_exp_C:<12.0f} {error_str:<10} [{r['ref_num']}]")

print('-'*100)

# Summary statistics
valid_results = [r for r in results_no_metals if r['T_predicted'] is not None]
if valid_results:
    avg_error = np.mean([abs(r['error_percent']) for r in valid_results])
    within_10 = sum(1 for r in valid_results if abs(r['error_percent']) < 10)
    within_15 = sum(1 for r in valid_results if abs(r['error_percent']) < 15)
    print(f"\nSummary: {len(valid_results)} materials | Avg error: {avg_error:.1f}% | "
          f"Within ±10%: {within_10}/{len(valid_results)} | Within ±15%: {within_15}/{len(valid_results)}")

print('\n' + '='*100)

# LaTeX table
print('\n\nLaTeX Table Format:')
print('-'*100)
print(r'\begin{table}[htbp]')
print(r'\centering')
print(r'\caption{Flash Balance model validation against experimental onset data (ceramics only)}')
print(r'\label{tab:flash_validation}')
print(r'\begin{tabular}{llccccc}')
print(r'\hline')
print(r'\textbf{Material} & \textbf{Family} & \textbf{E (V/cm)} & \textbf{$T_{pred}$ (°C)} & \textbf{$T_{exp}$ (°C)} & \textbf{Error (\%)} & \textbf{Ref.} \\')
print(r'\hline')

for family in families_order:
    family_results = [r for r in results_no_metals if r['family'] == family]
    if not family_results:
        continue

    for r in family_results:
        T_exp_C = r['T_experimental'] - 273
        mat_name = r['material'].replace('_', r'\_')

        if r['T_predicted'] is not None:
            T_pred_C = r['T_predicted'] - 273
            error_str = f"{r['error_percent']:+.1f}"
        else:
            T_pred_C = '--'
            error_str = '--'

        T_pred_str = f'{T_pred_C:.0f}' if isinstance(T_pred_C, (int, float)) else T_pred_C
        print(f"{mat_name} & {r['family']} & {r['E_field_Vcm']:.0f} & {T_pred_str} & {T_exp_C:.0f} & {error_str} & [{r['ref_num']}] \\\\")

print(r'\hline')
print(r'\end{tabular}')
print(r'\end{table}')

# CSV
print('\n\nCSV Format:')
print('-'*100)
print('Material,Family,E (V/cm),T_pred (K),T_pred (°C),T_exp (K),T_exp (°C),Error (%),Ref')

for r in results_no_metals:
    T_exp_C = r['T_experimental'] - 273
    if r['T_predicted'] is not None:
        T_pred_K = r['T_predicted']
        T_pred_C = T_pred_K - 273
        error = r['error_percent']
        print(f"{r['material']},{r['family']},{r['E_field_Vcm']},{T_pred_K:.0f},{T_pred_C:.0f},"
              f"{r['T_experimental']},{T_exp_C:.0f},{error:+.1f},[{r['ref_num']}]")
    else:
        print(f"{r['material']},{r['family']},{r['E_field_Vcm']},--,--,"
              f"{r['T_experimental']},{T_exp_C:.0f},--,[{r['ref_num']}]")
