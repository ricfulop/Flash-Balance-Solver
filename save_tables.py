#!/usr/bin/env python3
"""
Save validation tables to files (excluding metals) with full parameters and DOIs.
"""

from flash_balance_solver import validate_solver, ORIGINAL_REFERENCES
import numpy as np

# Get validation results and filter out metals
results = validate_solver()
results_no_metals = [r for r in results if r['family'] != 'metal']

# Write full text table with parameters
with open('validation_table_no_metals.txt', 'w') as f:
    f.write('\n' + '='*175 + '\n')
    f.write('FLASH BALANCE SOLVER - ACCURACY VALIDATION TABLE (NO METALS)\n')
    f.write('='*175 + '\n')
    f.write(f"{'Material':<8} {'Family':<10} {'Ea':<5} {'beta':<5} {'alpha':<6} {'gamma':<6} "
            f"{'r_eff':<8} {'k_soft':<6} {'E':<6} {'T_exp':<6} {'T_pred':<7} {'Err%':<7} "
            f"{'Status':<6} {'Reference':<20} {'DOI':<30}\n")
    f.write(f"{'':8} {'':10} {'(eV)':<5} {'':5} {'':6} {'':6} "
            f"{'(um)':<8} {'':6} {'V/cm':<6} {'(K)':<6} {'(K)':<7} {'':7} "
            f"{'':6} {'':20} {'':30}\n")
    f.write('-'*175 + '\n')

    total_error = 0
    count = 0
    accurate_count = 0
    family_errors = {}

    for r in results_no_metals:
        if r['T_predicted'] is not None:
            status = 'PASS' if r['accuracy_met'] else 'FAIL'
            symbol = '✓' if r['accuracy_met'] else '✗'
            
            # Format r_eff
            r_eff_str = f"{r['r_eff_um']:.1f}" if r['r_eff_um'] < 100 else f"{r['r_eff_um']:.0f}"
            
            # Truncate DOI for display
            doi_short = r['doi'][:28] + '..' if len(r['doi']) > 30 else r['doi']
            
            f.write(f"{r['material']:<8} {r['family']:<10} {r['Ea']:<5.2f} {r['beta']:<5.2f} "
                    f"{r['alpha_res']:<6.2f} {r['gamma']:<6.1f} {r_eff_str:<8} {r['k_soft']:<6.2f} "
                    f"{r['E_field_Vcm']:<6.0f} {r['T_experimental']:<6.0f} {r['T_predicted']:<7.0f} "
                    f"{r['error_percent']:<+7.1f} {symbol}{status:<5} {r['reference']:<20} {doi_short:<30}\n")

            total_error += abs(r['error_percent'])
            count += 1
            if r['accuracy_met']:
                accurate_count += 1

            fam = r['family']
            if fam not in family_errors:
                family_errors[fam] = []
            family_errors[fam].append(abs(r['error_percent']))
        else:
            f.write(f"{r['material']:<8} {r['family']:<10} {r['Ea']:<5.2f} {r['beta']:<5.2f} "
                    f"{r['alpha_res']:<6.2f} {r['gamma']:<6.1f} {'N/A':<8} {'N/A':<6} "
                    f"{r['E_field_Vcm']:<6.0f} {r['T_experimental']:<6.0f} {'N/A':<7} "
                    f"{'N/A':<7} ✗FAIL  {r['reference']:<20}\n")

    f.write('-'*175 + '\n')
    f.write(f"\n{'='*60}\n")
    f.write('SUMMARY\n')
    f.write(f"{'='*60}\n")
    if count > 0:
        avg_error = total_error / count
        accuracy_rate = 100 * accurate_count / count
        f.write(f'Total samples validated:     {count}\n')
        f.write(f'Average absolute error:      {avg_error:.1f}%\n')
        f.write(f'Samples within ±15%:         {accurate_count}/{count} ({accuracy_rate:.0f}%)\n')
        f.write(f"\n{'Error by Material Family:':<30}\n")
        f.write('-'*40 + '\n')
        for fam, errors in sorted(family_errors.items()):
            avg = np.mean(errors)
            f.write(f'  {fam:<20} {avg:>6.1f}% (n={len(errors)})\n')
    
    # References section
    f.write(f"\n{'='*80}\n")
    f.write('REFERENCES\n')
    f.write(f"{'='*80}\n")
    seen_refs = set()
    for r in results_no_metals:
        if r['ref_num'] not in seen_refs:
            seen_refs.add(r['ref_num'])
            ref_info = ORIGINAL_REFERENCES.get(r['ref_num'], {})
            f.write(f"  [{r['ref_num']:>2}] {ref_info.get('author', 'Unknown')} ({ref_info.get('year', '')})\n")
            f.write(f"       DOI: {ref_info.get('doi', 'N/A')}\n")

# Write LaTeX table with parameters
with open('validation_table_latex.tex', 'w') as f:
    f.write(r'\begin{table}[htbp]' + '\n')
    f.write(r'\centering' + '\n')
    f.write(r'\caption{Flash Balance model validation against experimental onset data (ceramics only)}' + '\n')
    f.write(r'\label{tab:flash_validation}' + '\n')
    f.write(r'\begin{tabular}{llcccccccccc}' + '\n')
    f.write(r'\hline' + '\n')
    f.write(r'\textbf{Material} & \textbf{Family} & \textbf{$E_a$} & \textbf{$\beta$} & '
            r'\textbf{$\alpha_{res}$} & \textbf{$r_{eff}$} & \textbf{E} & '
            r'\textbf{$T_{exp}$} & \textbf{$T_{pred}$} & \textbf{Error} & \textbf{Ref.} \\' + '\n')
    f.write(r' & & (eV) & & & ($\mu$m) & (V/cm) & (°C) & (°C) & (\%) & \\' + '\n')
    f.write(r'\hline' + '\n')

    families_order = ['fluorite', 'rutile', 'perovskite', 'spinel', 'nitride', 'carbide']
    for family in families_order:
        family_results = [r for r in results_no_metals if r['family'] == family]
        if not family_results:
            continue

        for r in family_results:
            T_exp_C = r['T_experimental'] - 273
            mat_name = r['material'].replace('_', r'\_')
            r_eff_str = f"{r['r_eff_um']:.1f}"

            if r['T_predicted'] is not None:
                T_pred_C = r['T_predicted'] - 273
                error_str = f"{r['error_percent']:+.1f}"
                T_pred_str = f'{T_pred_C:.0f}'
            else:
                T_pred_str = '--'
                error_str = '--'

            f.write(f"{mat_name} & {r['family']} & {r['Ea']:.2f} & {r['beta']:.2f} & "
                    f"{r['alpha_res']:.2f} & {r_eff_str} & {r['E_field_Vcm']:.0f} & "
                    f"{T_exp_C:.0f} & {T_pred_str} & {error_str} & [{r['ref_num']}] \\\\\n")

    f.write(r'\hline' + '\n')
    f.write(r'\end{tabular}' + '\n')
    f.write(r'\end{table}' + '\n')

# Write CSV table with ALL columns including DOI
with open('validation_table.csv', 'w') as f:
    # Header with all columns
    f.write('Material,Family,Ea_eV,beta,alpha_res,gamma,r_eff_um,k_soft,'
            'E_Vcm,T_exp_K,T_exp_C,T_pred_K,T_pred_C,Error_pct,Status,Reference,DOI\n')
    
    for r in results_no_metals:
        T_exp_C = r['T_experimental'] - 273
        status = 'PASS' if r['accuracy_met'] else 'FAIL'
        
        if r['T_predicted'] is not None:
            T_pred_K = r['T_predicted']
            T_pred_C = T_pred_K - 273
            error = r['error_percent']
            f.write(f"{r['material']},{r['family']},{r['Ea']:.3f},{r['beta']:.2f},"
                    f"{r['alpha_res']:.2f},{r['gamma']:.1f},{r['r_eff_um']:.2f},{r['k_soft']:.3f},"
                    f"{r['E_field_Vcm']:.0f},{r['T_experimental']:.0f},{T_exp_C:.0f},"
                    f"{T_pred_K:.0f},{T_pred_C:.0f},{error:+.1f},{status},"
                    f"\"{r['reference']}\",{r['doi']}\n")
        else:
            f.write(f"{r['material']},{r['family']},{r['Ea']:.3f},{r['beta']:.2f},"
                    f"{r['alpha_res']:.2f},{r['gamma']:.1f},{r['r_eff_um']:.2f},,"
                    f"{r['E_field_Vcm']:.0f},{r['T_experimental']:.0f},{T_exp_C:.0f},"
                    f",,,,\"{r['reference']}\",{r['doi']}\n")

print("Files created:")
print("  - validation_table_no_metals.txt (full text table with parameters)")
print("  - validation_table_latex.tex (LaTeX table with parameters)")
print("  - validation_table.csv (CSV format with all columns including DOI)")
