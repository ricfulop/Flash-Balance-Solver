#!/usr/bin/env python3
"""
Clean Validation Data - Mark Unusable Entries

This script marks data points that should be excluded from analysis:
1. RT Flash (T_onset ≤ 350K) - different mechanism than thermal flash
2. Missing critical data (no J_calc, Ea=0)
3. Touch-free flash experiments (different setup)
4. Metal flash (Ni, W, Re - different conduction mechanism)
"""

import pandas as pd
from pathlib import Path

def clean_validation_data():
    script_dir = Path(__file__).parent
    data_dir = script_dir.parent / 'data'
    
    # Load data
    df = pd.read_csv(data_dir / 'complete_validation_table.csv')
    df.columns = df.columns.str.strip()
    
    if df.columns[0] == '#':
        df = df.rename(columns={'#': 'ID'})
    
    print(f"Total entries: {len(df)}")
    
    # Add exclusion columns
    df['Exclude'] = False
    df['Exclude_Reason'] = ''
    
    # 1. RT Flash - T_onset ≤ 350K (different mechanism)
    mask_rt = df['T_onset(K)'] <= 350
    df.loc[mask_rt, 'Exclude'] = True
    df.loc[mask_rt, 'Exclude_Reason'] = 'RT_Flash'
    print(f"\nRT Flash (T ≤ 350K): {mask_rt.sum()} entries")
    for _, row in df[mask_rt].iterrows():
        print(f"  - ID {row['ID']}: {row['Material']} at {row['T_onset(K)']}K")
    
    # 2. Missing J_calc data
    mask_no_j = df['J_calc(mA/mm²)'].isna() | (df['J_calc(mA/mm²)'] == 0)
    df.loc[mask_no_j & ~df['Exclude'], 'Exclude'] = True
    df.loc[mask_no_j & (df['Exclude_Reason'] == ''), 'Exclude_Reason'] = 'Missing_J'
    new_missing = mask_no_j & (df['Exclude_Reason'] == 'Missing_J')
    print(f"\nMissing J_calc: {new_missing.sum()} entries")
    for _, row in df[new_missing].iterrows():
        print(f"  - ID {row['ID']}: {row['Material']}")
    
    # 3. Ea = 0 (anomalous/different mechanism)
    mask_ea0 = df['Ea(eV)'] == 0
    df.loc[mask_ea0 & ~df['Exclude'], 'Exclude'] = True
    df.loc[mask_ea0 & (df['Exclude_Reason'] == ''), 'Exclude_Reason'] = 'Ea_Zero'
    new_ea0 = mask_ea0 & (df['Exclude_Reason'] == 'Ea_Zero')
    print(f"\nEa = 0: {new_ea0.sum()} entries")
    for _, row in df[new_ea0].iterrows():
        print(f"  - ID {row['ID']}: {row['Material']}")
    
    # 4. Metal flash (Ni, W, Re) - different conduction mechanism
    metals = ['Ni', 'W', 'Re']
    mask_metal = df['Material'].isin(metals)
    df.loc[mask_metal & ~df['Exclude'], 'Exclude'] = True
    df.loc[mask_metal & (df['Exclude_Reason'] == ''), 'Exclude_Reason'] = 'Metal_Flash'
    new_metal = mask_metal & (df['Exclude_Reason'] == 'Metal_Flash')
    print(f"\nMetal flash: {new_metal.sum()} entries")
    for _, row in df[new_metal].iterrows():
        print(f"  - ID {row['ID']}: {row['Material']}")
    
    # 5. Single crystals (different from powder/green compact flash)
    mask_sc = df['Material'].str.contains('Single Crystal', case=False, na=False)
    df.loc[mask_sc & ~df['Exclude'], 'Exclude'] = True
    df.loc[mask_sc & (df['Exclude_Reason'] == ''), 'Exclude_Reason'] = 'Single_Crystal'
    new_sc = mask_sc & (df['Exclude_Reason'] == 'Single_Crystal')
    print(f"\nSingle Crystal: {new_sc.sum()} entries")
    for _, row in df[new_sc].iterrows():
        print(f"  - ID {row['ID']}: {row['Material']}")
    
    # Summary
    total_excluded = df['Exclude'].sum()
    total_valid = len(df) - total_excluded
    
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total entries:    {len(df)}")
    print(f"Excluded:         {total_excluded}")
    print(f"Valid for use:    {total_valid}")
    
    # Breakdown by reason
    print("\nExclusion reasons:")
    for reason in df[df['Exclude']]['Exclude_Reason'].unique():
        count = (df['Exclude_Reason'] == reason).sum()
        print(f"  - {reason}: {count}")
    
    # Save cleaned data
    output_file = data_dir / 'validation_table_cleaned.csv'
    df.to_csv(output_file, index=False)
    print(f"\n✓ Saved: {output_file}")
    
    # Also save a version with only valid entries
    df_valid = df[~df['Exclude']].drop(columns=['Exclude', 'Exclude_Reason'])
    valid_file = data_dir / 'validation_table_valid_only.csv'
    df_valid.to_csv(valid_file, index=False)
    print(f"✓ Saved: {valid_file}")
    
    # Print materials in valid set
    print("\n" + "=" * 60)
    print("VALID MATERIALS (by count)")
    print("=" * 60)
    valid_counts = df_valid.groupby('Material').size().sort_values(ascending=False)
    for mat, count in valid_counts.items():
        print(f"  {mat}: {count}")

if __name__ == '__main__':
    clean_validation_data()
