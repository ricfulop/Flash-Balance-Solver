# Flash Balance Solver

A Python implementation of the Flash Balance equation from the paper "Lattice resonant amplification in solids and Flash" by Ric Fulop (MIT, 2025).

## Overview

Flash is a distinct state of matter that emerges when an external electrical drive overwhelms intrinsic damping in a solid lattice. This solver implements the Flash Balance equation that predicts the onset conditions for this state across diverse material families.

### The Flash Balance Equation

```
ΔB = ksoft × |ΔG°(T)| - [Wph + Δμchem + nFEr]
```

Flash ignition occurs when **ΔB ≤ 0**

Where:
- **ksoft**: Lattice softening factor (from first-principles β parameter)
- **ΔG°(T)**: Temperature-dependent Gibbs free energy of formation
- **Wph**: Phonon pumping work (non-thermal energy injection)
- **Δμchem**: Chemical potential work (e.g., from atmosphere)
- **nFEr**: Localized electrical work (n=electrons, F=Faraday, E=field, r=localization length)

## v1.2 - Complete Recalibration (January 2026)

### Key Updates
- **121 datapoints** fully calibrated across 17 material families
- **First-principles β derivation**: `β = 6400 × Z*/T_D² + 1.30` for ionic ceramics
- **0% prediction error** for 117/121 ceramic datapoints (analytical calibration)
- **J_crit model** for metals validated (Ni: +2.3% error)
- Complete DOI references for all experimental data

### Calibration Statistics

| Metric | Value |
|--------|-------|
| Total datapoints | 121 |
| Successfully calibrated | 117 (96.7%) |
| r_eff range | 3.3 - 825.6 μm |
| r_eff median | 32.0 μm |
| Prediction error | 0% (analytical) |

## Installation

```bash
pip install numpy scipy
```

Or using the requirements file:
```bash
pip install -r requirements.txt
```

## Quick Start

```python
from flash_balance_solver import FlashBalanceSolver, MATERIAL_DATABASE

# Create solver for 8YSZ (Yttria-Stabilized Zirconia)
material = MATERIAL_DATABASE["8YSZ"]
solver = FlashBalanceSolver(material)

# Predict onset temperature at 100 V/cm
E_field = 100 * 100  # Convert V/cm to V/m
T_onset = solver.solve_onset_temperature(E_field)
print(f"Predicted onset: {T_onset:.0f} K ({T_onset-273:.0f} °C)")

# Solve for critical field at 850°C
T = 1123  # K
E_crit = solver.solve_critical_field(T)
print(f"Critical field: {E_crit/100:.0f} V/cm")
```

## Supported Material Families

| Family | Count | β | k_soft | r_eff median | Typical Onset |
|--------|-------|---|--------|--------------|---------------|
| Fluorite | 24 | 1.38 | 0.26 | 39 μm | 850-1200°C |
| Perovskite | 25 | 1.50 | 0.20 | 22 μm | 800-1200°C |
| Spinel | 16 | 1.34 | 0.29 | 19 μm | 400-1900°C |
| Carbide | 11 | 0.31 | 0.83 | 184 μm | 1100-2400°C |
| Corundum | 7 | 1.32 | 0.30 | 5 μm | 1200-1500°C |
| Wurtzite | 7 | 1.38 | 0.27 | 21 μm | 800-1000°C |
| Garnet | 5 | 1.41 | 0.25 | 32 μm | 900-1200°C |
| Rutile | 4 | 1.38 | 0.26 | 67 μm | 900-1400°C |
| Glass-ceramic | 4 | 1.49 | 0.21 | 5 μm | 1200-1400°C |
| Metal | 3 | 0.50 | 0.73 | 102-500 μm | J_crit controlled |
| Nitride | 2 | 0.34 | 0.82 | 330 μm | 1000-1100°C |
| Oxide | 6 | 1.40 | 0.25 | 34 μm | 800-1400°C |

## First-Principles β Derivation

The ridge parameter β is derived from material properties:

### For Ionic Ceramics
```
β = 6400 × Z*/T_D² + 1.30
```

Where:
- **Z***: Born effective charge (from DFT or literature)
- **T_D**: Debye temperature (K)

### For Other Material Classes
| Class | β Formula | Notes |
|-------|-----------|-------|
| Covalent (SiC, B4C) | 0.30-0.31 | Strong bonds, minimal softening |
| Nitrides | 0.33-0.35 | Mixed ionic-covalent |
| Metals | 0.50 | Frenkel pair nucleation physics |

### k_soft Calculation
```
k_soft = 1 - β × (q*/q_D)² = 1 - 0.533 × β
```

Where q*/q_D ≈ 0.73 is the universal phonon ridge value from Nature Physics (2025).

## Data Files

### Primary Dataset
- **`data/calibrated_dataset_with_doi.csv`** - 121 fully calibrated datapoints with:
  - Material name and family
  - β and k_soft (first-principles)
  - Ea, σ₀, r_eff (calibrated)
  - T_onset, E-field, J_crit
  - DOI references

### Supporting Data
- `data/complete_validation_table.csv` - Full experimental data
- `data/derived_beta_values.csv` - First-principles β calculations
- `data/collapse_analysis/` - Universality proof figures

## Usage Examples

### 1. Predict Flash Onset Temperature

```python
from flash_balance_solver import FlashBalanceSolver, MATERIAL_DATABASE

# For TiO2 at 150 V/cm
solver = FlashBalanceSolver(MATERIAL_DATABASE["TiO2"])
T_onset = solver.solve_onset_temperature(150 * 100)  # V/cm to V/m
print(f"TiO2 onset at 150 V/cm: {T_onset:.0f} K")
```

### 2. Calculate Critical Current Density

```python
# Critical J for SrTiO3 at 550°C
solver = FlashBalanceSolver(MATERIAL_DATABASE["SrTiO3"])
T = 823  # K (550°C)
result = solver.solve_critical_current_density(T)
if result:
    J_crit, E_crit = result
    print(f"Critical J: {J_crit:.2e} A/m²")
    print(f"Critical E: {E_crit/100:.0f} V/cm")
```

### 3. Analyze Flash Balance Components

```python
solver = FlashBalanceSolver(MATERIAL_DATABASE["8YSZ"])
components = solver.get_balance_components(T=1100, E=100*100)

print(f"Softened barrier: {components['barrier']/1000:.1f} kJ/mol")
print(f"Phonon work:      {components['W_phonon']/1000:.1f} kJ/mol")
print(f"Electrical work:  {components['W_electrical']/1000:.1f} kJ/mol")
print(f"Flash Balance ΔB: {components['delta_B']/1000:.1f} kJ/mol")
print(f"Flash accessible: {components['flash_accessible']}")
```

### 4. Quick E_crit Estimation (Without Calibrated r_eff)

```python
from flash_balance_solver import estimate_E_crit

# Estimate for a new fluorite material
result = estimate_E_crit(
    Ea=0.9,           # Activation energy (eV)
    sigma_0=3.4e4,    # Pre-exponential (S/m)
    T=1100,           # Temperature (K)
    family="Fluorite",
    d50_nm=100        # Optional: particle size
)
print(f"E_crit ≈ {result['E_crit']:.0f} V/cm")
print(f"Confidence: {result['confidence']}")
```

### 5. Metal J_crit Estimation

```python
from flash_balance_solver import estimate_J_crit_metal

# For Ni at 1000°C
result = estimate_J_crit_metal(T_target=1273, T_debye=450)
print(f"J_crit ≈ {result['J_crit']:.1f} A/mm²")
# Output: J_crit ≈ 20.5 A/mm² (matches experimental 20 A/mm²)
```

## Run Validation

```bash
python flash_balance_solver.py
```

This will run example calculations and print the full validation table.

## Key Parameters

### Material Parameters
- **Ea**: Activation energy for electrical conduction (eV)
- **sigma_0**: Pre-exponential conductivity factor (S/m)
- **beta (β)**: Ridge parameter (from first-principles or calibration)
- **k_soft**: Lattice softening factor = 1 - 0.533β
- **alpha_res**: Resonant coupling efficiency
- **gamma**: Damping exponent for phonon pumping
- **delta_H, delta_S**: Formation enthalpy and entropy (J/mol, J/mol·K)
- **n_electrons**: Electrons transferred per reaction unit
- **r_eff**: Effective localization length (m) - calibrated per material

### Universal Constants
- **q*/qD ≈ 0.73**: Universal phonon damping resonance (Nature Physics 2025)
- **ksoft = 1 - 0.533β**: Lattice softening factor

## References

1. Fulop, R. "Lattice resonant amplification in solids and Flash" (2025)
2. Ding, G. et al. "Unified theory of phonon in solids" Nature Physics 21, 1911-1919 (2025)
3. Cologna, M., Rashkova, B. & Raj, R. "Flash sintering of nanograin zirconia" J. Am. Ceram. Soc. 93, 3556-3559 (2010)
4. Francis, J.S.C. & Raj, R. "Influence of the field and current limit on flash sintering" J. Am. Ceram. Soc. 96, 2754-2758 (2013)

## Version History

| Version | Date | Changes |
|---------|------|---------|
| v1.2 | Jan 2026 | Complete r_eff recalibration, first-principles β, 121 datapoints |
| v1.1 | - | λ_flash interpolation system |
| v1.0 | - | Initial release |

## License

MIT License - see LICENSE file for details.
