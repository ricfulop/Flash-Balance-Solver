# Flash Balance Solver

A Python implementation of the Flash Balance equation from the paper "Lattice resonant amplification in solids and Flash" by Ric Fulop (MIT, 2025).

## Overview

Flash is a distinct state of matter that emerges when an external electrical drive overwhelms intrinsic damping in a solid lattice. This solver implements the Flash Balance equation that predicts the onset conditions for this state across diverse material families.

### The Flash Balance Equation

```
ΔB = ksoft * |ΔG°(T)| - [Wph + Δμchem + nFEr]
```

Flash ignition occurs when **ΔB ≤ 0**

Where:
- **ksoft**: Lattice softening factor (reduces the thermochemical barrier)
- **ΔG°(T)**: Temperature-dependent Gibbs free energy of formation
- **Wph**: Phonon pumping work (non-thermal energy injection)
- **Δμchem**: Chemical potential work (e.g., from atmosphere)
- **nFEr**: Localized electrical work (n=electrons, F=Faraday, E=field, r=localization length)

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

## Supported Materials

| Family | Materials | Typical Onset | ksoft |
|--------|-----------|---------------|-------|
| Fluorite | 8YSZ, 3YSZ, GDC10 | 850-950°C @ 50-150 V/cm | ~0.10 |
| Rutile | TiO2, SnO2 | 600-700°C @ 100-200 V/cm | ~0.24 |
| Perovskite | SrTiO3, BaTiO3 | 500-600°C @ 50-100 V/cm | ~0.20 |
| Spinel | Al2O3 | 1100-1200°C @ 150-250 V/cm | ~0.36 |
| Metals | Cu, Ni, W | Joule-heating dominated | ~1.0 |
| Nitrides | TiN, Si3N4 | Variable | ~0.79 |
| Carbides | SiC | >1200°C | ~0.35 |
| Hydrides | ZrH2 | Low fields | Variable |

## Validation Results

The solver has been validated against experimental Flash sintering data from literature:

```
============================================================
SUMMARY
============================================================
Total samples validated:     16
Average absolute error:      15.6%
Samples within ±15%:         9/16 (56%)

Error by Material Family:
----------------------------------------
  fluorite               16.4% (n=7)
  nitride                 5.5% (n=1)
  perovskite             15.6% (n=3)
  rutile                 16.9% (n=3)
  spinel                 16.0% (n=2)
```

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

### 4. Temperature Sweep

```python
solver = FlashBalanceSolver(MATERIAL_DATABASE["8YSZ"])
E = 100 * 100  # 100 V/cm

for T in [800, 900, 1000, 1100, 1200]:
    delta_B = solver.flash_balance(T, E)
    status = "FLASH" if delta_B <= 0 else "Normal"
    print(f"T={T}K: ΔB={delta_B/1000:+.1f} kJ/mol - {status}")
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
- **beta**: Ridge parameter for ksoft calculation
- **alpha_res**: Resonant coupling efficiency
- **gamma**: Damping exponent for phonon pumping
- **delta_H, delta_S**: Formation enthalpy and entropy (J/mol, J/mol·K)
- **n_electrons**: Electrons transferred per reaction unit
- **r_eff**: Effective localization length (m)

### Universal Constants
- **q*/qD ≈ 0.73**: Universal phonon damping resonance
- **ksoft = 1 - β(q*/qD)²**: Lattice softening factor

## References

1. Fulop, R. "Lattice resonant amplification in solids and Flash" (2025)
2. Cologna, M., Rashkova, B. & Raj, R. Flash sintering of nanograin zirconia. J. Am. Ceram. Soc. 93, 3556-3559 (2010)
3. Francis, J.S.C. & Raj, R. Influence of the field and current limit on flash sintering. J. Am. Ceram. Soc. 96, 2754-2758 (2013)
4. Ding, G. et al. Unified theory of phonon in solids. Nature Physics 21, 1911-1919 (2025)

## License

MIT License - see LICENSE file for details.
