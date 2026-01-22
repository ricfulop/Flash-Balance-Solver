# Flash Balance First Principles Parameter Extraction Workflow

## Supplementary Material for Flash Balance Solver

This document describes the workflow for extracting Flash Balance model parameters from first principles for new materials. This enables researchers to apply the Flash Balance equation to materials not included in the original validation database.

---

## Table of Contents

1. [Overview](#overview)
2. [Required Input Data](#required-input-data)
3. [Parameter Extraction Procedure](#parameter-extraction-procedure)
4. [Data Sources](#data-sources)
5. [Parameter Guidelines by Material Family](#parameter-guidelines-by-material-family)
6. [Validation Criteria](#validation-criteria)
7. [Worked Example: BiFeO3](#worked-example-bifeo3)
8. [Troubleshooting](#troubleshooting)
9. [References](#references)

---

## Overview

The Flash Balance equation predicts Flash sintering onset:

```
ΔB = k_soft × |ΔG°(T)| - [W_ph + Δμ_chem + nFEr_eff]
```

Flash occurs when ΔB ≤ 0.

### Required Parameters

| Parameter | Symbol | Units | Description |
|-----------|--------|-------|-------------|
| Activation energy | Ea | eV | Arrhenius conductivity activation |
| Pre-exponential | σ₀ | S/m | Arrhenius pre-exponential |
| Ridge parameter | β | - | Phonon softening parameter |
| Resonant coupling | α_res | - | Electric field-phonon coupling |
| Damping exponent | γ | - | Temperature damping factor |
| Formation enthalpy | ΔH | J/mol | Standard formation enthalpy |
| Formation entropy | ΔS | J/(mol·K) | Standard formation entropy |
| Electron count | n | - | Electrons per reaction unit |
| Localization length | r_eff | m | Effective energy localization |

---

## Required Input Data

### Minimum Required Data

1. **Material identification**
   - Chemical formula
   - Crystal structure / material family

2. **One experimental Flash onset point**
   - Onset temperature T_onset (K)
   - Applied electric field E (V/cm)

3. **Thermodynamic data** (at least one source)
   - Formation enthalpy ΔH_f
   - Formation entropy ΔS_f (or can be estimated)

4. **Electrical conductivity data**
   - Arrhenius activation energy Ea
   - Pre-exponential factor σ₀

### Optional Data (improves accuracy)

- Multiple experimental onset points at different fields
- Temperature-dependent conductivity measurements
- Phonon spectroscopy data (for β refinement)

---

## Parameter Extraction Procedure

### Step 1: Material Classification

Identify the material family based on crystal structure:

| Family | Structure | Examples |
|--------|-----------|----------|
| Fluorite | Cubic fluorite | YSZ, CeO2, GDC |
| Rutile | Tetragonal rutile | TiO2, SnO2 |
| Perovskite | ABO3 perovskite | BaTiO3, SrTiO3, BiFeO3 |
| Spinel | AB2O4 spinel | MgAl2O4, NiFe2O4 |
| Wurtzite | Hexagonal wurtzite | ZnO, AlN |
| Carbide | Various | SiC, WC, TiC |
| Nitride | Various | Si3N4, TiN |

### Step 2: Thermodynamic Data Extraction

**Formation Enthalpy (ΔH_f)**

1. Look up in thermodynamic databases:
   - NIST-JANAF Thermochemical Tables (preferred)
   - HSC Chemistry Database
   - FactSage / SGTE databases

2. Identify the formation reaction:
   - For oxides: Metal(s) + O2(g) → Oxide(s)
   - For perovskites: AO + BO2 → ABO3
   - For carbides: Metal(s) + C(s) → Carbide(s)

3. Use standard state values (298 K, 1 atm)

4. Normalize per reaction unit:
   - Oxides: per O2 (n=4 electrons)
   - Nitrides: per N2 (n=6 electrons)
   - Carbides: per formula unit

**Formation Entropy (ΔS_f)**

If not available, estimate from:

```
ΔS_f ≈ -k × |ΔH_f|
```

where k depends on material family:

| Family | k (J/(mol·K) per kJ/mol) |
|--------|--------------------------|
| Fluorite | 0.16 |
| Rutile | 0.20 |
| Perovskite | 0.12 |
| Spinel | 0.09 |
| Carbide | 0.25 |

### Step 3: Electrical Properties

**From Arrhenius Plot:**

Fit conductivity data to:

```
σ(T) = σ₀ × exp(-Ea / kT)
```

Where:
- Ea = slope × k_B (in eV)
- σ₀ = y-intercept (in S/m)

**Data Sources:**
- Published conductivity measurements
- Impedance spectroscopy data
- Four-point probe measurements

**Quality Checks:**
- Linear Arrhenius behavior over measurement range
- Temperature range includes expected onset region
- Consistent with conduction mechanism (ionic/electronic/mixed)

### Step 4: Family-Based Parameter Estimation

Use typical values for material family:

#### Fluorite Oxides
```
β = 1.69 ± 0.10
α_res = 0.15 ± 0.03
γ = 2.0 ± 0.2
k_soft ≈ 0.10
```

#### Rutile Oxides
```
β = 1.43 ± 0.08
α_res = 0.22 ± 0.05
γ = 1.8 ± 0.2
k_soft ≈ 0.24
```

#### Perovskites
```
β = 1.50 ± 0.10
α_res = 0.28 ± 0.05 (0.30-0.40 for multiferroics)
γ = 1.4 ± 0.2
k_soft ≈ 0.20
```

#### Spinel/Corundum
```
β = 1.20 ± 0.15
α_res = 0.20 ± 0.04
γ = 1.8 ± 0.2
k_soft ≈ 0.36
```

### Step 5: r_eff Calibration

The effective localization length r_eff is calibrated using one experimental onset point:

1. Set initial r_eff to family typical value
2. Use numerical optimization to find r_eff that gives:
   ```
   T_predicted(E_exp) = T_onset_exp
   ```
3. Check that calibrated r_eff is within physical range

**Typical r_eff Ranges:**

| Family | r_eff (μm) |
|--------|------------|
| Fluorite | 10 - 30 |
| Rutile | 12 - 25 |
| Perovskite | 15 - 55 |
| Spinel | 20 - 120 |
| Carbide | 4 - 15 |

### Step 6: Validation

Test with additional experimental data:

1. Predict onset at different electric fields
2. Compare with experimental values
3. Target: < 15% error

---

## Data Sources

### Recommended Thermodynamic Databases

| Source | URL | Coverage |
|--------|-----|----------|
| NIST-JANAF | https://janaf.nist.gov/ | Common oxides, carbides |
| HSC Chemistry | https://www.hsc-chemistry.com/ | Comprehensive inorganics |
| FactSage | https://www.factsage.com/ | Complex systems |
| SGTE | https://www.sgte.net/ | Pure substances |

### Electrical Property Sources

- Published literature (check DOI)
- Materials property databases (MatWeb, etc.)
- Your own measurements (preferred)

---

## Parameter Guidelines by Material Family

### Fluorite Oxides (YSZ, CeO2, GDC)

**Typical Parameters:**
- Ea: 0.7 - 1.1 eV
- σ₀: 10⁴ - 10⁵ S/m
- β: 1.69
- α_res: 0.15
- γ: 2.0
- n_electrons: 4 (per O2)

**Notes:**
- Ionic conduction dominates
- Well-characterized family
- Most reliable predictions

### Perovskites (BaTiO3, SrTiO3, BiFeO3)

**Typical Parameters:**
- Ea: 0.3 - 0.8 eV
- σ₀: 5×10² - 10⁴ S/m
- β: 1.50
- α_res: 0.28 (0.35 for multiferroics)
- γ: 1.4
- n_electrons: 4 (per O2)

**Notes:**
- Large r_eff due to polar nanodomains
- Multiferroics have enhanced coupling
- Mixed ionic/electronic conduction common

### Carbides (SiC, WC)

**Typical Parameters:**
- Ea: 0.3 - 1.0 eV (varies widely)
- σ₀: 10⁴ - 10⁶ S/m
- β: 1.20
- α_res: 0.05 - 0.10
- γ: 1.8
- n_electrons: 4

**Notes:**
- WC is metallic - may not follow model
- SiC: use k_soft = 0.88 override
- Low phonon coupling (covalent bonds)

---

## Validation Criteria

### Acceptable Accuracy

| Category | Error Threshold |
|----------|-----------------|
| Excellent | < 5% |
| Good | 5 - 10% |
| Acceptable | 10 - 15% |
| Poor | > 15% |

### Warning Signs

1. **r_eff hits calibration bounds**
   - Too small (< 1 μm): Other parameters need adjustment
   - Too large (> 200 μm): Check thermodynamic data

2. **Large errors at different fields**
   - Model may need field-dependent parameters
   - Check if material has phase transitions

3. **Negative T_predicted**
   - Thermodynamic data may be incorrect
   - Check sign of ΔH

---

## Worked Example: BiFeO3

### Step 1: Material Identification
- Formula: BiFeO3
- Family: Perovskite
- Note: Multiferroic (ferroelectric + antiferromagnetic)

### Step 2: Thermodynamic Data
- Reaction: ½Bi2O3 + ½Fe2O3 → BiFeO3
- ΔH_f = -945 kJ/mol (NIST)
- ΔS_f ≈ -175 J/(mol·K) (estimated)

### Step 3: Electrical Properties
- Ea = 0.70 eV (higher due to p-type + ionic)
- σ₀ = 1.0×10⁴ S/m (mixed conduction)

### Step 4: Family Parameters
- β = 1.50 (perovskite typical)
- α_res = 0.35 (enhanced for multiferroic)
- γ = 1.4

### Step 5: Calibration
- Experimental: T_onset = 923 K at E = 100 V/cm
- Calibrated r_eff = 15.7 μm

### Step 6: Validation
- At 100 V/cm: Error = 0.0% ✓
- At 150 V/cm: Error = -6.0% ✓

### Final Parameters
```python
"BiFeO3": MaterialParameters(
    name="Bismuth Ferrite",
    family=MaterialFamily.PEROVSKITE,
    Ea=0.70,
    sigma_0=1.0e4,
    beta=1.50,
    alpha_res=0.35,
    gamma=1.4,
    delta_H=-945000,
    delta_S=-175,
    n_electrons=4,
    r_eff=1.57e-05,
),
```

---

## Troubleshooting

### Problem: r_eff hits lower bound

**Possible causes:**
1. Ea too low → Increase Ea
2. α_res too low → Increase for polar materials
3. |ΔH| too small → Verify thermodynamic data
4. Material may not follow standard model (metallic)

### Problem: r_eff hits upper bound

**Possible causes:**
1. Ea too high → Verify Arrhenius fit
2. |ΔH| too large → Check reaction stoichiometry
3. σ₀ too low → Check conductivity data

### Problem: Good at one field, poor at others

**Possible causes:**
1. Phase transition in temperature range
2. Field-dependent parameters needed
3. Different conduction mechanism at high field

---

## Using the Extraction Script

Run the interactive tool:

```bash
python scripts/extract_material_params.py
```

The script will guide you through each step and output:
1. MaterialParameters Python code
2. JSON file with all parameters and metadata
3. Validation results

---

## References

1. Flash Balance paper (this work)
2. NIST-JANAF Thermochemical Tables
3. Cologna et al. (2010) - Flash sintering of YSZ
4. Raj (2012) - Flash sintering mechanism
5. Materials-specific literature (see DOIs in output)

---

*This workflow is part of the Flash Balance Solver supplementary material.*
*For questions or contributions, see the project repository.*
