# Flash Balance Data Extraction Guide

This guide explains how to extract Flash Balance parameters from literature papers for the Flash Balance Solver.

## Table of Contents

1. [Paper Screening](#1-paper-screening)
2. [Required Parameters](#2-required-parameters)
3. [Where to Find Parameters](#3-where-to-find-parameters)
4. [Graph Digitization](#4-graph-digitization)
5. [Unit Conversions](#5-unit-conversions)
6. [Confidence Levels](#6-confidence-levels)
7. [Common Pitfalls](#7-common-pitfalls)

---

## 1. Paper Screening

Before extracting data, categorize each paper using `paper_screening.csv`:

### Categories

| Category | Criteria | Action |
|----------|----------|--------|
| `FULL_DATA` | Has T_onset, E_field, AND conductivity (Arrhenius) data | Full extraction |
| `ONSET_ONLY` | Has T_onset and E_field, but no Arrhenius plot | Extract onset, estimate Ea from family |
| `CONDUCTIVITY_ONLY` | Has Arrhenius/conductivity data, no flash onset | Extract for reference only |
| `REVIEW` | Survey/review paper, no original experimental data | Skip, cite as reference |
| `THEORY` | Modeling/simulation, no experiments | Skip, note methodology |
| `DIFFERENT_TECHNIQUE` | Uses SPS, microwave, or other non-DC-field method | Flag technique, extract if comparable |
| `INSUFFICIENT` | Missing critical data | Document what's available |

### Screening Checklist

Answer these questions for each paper:

1. ☐ Does the paper report ORIGINAL flash sintering experiments?
2. ☐ Does it use DC electric field (not SPS/microwave)?
3. ☐ Is the onset temperature clearly reported or extractable?
4. ☐ Is the applied electric field clearly stated?
5. ☐ Are Arrhenius plots or conductivity data included?
6. ☐ Is current density at onset reported?
7. ☐ Are sample details provided (grain size, density)?

**Scoring:**
- Q1-4 all YES → Minimum `ONSET_ONLY`
- Q1-6 all YES → `FULL_DATA`
- Q1 NO → `REVIEW` or `THEORY`
- Q2 NO → `DIFFERENT_TECHNIQUE`

---

## 2. Required Parameters

### Critical (must have for solver validation)

| Parameter | Symbol | Unit | Description |
|-----------|--------|------|-------------|
| Onset Temperature | T_onset | K | Temperature at flash ignition |
| Electric Field | E | V/cm | Applied field at onset |

### Important (needed for accurate predictions)

| Parameter | Symbol | Unit | Description |
|-----------|--------|------|-------------|
| Activation Energy | Ea | eV | From Arrhenius plot slope |
| Pre-exponential | σ₀ | S/m | From Arrhenius plot intercept |
| Current Density | J_crit | A/cm² | Current at flash onset |

### Supplementary (can use defaults if missing)

| Parameter | Symbol | Unit | Source if Missing |
|-----------|--------|------|-------------------|
| Formation Enthalpy | ΔH | J/mol | NIST-JANAF tables |
| Formation Entropy | ΔS | J/(mol·K) | NIST-JANAF tables |
| Electrons per reaction | n | - | Stoichiometry (typically 4 for oxides) |

---

## 3. Where to Find Parameters

### 3.1 Onset Temperature (T_onset)

**Look in:**
- Abstract/Introduction: "Flash onset occurred at 850°C"
- Results section: Temperature vs current/power plots
- Tables: Experimental conditions summary
- Figure captions: "...showing flash onset at..."

**From graphs:** 
- Current density vs temperature: find sharp increase
- Power density vs temperature: find rapid rise
- Shrinkage vs temperature: find densification onset

**Watch out for:**
- Furnace temperature vs sample temperature (may differ!)
- Isothermal vs constant heating rate experiments
- Multiple onset conditions at different fields

### 3.2 Electric Field (E)

**Look in:**
- Experimental/Methods section: "Applied field of 100 V/cm"
- Tables: Experimental conditions
- Figure captions

**Calculate from:**
- Voltage and sample dimensions: E = V / L
- Ensure you have gauge length, not total sample length

**Common units:**
- V/cm (standard)
- V/m (divide by 100)
- kV/cm (multiply by 1000)

### 3.3 Activation Energy (Ea) and Pre-exponential (σ₀)

**Best source:** Arrhenius plot (ln(σ) or log(σ) vs 1/T or 1000/T)

**From the plot:**
```
ln(σ) = ln(σ₀) - Ea/(kB·T)

Slope = -Ea/kB
Intercept = ln(σ₀)

Therefore:
Ea (eV) = -slope × kB × 1000    [if x-axis is 1000/T]
Ea (eV) = -slope × kB           [if x-axis is 1/T]

σ₀ (S/m) = exp(intercept)       [if y-axis is ln(σ)]
σ₀ (S/m) = 10^intercept         [if y-axis is log₁₀(σ)]
```

Where kB = 8.617 × 10⁻⁵ eV/K

**Alternative sources:**
- Impedance spectroscopy data
- DC conductivity measurements at multiple temperatures
- Cited values from previous studies (verify!)

### 3.4 Current Density (J_crit)

**Look in:**
- Results section: "Current density reached X mA/mm² at onset"
- Current vs temperature plots (at onset point)
- Tables of experimental results

**Calculate from:**
- I and sample cross-section: J = I / A
- E and conductivity at onset: J = σ(T_onset) × E

### 3.5 Thermodynamic Data (ΔH, ΔS)

**Priority sources:**

1. **Paper-reported** (rare) - If authors measured/calculated
2. **NIST-JANAF Tables** - https://janaf.nist.gov/
3. **HSC Chemistry** - Software database
4. **FactSage/SGTE** - For complex compounds

**For oxides, use formation reaction:**
```
Example: ZrO2
Zr(s) + O2(g) → ZrO2(s)
ΔH°f = -1085 kJ/mol
ΔS°f = -178 J/(mol·K)
```

---

## 4. Graph Digitization

### Recommended Tool

**WebPlotDigitizer** (free): https://automeris.io/WebPlotDigitizer/

### Workflow

1. **Screenshot the graph** from PDF (high resolution)
2. **Open in WebPlotDigitizer**
3. **Calibrate axes:**
   - Click "2D (X-Y) Plot"
   - Define 4 calibration points (2 on X-axis, 2 on Y-axis)
   - Enter known values for each point
4. **Extract data points:**
   - Use "Automatic" mode for clean curves
   - Use "Manual" mode for noisy data
5. **Export as CSV**
6. **Save to** `data/digitized/[MaterialName]_[Author]/`

### Key Graph Types

#### Arrhenius Plot
- X-axis: 1000/T (K⁻¹) or 1/T (K⁻¹)
- Y-axis: ln(σ) or log₁₀(σ) (S/cm or S/m)
- Extract: All data points for linear fit
- Use: `scripts/analyze_arrhenius.py`

#### Onset Curve (J vs T)
- X-axis: Temperature (°C or K)
- Y-axis: Current density (mA/mm², A/cm²) or Current (mA)
- Extract: Points around the sharp rise
- Use: `scripts/find_onset.py`

#### I-V Characteristics
- X-axis: Voltage (V)
- Y-axis: Current (mA or A)
- Extract: Linear region slope at each temperature
- Calculate: σ = (I/V) × (L/A)

### Naming Convention

```
data/digitized/
├── BiFeO3_Raj/
│   ├── fig3_onset_JvsT.csv
│   ├── fig5_arrhenius.csv
│   └── extraction_notes.txt
├── 8YSZ_Francis/
│   └── ...
```

---

## 5. Unit Conversions

### Standard Units for Database

| Parameter | Standard Unit | 
|-----------|---------------|
| Temperature | K |
| Electric field | V/cm |
| Conductivity | S/m |
| Current density | A/cm² |
| Activation energy | eV |
| Enthalpy | J/mol |
| Entropy | J/(mol·K) |
| Grain size | μm |

### Common Conversions

```
Temperature:
  °C → K:  T(K) = T(°C) + 273.15

Electric field:
  V/m → V/cm:  E(V/cm) = E(V/m) / 100
  kV/cm → V/cm:  E(V/cm) = E(kV/cm) × 1000

Conductivity:
  S/cm → S/m:  σ(S/m) = σ(S/cm) × 100
  (Ω·cm)⁻¹ → S/m:  σ(S/m) = σ((Ω·cm)⁻¹) × 100
  (Ω·m)⁻¹ = S/m (same)

Current density:
  mA/mm² → A/cm²:  J(A/cm²) = J(mA/mm²) / 100
  A/m² → A/cm²:  J(A/cm²) = J(A/m²) / 10000

Activation energy:
  kJ/mol → eV:  Ea(eV) = Ea(kJ/mol) / 96.485
  J/mol → eV:  Ea(eV) = Ea(J/mol) / 96485
  
Enthalpy:
  kJ/mol → J/mol:  ΔH(J/mol) = ΔH(kJ/mol) × 1000
```

---

## 6. Confidence Levels

Assign a confidence code to each extracted value:

| Code | Level | Description |
|------|-------|-------------|
| `H` | High | Directly stated in table or text |
| `M` | Medium | Digitized from clear graph, R² > 0.98 |
| `L` | Low | Poor graph, estimated, or inferred |
| `D` | Default | From NIST-JANAF or standard tables |
| `C` | Calculated | Derived from other extracted values |

### Examples

```csv
T_onset_K,T_onset_source,T_onset_confidence
1123,Table 1,H
923,digitized_fig3,M
850,estimated_from_similar,L
```

---

## 7. Common Pitfalls

### Temperature Issues
- ❌ Using furnace temperature instead of sample temperature
- ❌ Confusing isothermal hold temperature with onset temperature
- ✅ Check if paper distinguishes furnace vs specimen temperature

### Field/Voltage Issues
- ❌ Using total voltage instead of field (V/cm)
- ❌ Ignoring electrode contact resistance
- ✅ Verify gauge length used for E = V/L calculation

### Conductivity Issues
- ❌ Mixing ionic and electronic conductivity
- ❌ Using AC impedance values without frequency context
- ✅ Note if conductivity is DC, AC, or total

### Arrhenius Plot Issues
- ❌ Wrong sign on slope (should be negative)
- ❌ Confusing ln and log₁₀ scales
- ❌ Fitting non-linear region
- ✅ Check R² > 0.95 for good fit
- ✅ Only fit linear Arrhenius region

### Unit Issues
- ❌ Forgetting unit conversions
- ❌ Assuming S/cm when paper uses S/m
- ✅ Always check axis labels and text descriptions
- ✅ Use `scripts/unit_converter.py` to validate

### Composite Materials
- ❌ Treating composite as single-phase
- ✅ Set `is_composite = TRUE`
- ✅ List all constituents and compositions

---

## Quick Reference Card

### Minimum Data for Solver

```
REQUIRED:
□ T_onset (K) - Flash onset temperature
□ E_field (V/cm) - Applied electric field

HIGHLY RECOMMENDED:
□ Ea (eV) - Activation energy
□ σ₀ (S/m) - Pre-exponential conductivity

CAN USE DEFAULTS:
□ ΔH (J/mol) - Use NIST-JANAF
□ ΔS (J/(mol·K)) - Use NIST-JANAF
□ n_electrons - Use stoichiometry (4 for oxides per O₂)
```

### Completeness Score Calculation

```
Score = sum of available parameters:
  T_onset:    25 points (critical)
  E_field:    25 points (critical)
  Ea:         15 points
  σ₀:         15 points
  ΔH:         10 points
  ΔS:         10 points
  
Total:       100 points
Ready for solver: Score ≥ 50 with T_onset and E_field present
```
