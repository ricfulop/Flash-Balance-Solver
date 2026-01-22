# Complete Flash Sintering Onset Prediction Formula

## Quick Reference Card

### Universal Equation (No Calibration Required)

```
E_crit = C / √σ(T)
```

Where:
- **E_crit** = Critical electric field (V/cm)
- **C** = Family constant (see table below)
- **σ(T) = σ₀·exp(-Ea/kT)** = Conductivity at temperature T (S/m)
- **k** = 8.617×10⁻⁵ eV/K (Boltzmann constant)

---

## Family Constants

| Family | C [V·(S/m)^0.5/cm] | Typical Materials |
|--------|---------------------|-------------------|
| Fluorite | 108 | YSZ, CeO₂, GDC, UO₂ |
| Perovskite | 145 | SrTiO₃, BaTiO₃, KNN |
| Spinel | 168 | MgAl₂O₄, MnCo₂O₄ |
| Wurtzite | 205 | ZnO |
| Rutile | 164 | TiO₂, SnO₂ |
| Corundum | 261 | α-Al₂O₃ |
| Carbide | 185 | SiC, WC, B₄C |
| Nitride | 140 | ZrN, (Nb,Ta,Ti)N |
| Garnet | 120 | LLZO |
| Glass | 150 | Porous SiO₂ |
| Metal | 19 | Ni, W, Re (use J_crit method instead) |
| **Universal** | **176** | Any unknown material |

---

## Correction Factors (Optional)

When microstructure data is available, apply corrections:

```
E_crit = (C / √σ(T)) × (d₅₀/500nm)^-0.10 × (φ/0.45)^0.45 × (HR/10)^0.05
```

| Parameter | Symbol | Reference | Effect on E_crit |
|-----------|--------|-----------|------------------|
| Particle size | d₅₀ | 500 nm | Smaller particles → Higher E |
| Porosity | φ | 0.45 (45%) | More porous → Higher E |
| Heating rate | HR | 10 °C/min | Faster heating → Higher E |

**Important:** 
- φ = 1 - ρ_rel (e.g., 55% dense → φ = 0.45)
- HR = Heating Rate in °C/min (NOT β, to avoid confusion with ridge parameter)

---

## Heating Rate Guidance

| HR (°C/min) | T_onset Shift | E Adjustment |
|-------------|---------------|--------------|
| 1 | -70 K | × 0.93 |
| 5 | -21 K | × 0.98 |
| **10** | **0 K (reference)** | **× 1.00** |
| 25 | +27 K | × 1.03 |
| 100 | +69 K | × 1.08 |

---

## Dense Material Adjustments

| Material State | λ_flash Factor | E_crit Factor |
|----------------|----------------|---------------|
| Green compact (50-65% dense) | × 1.0 | × 1.0 |
| Dense polycrystal (>95% dense) | × 0.3-0.5 | × 2 |
| Single crystal | × 0.2-0.5 | × 3-5 |

---

## Accuracy Summary

| Approach | Median Error | Within Factor of 2 |
|----------|--------------|-------------------|
| Universal formula only | 58% | 77% |
| With microstructure corrections | ~50% | 75% |
| With λ_flash interpolation | ~25% | 85% |
| **Full solver (calibrated λ_flash)** | **5%** | **99%** |

---

## λ_flash Interpolation (For New Materials)

If you don't have calibrated λ_flash, estimate it:

```
log₁₀(λ_flash) = a + b·Ea + c·log₁₀(σ₀)
```

| Family | a | b | c | Mean λ (µm) |
|--------|-----|------|------|-------------|
| Fluorite | 0.93 | -0.38 | +0.08 | 8.2 |
| Perovskite | 2.06 | -1.18 | -0.21 | 6.7 |
| Spinel | 1.00 | -1.88 | +0.36 | 14.5 |
| Carbide | 1.88 | +0.42 | -0.39 | 9.4 |
| Corundum | -0.69 | +0.65 | +0.18 | 1.7 |
| **Universal** | **0.50** | **-1.05** | **+0.30** | 10 |

---

## Example Calculation

### 3YSZ at 1100 K (827°C)

**Given:**
- σ₀ = 3.4×10⁴ S/m
- Ea = 0.9 eV
- d₅₀ = 100 nm (fine powder)
- φ = 0.45 (55% dense green compact)
- HR = 10 °C/min

**Step 1: Calculate conductivity**
```
σ(1100K) = 3.4×10⁴ × exp(-0.9 / (8.617×10⁻⁵ × 1100))
         = 3.4×10⁴ × exp(-9.49)
         = 5.7 S/m
```

**Step 2: Apply universal formula (Fluorite, C = 108)**
```
E_crit = 108 / √5.7 = 45 V/cm
```

**Step 3: Apply microstructure corrections**
```
E_crit = 45 × (100/500)^-0.10 × (0.45/0.45)^0.45 × (10/10)^0.05
       = 45 × 1.17 × 1.0 × 1.0
       = 53 V/cm
```

**Result:** E_crit ≈ 53 V/cm  
**Experimental:** 70 V/cm  
**Error:** 24%

---

## Metal Flash Sintering (Current-Rate)

For metals, use critical current density instead:

```
J_crit = 4.3 × (T/T_debye)^1.5  [A/mm²]
```

| Metal | T_debye (K) | Typical J_crit (A/mm²) |
|-------|-------------|----------------------|
| Ni | 450 | ~20 at 1000°C |
| W | 400 | ~31 at 1200°C |
| Re | 416 | ~10 at 900°C |

**Requirement:** Flash in metals requires T ≥ T_debye (Frenkel pair nucleation).

---

## Nomenclature Note

| Symbol | Meaning | Units |
|--------|---------|-------|
| E_crit | Critical electric field | V/cm |
| σ(T) | Conductivity at T | S/m |
| σ₀ | Pre-exponential conductivity | S/m |
| Ea | Activation energy | eV |
| λ_flash | Flash activation length | µm |
| d₅₀ | Median particle size | nm |
| φ | Porosity (= 1 - ρ_rel) | - |
| HR | Heating Rate | °C/min |
| β | Ridge parameter (phonon theory) | - |
| T_D | Debye temperature | K |
| J_crit | Critical current density | A/mm² |

**Important:** HR (Heating Rate) is used instead of β to avoid confusion with the ridge parameter β in the Flash Balance equation.

---

## References

See `data/complete_validation_table.csv` for the full validation dataset (135 datapoints, 49 publications).

See `data/lambda_flash_interpolation_map.png` for visual λ_flash estimation chart.
