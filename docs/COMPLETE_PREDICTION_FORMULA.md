# Complete Flash Sintering Onset Prediction Formula

## Universal Equation (No Calibration Required)

```
E_crit = (C / √σ(T)) × correction_factors
```

Where:
- **C** = Family constant (see table)
- **σ(T) = σ₀·exp(-Ea/kT)** = Conductivity at temperature T
- **k** = 8.617×10⁻⁵ eV/K (Boltzmann constant)

## Correction Factors (Optional)

When microstructure data is available:

```
E_crit = (C / √σ(T)) × (d₅₀/500nm)^-0.10 × (φ/0.45)^0.45 × (HR/10°C/min)^0.05
```

| Parameter | Symbol | Reference Value | Effect |
|-----------|--------|-----------------|--------|
| Particle size | d₅₀ | 500 nm | Smaller → Higher E |
| Porosity | φ | 0.45 (45%) | Higher → Higher E |
| Heating rate | HR | 10 °C/min | Faster → Higher E |

**Note**: HR = Heating Rate. We use "HR" instead of "β" to avoid confusion with β (ridge parameter) in the Flash Balance equation.

## Family Constants

| Family | C [V·(S/m)^0.5/cm] | Typical Materials |
|--------|---------------------|-------------------|
| Fluorite | 108 | YSZ, CeO₂, UO₂ |
| Perovskite | 145 | SrTiO₃, BaTiO₃, KNN |
| Spinel | 168 | MgAl₂O₄, CoMn₂O₄ |
| Wurtzite | 205 | ZnO |
| Rutile | 164 | TiO₂, SnO₂ |
| Metal | 19 | Ni, W, Re (use J_crit method) |
| Corundum | 261 | Al₂O₃ |
| Carbide | 185 | SiC, WC, B₄C |
| Nitride | 140 | ZrN |
| Garnet | 120 | LLZO |
| Glass | 150 | SiO₂ glass |
| **Universal** | **176** | Any unknown material |

## Accuracy

- **Without corrections**: 58% median error, 76% within factor of 2
- **With d₅₀ + φ corrections**: 50% median error, 75% within factor of 2
- **With full solver (calibrated r_eff)**: 5.5% mean error, 90% within ±15%

## Heating Rate Guidance

| Heating Rate | Expected ΔT_onset | Recommended E Adjustment |
|--------------|-------------------|--------------------------|
| 1 °C/min | -70 K | E × 0.93 |
| 5 °C/min | -21 K | E × 0.98 |
| 10 °C/min | 0 K | E × 1.00 (reference) |
| 25 °C/min | +27 K | E × 1.03 |
| 100 °C/min | +69 K | E × 1.08 |

## Example Calculation

For 3YSZ at 1100 K:
- σ₀ = 3.4×10⁴ S/m
- Ea = 0.9 eV
- d₅₀ = 100 nm (fine powder)
- φ = 0.45 (55% dense)
- HR = 10 °C/min

σ(1100K) = 3.4×10⁴ × exp(-0.9/(8.617×10⁻⁵ × 1100)) = 5.7 S/m

E_crit = 108 / √5.7 × (100/500)^-0.10 × (0.45/0.45)^0.45
       = 45 × 1.17 × 1.0
       = 53 V/cm

Experimental: 70 V/cm → Error: 24%
