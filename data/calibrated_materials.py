"""Calibrated MaterialParameters from Full Data Papers."""

from flash_balance_solver import MaterialParameters, MaterialFamily

CALIBRATED_MATERIALS = {
    "Al2O3-ZrO2": MaterialParameters(
        name="Alumina-Zirconia Composite",
        family=MaterialFamily.SPINEL,
        Ea=0.782,
        sigma_0=4.36e+06,
        beta=1.20,
        alpha_res=0.20,
        gamma=1.8,
        delta_H=-1676000,
        delta_S=-210,
        n_electrons=4,
        r_eff=2.00e-05,  # Calibrated at T=1173K, E=100V/cm
    ),  # DOI: 10.1016/j.jeurceramsoc.2015.02.033
    "BiFeO3": MaterialParameters(
        name="Bismuth Ferrite",
        family=MaterialFamily.PEROVSKITE,
        Ea=0.500,
        sigma_0=1.00e+04,
        beta=1.50,
        alpha_res=0.28,
        gamma=1.4,
        delta_H=-945000,
        delta_S=-175,
        n_electrons=4,
        r_eff=2.00e-05,  # Calibrated at T=923K, E=100V/cm
    ),  # DOI: 10.1111/jace.14990
    "Li7La3Zr2O12": MaterialParameters(
        name="LLZO Garnet (Ta-doped)",
        family=MaterialFamily.GARNET,
        Ea=0.304,
        sigma_0=7.65e-01,
        beta=1.50,
        alpha_res=0.20,
        gamma=1.6,
        delta_H=-1750000,
        delta_S=-200,
        n_electrons=4,
        r_eff=8.92e-05,  # Calibrated at T=973K, E=90V/cm
    ),  # DOI: 10.1007/s11581-024-06046-7
    "Li7La3Zr2O12": MaterialParameters(
        name="LLZO Garnet (RFS)",
        family=MaterialFamily.GARNET,
        Ea=0.946,
        sigma_0=1.54e+08,
        beta=1.50,
        alpha_res=0.20,
        gamma=1.6,
        delta_H=-1750000,
        delta_S=-200,
        n_electrons=4,
        r_eff=2.00e-05,  # Calibrated at T=955K, E=40V/cm
    ),  # DOI: 10.1111/jace.16625
    "Ce0.8Sm0.2O1.9": MaterialParameters(
        name="SDC (Samarium-doped Ceria)",
        family=MaterialFamily.FLUORITE,
        Ea=0.570,
        sigma_0=2.12e+02,
        beta=1.69,
        alpha_res=0.15,
        gamma=2.0,
        delta_H=-1024000,
        delta_S=-195,
        n_electrons=4,
        r_eff=7.28e-05,  # Calibrated at T=851K, E=30V/cm
    ),  # DOI: 10.1016/j.scriptamat.2018.09.014
}
