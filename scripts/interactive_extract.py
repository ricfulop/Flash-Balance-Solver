#!/usr/bin/env python3
"""
Interactive Paper Extraction for Flash Balance Solver - Full Datapoint Version

Extracts ALL datapoints from each paper, classifies them by type,
and stores all conditions in the database.

Usage:
    python scripts/interactive_extract.py           # Start from current progress
    python scripts/interactive_extract.py --paper 5 # Start at paper #5
    python scripts/interactive_extract.py --reset   # Reset progress

Datapoint Types:
    PRE_FLASH    - Before ignition (Arrhenius/conductivity data)
    AT_ONSET     - Flash ignition point (primary for solver)
    STEADY_STATE - During flash hold
    POST_FLASH   - After flash completed
    ARRHENIUS    - Conductivity measurement (no field/flash)
"""

import sys
import json
import re
import argparse
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, asdict, field
import warnings

warnings.filterwarnings('ignore')

sys.path.insert(0, str(Path(__file__).parent.parent))

import fitz  # PyMuPDF
import numpy as np

from flash_balance_solver import MaterialParameters, MaterialFamily, FlashBalanceSolver
from scipy.optimize import brentq


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass
class Datapoint:
    """A single experimental datapoint from a paper."""
    condition_id: int
    datapoint_type: str  # PRE_FLASH, AT_ONSET, STEADY_STATE, POST_FLASH, ARRHENIUS
    
    # Core flash data
    T_K: Optional[float] = None
    T_C: Optional[float] = None
    E_field_Vcm: Optional[float] = None
    J_Acm2: Optional[float] = None
    J_mAmm2: Optional[float] = None
    P_Wcm3: Optional[float] = None
    V_applied: Optional[float] = None
    I_mA: Optional[float] = None
    
    # Experimental conditions
    heating_rate_Cmin: Optional[float] = None
    hold_time_s: Optional[float] = None
    atmosphere: Optional[str] = None
    
    # Sample info
    sample_geometry: Optional[str] = None
    gauge_length_mm: Optional[float] = None
    cross_section_mm2: Optional[float] = None
    green_density_pct: Optional[float] = None
    final_density_pct: Optional[float] = None
    grain_size_um: Optional[float] = None
    
    # Conductivity data (for Arrhenius)
    sigma_Sm: Optional[float] = None
    sigma_Scm: Optional[float] = None
    
    # Source tracking
    source: str = ""  # "Table 1", "Fig 3", "text p.4"
    page: Optional[int] = None
    confidence: str = "M"  # H, M, L
    
    # Classification
    use_for_calibration: bool = False
    notes: str = ""


@dataclass 
class PaperExtraction:
    """Complete extraction from a single paper."""
    paper_num: int
    filename: str
    doi: str
    
    # Material info
    material: str
    material_formula: str
    material_family: str
    is_composite: bool = False
    constituents: List[str] = field(default_factory=list)
    
    # Paper classification
    category: str = ""  # FULL_DATA, ONSET_ONLY, REVIEW, etc.
    completeness_score: int = 0
    
    # Arrhenius parameters (if extracted)
    Ea_eV: Optional[float] = None
    Ea_source: str = ""
    Ea_confidence: str = "L"
    sigma_0_Sm: Optional[float] = None
    sigma_0_source: str = ""
    sigma_0_confidence: str = "L"
    
    # Thermodynamic data
    delta_H: Optional[float] = None
    delta_S: Optional[float] = None
    delta_H_source: str = "NIST"
    
    # All datapoints
    datapoints: List[Datapoint] = field(default_factory=list)
    
    # Validation
    r_eff_m: Optional[float] = None
    validation_error_pct: Optional[float] = None
    
    # Metadata
    extraction_date: str = ""
    notes: str = ""


# =============================================================================
# PATHS AND CONSTANTS
# =============================================================================

DATA_DIR = Path(__file__).parent.parent / "data"
PAPERS_DIR = Path(__file__).parent.parent / "Expanded Solver Papers"
APPROVED_DB = DATA_DIR / "approved_materials.json"
PROGRESS_FILE = DATA_DIR / "extraction_progress.json"

FAMILY_DEFAULTS = {
    "fluorite": {"beta": 1.69, "alpha_res": 0.15, "gamma": 2.0, "delta_H": -1085000, "delta_S": -178},
    "perovskite": {"beta": 1.50, "alpha_res": 0.28, "gamma": 1.4, "delta_H": -1660000, "delta_S": -192},
    "spinel": {"beta": 1.20, "alpha_res": 0.20, "gamma": 1.8, "delta_H": -1676000, "delta_S": -210},
    "garnet": {"beta": 1.50, "alpha_res": 0.20, "gamma": 1.6, "delta_H": -1750000, "delta_S": -200},
    "rutile": {"beta": 1.43, "alpha_res": 0.22, "gamma": 1.8, "delta_H": -944000, "delta_S": -186},
    "wurtzite": {"beta": 1.30, "alpha_res": 0.18, "gamma": 1.5, "delta_H": -350000, "delta_S": -44},
    "carbide": {"beta": 1.10, "alpha_res": 0.25, "gamma": 1.2, "delta_H": -40000, "delta_S": -10},
    "nitride": {"beta": 1.15, "alpha_res": 0.22, "gamma": 1.3, "delta_H": -300000, "delta_S": -100},
}


# =============================================================================
# DATABASE FUNCTIONS
# =============================================================================

def load_database() -> Dict:
    """Load the approved materials database."""
    if APPROVED_DB.exists():
        with open(APPROVED_DB) as f:
            return json.load(f)
    return {"metadata": {}, "papers": [], "rejected": []}


def save_database(db: Dict):
    """Save the approved materials database."""
    db["metadata"]["last_updated"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Count total datapoints
    total_dp = sum(len(p.get('datapoints', [])) for p in db.get('papers', []))
    db["metadata"]["total_datapoints"] = total_dp
    with open(APPROVED_DB, 'w') as f:
        json.dump(db, f, indent=2)


def load_progress() -> Dict:
    if PROGRESS_FILE.exists():
        with open(PROGRESS_FILE) as f:
            return json.load(f)
    return {"current_paper": 1, "total_papers": 83}


def save_progress(progress: Dict):
    with open(PROGRESS_FILE, 'w') as f:
        json.dump(progress, f, indent=2)


# =============================================================================
# EXTRACTION FUNCTIONS
# =============================================================================

def get_paper_list() -> List[str]:
    """Get sorted list of PDF files."""
    return sorted([f.name for f in PAPERS_DIR.glob("*.pdf")])


def extract_doi(pdf_path: Path) -> str:
    """Extract DOI from PDF - returns FULL DOI."""
    try:
        with fitz.open(pdf_path) as doc:
            for page_num in range(min(3, len(doc))):
                text = doc[page_num].get_text()
                # Look for DOI pattern - get the FULL match
                doi_match = re.search(r'(10\.\d{4,}/[^\s\]>\)\"\'\,]+)', text)
                if doi_match:
                    doi = doi_match.group(1)
                    # Clean trailing punctuation
                    doi = doi.rstrip('.,;:')
                    return doi
    except:
        pass
    return ""


def extract_material_info(filename: str) -> Tuple[str, str, str]:
    """
    Extract material name, formula, and family from filename.
    Returns (name, formula, family).
    """
    name = filename.replace('.pdf', '')
    
    patterns = {
        r'(\d*Y?SZ|YSZ|ZrO2)': ('YSZ', 'ZrO2', 'fluorite'),
        r'(LLZO|Li7La3Zr2O12)': ('LLZO', 'Li7La3Zr2O12', 'garnet'),
        r'(GDC|Ce.*O2|SDC|Sm.*Ce)': ('Doped Ceria', 'CeO2', 'fluorite'),
        r'(BaTiO3)': ('BaTiO3', 'BaTiO3', 'perovskite'),
        r'(SrTiO3)': ('SrTiO3', 'SrTiO3', 'perovskite'),
        r'(BiFeO3)': ('BiFeO3', 'BiFeO3', 'perovskite'),
        r'(KNN|K.*Na.*Nb)': ('KNN', 'K0.5Na0.5NbO3', 'perovskite'),
        r'(NaNbO3)': ('NaNbO3', 'NaNbO3', 'perovskite'),
        r'(Al2O3|Alumina)': ('Al2O3', 'Al2O3', 'spinel'),
        r'(MgAl2O4)': ('MgAl2O4', 'MgAl2O4', 'spinel'),
        r'(TiO2|Titania)': ('TiO2', 'TiO2', 'rutile'),
        r'(SnO2)': ('SnO2', 'SnO2', 'rutile'),
        r'(ZnO)': ('ZnO', 'ZnO', 'wurtzite'),
        r'(SiC)': ('SiC', 'SiC', 'carbide'),
        r'(WC)': ('WC', 'WC', 'carbide'),
        r'(B4C)': ('B4C', 'B4C', 'carbide'),
        r'(TiN|ZrN)': ('Nitride', 'TiN', 'nitride'),
    }
    
    for pattern, (mat_name, formula, family) in patterns.items():
        if re.search(pattern, name, re.IGNORECASE):
            # Get more specific name from filename
            specific = name.split()[0] if ' ' in name else name
            return specific, formula, family
    
    parts = name.split()
    return parts[0] if parts else name, "", 'fluorite'


def classify_datapoint(T_K: Optional[float], E_Vcm: Optional[float], 
                       J_Acm2: Optional[float] = None,
                       context: str = "") -> str:
    """
    Classify a datapoint by type based on conditions.
    """
    context_lower = context.lower()
    
    # Check context clues first
    if 'onset' in context_lower or 'flash' in context_lower:
        if 'before' in context_lower or 'pre' in context_lower:
            return "PRE_FLASH"
        elif 'after' in context_lower or 'post' in context_lower:
            return "POST_FLASH"
        else:
            return "AT_ONSET"
    
    if 'arrhenius' in context_lower or 'conductivity' in context_lower:
        return "ARRHENIUS"
    
    if 'steady' in context_lower or 'hold' in context_lower:
        return "STEADY_STATE"
    
    # Use heuristics based on values
    if E_Vcm is None or E_Vcm < 10:
        return "ARRHENIUS"
    
    if T_K:
        # Very high temperatures often post-flash
        if T_K > 1500 and E_Vcm and E_Vcm < 200:
            return "POST_FLASH"
        # Moderate temps with moderate field = likely onset
        if 700 < T_K < 1400 and E_Vcm and 30 < E_Vcm < 500:
            return "AT_ONSET"
    
    if J_Acm2:
        # Very high current density = steady state or post
        if J_Acm2 > 50:
            return "STEADY_STATE"
    
    # Default
    return "AT_ONSET"


def extract_all_datapoints(pdf_path: Path) -> Tuple[List[Datapoint], Dict]:
    """
    Extract ALL datapoints from a PDF.
    Returns (list of datapoints, metadata dict).
    """
    datapoints = []
    metadata = {
        'T_candidates': [],
        'E_candidates': [],
        'J_candidates': [],
        'has_arrhenius': False,
        'has_tables': False,
        'text_snippets': [],
        'heating_rates': [],
        'atmospheres': [],
        'densities': [],
        'grain_sizes': [],
    }
    
    try:
        with fitz.open(pdf_path) as doc:
            full_text = ""
            for page_num, page in enumerate(doc):
                page_text = page.get_text()
                full_text += page_text + "\n"
                
                # Extract temperatures with context - multiple patterns
                # Handle various degree symbols: ° º ◦ o ˚
                temp_c_patterns = [
                    r'(\d{3,4})\s*[°º◦˚o]?\s*C(?:elsius)?(?![a-zA-Z])',  # 900°C, 900 C
                    r'(\d{3,4})\s*[°º◦˚]\s*C',  # 900 ° C with various degree symbols
                    r'[Tt](?:emperature)?\s*[=:]\s*(\d{3,4})',  # T = 900, Temperature: 900
                    r'at\s+(\d{3,4})\s*[°º◦˚]?\s*C',  # at 900°C
                    r'(\d{3,4})\s*◦\s*C',  # Specifically for ◦ character
                ]
                
                for pattern in temp_c_patterns:
                    temp_c_matches = list(re.finditer(pattern, page_text, re.IGNORECASE))
                    for m in temp_c_matches:
                        t_val = int(m.group(1))
                        if 200 <= t_val <= 1800:
                            # Get context (50 chars before and after)
                            start = max(0, m.start() - 50)
                            end = min(len(page_text), m.end() + 50)
                            context = page_text[start:end]
                            
                            metadata['T_candidates'].append({
                                'value': t_val,
                                'unit': 'C',
                                'K': t_val + 273,
                                'page': page_num + 1,
                                'context': context.replace('\n', ' ')[:100]
                            })
                
                # Extract electric fields with context
                field_matches = list(re.finditer(
                    r'(\d{1,4})\s*V\s*[/·]?\s*cm',
                    page_text, re.IGNORECASE
                ))
                for m in field_matches:
                    e_val = int(m.group(1))
                    if 5 <= e_val <= 5000:
                        start = max(0, m.start() - 50)
                        end = min(len(page_text), m.end() + 50)
                        context = page_text[start:end]
                        
                        metadata['E_candidates'].append({
                            'value': e_val,
                            'page': page_num + 1,
                            'context': context.replace('\n', ' ')[:100]
                        })
                
                # Extract current density
                j_matches = list(re.finditer(
                    r'(\d+\.?\d*)\s*(?:mA\s*/?\s*mm|A\s*/?\s*cm)',
                    page_text, re.IGNORECASE
                ))
                for m in j_matches:
                    j_val = float(m.group(1))
                    start = max(0, m.start() - 50)
                    end = min(len(page_text), m.end() + 50)
                    context = page_text[start:end]
                    
                    metadata['J_candidates'].append({
                        'value': j_val,
                        'unit': 'mA/mm2' if 'mA' in m.group() else 'A/cm2',
                        'page': page_num + 1,
                        'context': context.replace('\n', ' ')[:100]
                    })
                
                # Extract heating rates
                hr_matches = re.findall(r'(\d+)\s*[°º]?\s*C\s*/\s*min', page_text, re.IGNORECASE)
                for hr in hr_matches:
                    hr_val = int(hr)
                    if 1 <= hr_val <= 100:
                        metadata['heating_rates'].append(hr_val)
                
                # Extract densities
                dens_matches = re.findall(r'(\d{2,3})\.?\d*\s*%\s*(?:TD|theoretical|relative|dense)', 
                                         page_text, re.IGNORECASE)
                for d in dens_matches:
                    d_val = int(d)
                    if 40 <= d_val <= 100:
                        metadata['densities'].append(d_val)
                
                # Extract grain sizes
                gs_matches = re.findall(r'(\d+\.?\d*)\s*[μµ]m', page_text)
                for gs in gs_matches:
                    gs_val = float(gs)
                    if 0.01 <= gs_val <= 100:
                        metadata['grain_sizes'].append(gs_val)
            
            # Check for Arrhenius
            arrhenius_kw = ['arrhenius', 'activation energy', 'ln(σ', 'log(σ', 'e_a', 'ea =']
            metadata['has_arrhenius'] = any(kw in full_text.lower() for kw in arrhenius_kw)
            
            # Check for tables
            metadata['has_tables'] = 'table' in full_text.lower()
            
            # Extract Ea values directly
            ea_matches = re.findall(r'[Ee][_a]?\s*[=:≈]\s*([\d.]+)\s*eV', full_text)
            metadata['Ea_values'] = [float(e) for e in ea_matches if 0.1 < float(e) < 3.0]
            
            # Extract atmosphere
            if 'air' in full_text.lower():
                metadata['atmospheres'].append('air')
            if 'argon' in full_text.lower() or 'ar' in full_text.lower():
                metadata['atmospheres'].append('Ar')
            if 'nitrogen' in full_text.lower() or 'n2' in full_text.lower():
                metadata['atmospheres'].append('N2')
            if 'vacuum' in full_text.lower():
                metadata['atmospheres'].append('vacuum')
    
    except Exception as e:
        metadata['error'] = str(e)
    
    # Remove duplicates and sort
    seen_T = set()
    unique_T = []
    for t in metadata['T_candidates']:
        if t['K'] not in seen_T:
            seen_T.add(t['K'])
            unique_T.append(t)
    metadata['T_candidates'] = sorted(unique_T, key=lambda x: x['K'])
    
    seen_E = set()
    unique_E = []
    for e in metadata['E_candidates']:
        if e['value'] not in seen_E:
            seen_E.add(e['value'])
            unique_E.append(e)
    metadata['E_candidates'] = sorted(unique_E, key=lambda x: x['value'])
    
    # Create datapoints from combinations
    condition_id = 1
    
    # Create datapoints for each T-E combination that appears together
    for t_data in metadata['T_candidates'][:20]:  # Limit to avoid explosion
        for e_data in metadata['E_candidates'][:10]:
            # Check if they appear on same page or in similar context
            same_page = t_data['page'] == e_data['page']
            
            if same_page:
                dp_type = classify_datapoint(
                    t_data['K'], 
                    e_data['value'],
                    context=t_data['context'] + ' ' + e_data['context']
                )
                
                dp = Datapoint(
                    condition_id=condition_id,
                    datapoint_type=dp_type,
                    T_K=t_data['K'],
                    T_C=t_data['value'],
                    E_field_Vcm=e_data['value'],
                    source=f"Page {t_data['page']}",
                    page=t_data['page'],
                    confidence="M" if same_page else "L",
                    notes=f"T context: {t_data['context'][:50]}..."
                )
                
                # Add current density if found on same page
                for j in metadata['J_candidates']:
                    if j['page'] == t_data['page']:
                        if j['unit'] == 'mA/mm2':
                            dp.J_mAmm2 = j['value']
                            dp.J_Acm2 = j['value'] / 100  # Convert
                        else:
                            dp.J_Acm2 = j['value']
                        break
                
                # Add experimental conditions
                if metadata['heating_rates']:
                    dp.heating_rate_Cmin = metadata['heating_rates'][0]
                if metadata['atmospheres']:
                    dp.atmosphere = metadata['atmospheres'][0]
                if metadata['densities']:
                    dp.green_density_pct = min(metadata['densities'])
                    dp.final_density_pct = max(metadata['densities'])
                if metadata['grain_sizes']:
                    dp.grain_size_um = min(metadata['grain_sizes'])
                
                datapoints.append(dp)
                condition_id += 1
    
    return datapoints, metadata


def categorize_paper(metadata: Dict, datapoints: List[Datapoint]) -> Tuple[str, int]:
    """Categorize paper and calculate completeness score."""
    score = 0
    
    has_T = len(metadata.get('T_candidates', [])) > 0
    has_E = len(metadata.get('E_candidates', [])) > 0
    has_J = len(metadata.get('J_candidates', [])) > 0
    has_arrhenius = metadata.get('has_arrhenius', False)
    has_onset_points = any(dp.datapoint_type == 'AT_ONSET' for dp in datapoints)
    
    if has_T: score += 25
    if has_E: score += 25
    if has_J: score += 10
    if has_arrhenius: score += 20
    if has_onset_points: score += 20
    
    if has_T and has_E and has_arrhenius and has_onset_points:
        category = "FULL_DATA"
    elif has_T and has_E and has_onset_points:
        category = "ONSET_ONLY"
    elif has_arrhenius and not has_onset_points:
        category = "CONDUCTIVITY_ONLY"
    elif score < 30:
        category = "INSUFFICIENT"
    else:
        category = "REVIEW"
    
    return category, score


# =============================================================================
# DISPLAY FUNCTIONS
# =============================================================================

def print_header(paper_num: int, total: int, filename: str):
    print("\n" + "=" * 90)
    print(f"PAPER #{paper_num}/{total}: {filename}")
    print("=" * 90)


def print_datapoints_table(datapoints: List[Datapoint], extraction: 'PaperExtraction'):
    """Print all datapoints in comprehensive validation table format."""
    print("\n" + "-" * 180)
    print("EXTRACTED DATAPOINTS - FULL VALIDATION FORMAT")
    print("-" * 180)
    
    # Header
    print(f"{'#':<2} | {'Family':<10} | {'Category':<12} | {'Material':<15} | {'Ea':<5} | {'β':<5} | {'α_res':<5} | "
          f"{'r_eff':<7} | {'k_soft':<6} | {'n':<2} | {'T(K)':<6} | {'E':<6} | {'T_pred':<6} | {'Err%':<6} | DOI")
    print("-" * 180)
    
    # Get family defaults for parameters
    family = extraction.material_family.lower()
    defaults = FAMILY_DEFAULTS.get(family, FAMILY_DEFAULTS['fluorite'])
    
    Ea = extraction.Ea_eV or 0.9
    beta = defaults['beta']
    alpha_res = defaults['alpha_res']
    n = 4
    
    # Calculate k_soft
    q_ratio = 0.73
    k_soft = max(0.01, 1 - beta * (q_ratio ** 2))
    
    # Short DOI (but longer than before)
    doi_short = extraction.doi[:35] + "..." if len(extraction.doi) > 35 else extraction.doi
    
    for dp in datapoints:
        t_k = f"{dp.T_K:.0f}" if dp.T_K else "-"
        e = f"{dp.E_field_Vcm:.0f}" if dp.E_field_Vcm else "-"
        
        # Calculate T_pred and error if we have data
        t_pred = "-"
        err = "-"
        r_eff_str = "-"
        
        if dp.T_K and dp.E_field_Vcm and Ea:
            # Try to calculate prediction
            try:
                family_map = {
                    'fluorite': MaterialFamily.FLUORITE,
                    'perovskite': MaterialFamily.PEROVSKITE,
                    'spinel': MaterialFamily.SPINEL,
                    'garnet': MaterialFamily.FLUORITE,
                    'rutile': MaterialFamily.RUTILE,
                    'wurtzite': MaterialFamily.WURTZITE,
                    'carbide': MaterialFamily.CARBIDE,
                    'nitride': MaterialFamily.NITRIDE,
                }
                
                sigma_0 = extraction.sigma_0_Sm or 1e4
                
                # Calibrate r_eff for this point
                params = MaterialParameters(
                    name=extraction.material,
                    family=family_map.get(family, MaterialFamily.FLUORITE),
                    Ea=Ea,
                    sigma_0=sigma_0,
                    beta=beta,
                    alpha_res=alpha_res,
                    gamma=defaults['gamma'],
                    delta_H=extraction.delta_H or defaults['delta_H'],
                    delta_S=extraction.delta_S or defaults['delta_S'],
                    n_electrons=n,
                    r_eff=20e-6,
                )
                
                # Calibrate
                E_field = dp.E_field_Vcm * 100
                
                def objective(r_eff):
                    test_p = MaterialParameters(
                        name=params.name, family=params.family, Ea=params.Ea,
                        sigma_0=params.sigma_0, beta=params.beta, alpha_res=params.alpha_res,
                        gamma=params.gamma, delta_H=params.delta_H, delta_S=params.delta_S,
                        n_electrons=params.n_electrons, r_eff=r_eff,
                    )
                    solver = FlashBalanceSolver(test_p)
                    T_p = solver.solve_onset_temperature(E_field)
                    return (T_p - dp.T_K) if T_p else 1e6
                
                try:
                    r_eff_cal = brentq(objective, 1e-9, 1e-2, rtol=1e-4)
                    r_eff_str = f"{r_eff_cal*1e6:.1f}"
                    
                    # Validate
                    params_cal = MaterialParameters(
                        name=params.name, family=params.family, Ea=params.Ea,
                        sigma_0=params.sigma_0, beta=params.beta, alpha_res=params.alpha_res,
                        gamma=params.gamma, delta_H=params.delta_H, delta_S=params.delta_S,
                        n_electrons=params.n_electrons, r_eff=r_eff_cal,
                    )
                    solver = FlashBalanceSolver(params_cal)
                    T_predicted = solver.solve_onset_temperature(E_field)
                    
                    if T_predicted:
                        t_pred = f"{T_predicted:.0f}"
                        error_pct = 100 * (T_predicted - dp.T_K) / dp.T_K
                        err = f"{error_pct:+.1f}%"
                        dp.use_for_calibration = abs(error_pct) < 15
                except:
                    r_eff_str = "20.0*"  # Failed calibration
            except Exception as e:
                pass
        
        # Print row
        print(f"{dp.condition_id:<2} | {extraction.material_family:<10} | {extraction.category:<12} | "
              f"{extraction.material_formula:<15} | {Ea:<5.2f} | {beta:<5.2f} | {alpha_res:<5.2f} | "
              f"{r_eff_str:<7} | {k_soft:<6.3f} | {n:<2} | {t_k:<6} | {e:<6} | {t_pred:<6} | {err:<6} | {doi_short}")
    
    print("-" * 180)
    print(f"Total: {len(datapoints)} datapoints")


def print_extraction_summary(extraction: PaperExtraction):
    """Print full extraction summary."""
    print("\n" + "=" * 90)
    print("EXTRACTION SUMMARY")
    print("=" * 90)
    print(f"  Paper #:        {extraction.paper_num}")
    print(f"  Filename:       {extraction.filename}")
    print(f"  DOI:            {extraction.doi}")
    print(f"  Material:       {extraction.material}")
    print(f"  Formula:        {extraction.material_formula}")
    print(f"  Family:         {extraction.material_family}")
    print(f"  Category:       {extraction.category}")
    print(f"  Score:          {extraction.completeness_score}/100")
    print("-" * 90)
    
    # Count datapoints by type
    type_counts = {}
    for dp in extraction.datapoints:
        type_counts[dp.datapoint_type] = type_counts.get(dp.datapoint_type, 0) + 1
    
    print("  Datapoints by type:")
    for dtype, count in sorted(type_counts.items()):
        print(f"    {dtype}: {count}")
    
    print("-" * 90)
    if extraction.Ea_eV:
        print(f"  Ea:             {extraction.Ea_eV:.3f} eV ({extraction.Ea_confidence})")
    if extraction.sigma_0_Sm:
        print(f"  σ₀:             {extraction.sigma_0_Sm:.2e} S/m")
    print("=" * 90)


# =============================================================================
# MAIN WORKFLOW
# =============================================================================

def process_single_paper(paper_num: int, filename: str, total: int) -> PaperExtraction:
    """Process a single paper - extract all datapoints."""
    pdf_path = PAPERS_DIR / filename
    
    print_header(paper_num, total, filename)
    
    # Initialize extraction
    extraction = PaperExtraction(
        paper_num=paper_num,
        filename=filename,
        doi="",
        material="",
        material_formula="",
        material_family="",
        extraction_date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    )
    
    # Step 1: Extract DOI
    print("\n--- STEP 1: DOI ---")
    doi = extract_doi(pdf_path)
    extraction.doi = doi
    print(f"  DOI: {doi if doi else 'Not found'}")
    
    # Step 2: Material info
    print("\n--- STEP 2: MATERIAL ---")
    name, formula, family = extract_material_info(filename)
    extraction.material = name
    extraction.material_formula = formula
    extraction.material_family = family
    print(f"  Name: {name}")
    print(f"  Formula: {formula}")
    print(f"  Family: {family}")
    
    # Step 3: Extract all datapoints
    print("\n--- STEP 3: EXTRACT ALL DATAPOINTS ---")
    datapoints, metadata = extract_all_datapoints(pdf_path)
    extraction.datapoints = datapoints
    
    print(f"  Temperatures found: {len(metadata['T_candidates'])}")
    print(f"  E-fields found: {len(metadata['E_candidates'])}")
    print(f"  Current densities found: {len(metadata['J_candidates'])}")
    print(f"  Has Arrhenius data: {metadata['has_arrhenius']}")
    
    # Print sample candidates
    if metadata['T_candidates']:
        print("\n  Temperature candidates (first 5):")
        for t in metadata['T_candidates'][:5]:
            print(f"    {t['value']}°C = {t['K']}K (p.{t['page']})")
    
    if metadata['E_candidates']:
        print("\n  E-field candidates (first 5):")
        for e in metadata['E_candidates'][:5]:
            print(f"    {e['value']} V/cm (p.{e['page']})")
    
    # Step 4: Categorize
    print("\n--- STEP 4: CATEGORIZE ---")
    category, score = categorize_paper(metadata, datapoints)
    extraction.category = category
    extraction.completeness_score = score
    print(f"  Category: {category}")
    print(f"  Completeness: {score}/100")
    
    # Step 5: Arrhenius
    print("\n--- STEP 5: ARRHENIUS ---")
    if metadata.get('Ea_values'):
        extraction.Ea_eV = metadata['Ea_values'][0]
        extraction.Ea_source = "text"
        extraction.Ea_confidence = "M"
        print(f"  Ea found in text: {extraction.Ea_eV} eV")
    elif metadata['has_arrhenius']:
        # Use family default
        defaults = FAMILY_DEFAULTS.get(family, FAMILY_DEFAULTS['fluorite'])
        extraction.Ea_eV = 0.9 if family == 'fluorite' else 0.5
        extraction.Ea_source = "family_default"
        extraction.Ea_confidence = "L"
        print(f"  Ea (default): {extraction.Ea_eV} eV")
    
    # Set thermodynamic defaults
    defaults = FAMILY_DEFAULTS.get(family, FAMILY_DEFAULTS['fluorite'])
    extraction.delta_H = defaults['delta_H']
    extraction.delta_S = defaults['delta_S']
    
    # Print datapoints table
    if datapoints:
        print_datapoints_table(datapoints, extraction)
    else:
        print("\n  No datapoints could be paired (T + E on same page)")
    
    # Print summary
    print_extraction_summary(extraction)
    
    return extraction


def approve_extraction(extraction: PaperExtraction, db: Dict) -> Dict:
    """Add approved extraction to database."""
    # Convert to dict for JSON storage
    entry = {
        "paper_num": extraction.paper_num,
        "filename": extraction.filename,
        "doi": extraction.doi,
        "material": extraction.material,
        "material_formula": extraction.material_formula,
        "material_family": extraction.material_family,
        "category": extraction.category,
        "completeness_score": extraction.completeness_score,
        "Ea_eV": extraction.Ea_eV,
        "Ea_source": extraction.Ea_source,
        "sigma_0_Sm": extraction.sigma_0_Sm,
        "delta_H": extraction.delta_H,
        "delta_S": extraction.delta_S,
        "datapoints": [asdict(dp) for dp in extraction.datapoints],
        "extraction_date": extraction.extraction_date,
        "approved": True,
    }
    
    db['papers'].append(entry)
    db['metadata']['papers_approved'] = len(db['papers'])
    db['metadata']['papers_processed'] = len(db['papers']) + len(db.get('rejected', []))
    
    save_database(db)
    print(f"\n✓ Paper #{extraction.paper_num} approved and saved to database")
    print(f"  Total datapoints: {len(extraction.datapoints)}")
    return db


def reject_extraction(extraction: PaperExtraction, reason: str, db: Dict) -> Dict:
    """Add rejected paper to database."""
    entry = {
        "paper_num": extraction.paper_num,
        "filename": extraction.filename,
        "doi": extraction.doi,
        "material": extraction.material,
        "category": extraction.category,
        "reason": reason,
        "rejection_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    }
    
    db['rejected'].append(entry)
    db['metadata']['papers_rejected'] = len(db['rejected'])
    db['metadata']['papers_processed'] = len(db['papers']) + len(db['rejected'])
    
    save_database(db)
    print(f"\n✗ Paper #{extraction.paper_num} rejected: {reason}")
    return db


def run_interactive(start_paper: int = 1):
    """Run the interactive extraction workflow."""
    papers = get_paper_list()
    total = len(papers)
    db = load_database()
    
    print("\n" + "=" * 90)
    print("FLASH BALANCE SOLVER - INTERACTIVE FULL EXTRACTION")
    print("=" * 90)
    print(f"Total papers: {total}")
    print(f"Starting at: #{start_paper}")
    print(f"Already approved: {len(db.get('papers', []))}")
    print(f"Already rejected: {len(db.get('rejected', []))}")
    print("=" * 90)
    
    current = start_paper
    filename = papers[current - 1]
    
    # Check if already done
    already_done = any(
        p.get('paper_num') == current or p.get('filename') == filename
        for p in db.get('papers', []) + db.get('rejected', [])
    )
    
    if already_done:
        print(f"\n⏭ Paper #{current} already processed")
        return None, current, papers, db
    
    extraction = process_single_paper(current, filename, total)
    
    # Save progress
    progress = {
        "current_paper": current,
        "total_papers": total,
        "last_processed": filename,
        "session_time": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    }
    save_progress(progress)
    
    return extraction, current, papers, db


def main():
    parser = argparse.ArgumentParser(description='Interactive full extraction')
    parser.add_argument('--paper', '-p', type=int, default=None)
    parser.add_argument('--reset', action='store_true')
    args = parser.parse_args()
    
    if args.reset:
        save_progress({"current_paper": 1, "total_papers": 83})
        print("Progress reset")
        return
    
    progress = load_progress()
    start = args.paper if args.paper else progress.get('current_paper', 1)
    
    extraction, current, papers, db = run_interactive(start)
    
    if extraction:
        print("\n" + "-" * 90)
        print("APPROVAL REQUIRED")
        print("-" * 90)
        print("  [A] Approve all datapoints")
        print("  [R] Reject paper (with reason)")
        print("  [S] Skip for now")
        print("-" * 90)


if __name__ == "__main__":
    main()
