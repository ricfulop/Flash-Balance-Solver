#!/usr/bin/env python3
"""
PDF Data Extraction Pipeline for Flash Balance Solver

This script extracts DOIs, onset temperatures, electric fields, and other
experimental data from PDF papers in the Expanded Solver Papers folder.

Incorporates:
- Paper screening criteria (FULL_DATA, ONSET_ONLY, etc.)
- Confidence levels (H, M, L, D, C)
- Completeness scoring
- Graph digitization flagging

Requirements:
    pip install pdfplumber

Usage:
    python scripts/extract_from_pdfs.py

Output:
    data/paper_screening.csv - Updated with extracted data and categories
    data/pdf_extraction_results.csv - Full extraction details
    data/extracted_dois.csv - DOI list only
    data/papers_needing_digitization.csv - Papers that need graph digitization
"""

import os
import re
import csv
import pdfplumber
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, asdict, field


# =============================================================================
# PAPER CATEGORIES (from EXTRACTION_GUIDE.md)
# =============================================================================

PAPER_CATEGORIES = {
    "FULL_DATA": "Has T_onset, E_field, AND conductivity (Arrhenius) data",
    "ONSET_ONLY": "Has T_onset and E_field, but no Arrhenius plot",
    "CONDUCTIVITY_ONLY": "Has Arrhenius/conductivity data, no flash onset",
    "REVIEW": "Survey/review paper, no original experimental data",
    "THEORY": "Modeling/simulation, no experiments",
    "DIFFERENT_TECHNIQUE": "Uses SPS, microwave, or other non-DC-field method",
    "INSUFFICIENT": "Missing critical data",
}

# Keywords that indicate different paper types
REVIEW_KEYWORDS = ["review", "survey", "progress", "advances in", "state of the art", "overview"]
THEORY_KEYWORDS = ["simulation", "modeling", "theoretical", "first-principles", "DFT", "molecular dynamics"]
SPS_KEYWORDS = ["spark plasma", "SPS", "microwave", "FAST sintering", "current-assisted"]
ARRHENIUS_KEYWORDS = ["arrhenius", "activation energy", "conductivity", "ln(σ)", "log σ", "1000/T", "impedance"]


@dataclass
class PaperData:
    """Data extracted from a single paper."""
    filename: str
    doi: str = ""
    year: str = ""
    first_author: str = ""
    title: str = ""
    
    # Paper category
    paper_category: str = ""
    technique_used: str = "DC_FLASH"
    
    # Data availability flags
    has_onset_data: bool = False
    has_conductivity: bool = False
    has_thermodynamics: bool = False
    has_arrhenius_plot: bool = False
    
    # Extracted values
    T_onset_C: str = ""
    T_onset_K: str = ""
    E_field_Vcm: str = ""
    current_density: str = ""
    heating_rate: str = ""
    
    # Confidence
    T_confidence: str = ""  # H, M, L
    E_confidence: str = ""  # H, M, L
    
    # Completeness
    completeness_score: int = 0
    is_applicable: bool = False
    needs_digitization: bool = False
    
    # Materials
    material: str = ""
    materials_list: str = ""
    num_conditions: int = 0
    
    # Notes
    screening_notes: str = ""
    extraction_notes: str = ""


# =============================================================================
# REGEX PATTERNS FOR DATA EXTRACTION
# =============================================================================

DOI_PATTERN = r'10\.\d{4,9}/[^\s,\)\]\"\'<>]+'

# Temperature patterns
TEMP_PATTERNS = [
    (r'(?:flash|onset)[^\d]{0,30}(\d{3,4})\s*[°º]?\s*C', 'onset_text'),
    (r'T[_\s]?onset\s*[=:≈]\s*(\d{3,4})\s*K', 'T_onset_var'),
    (r'onset\s+temperature[^\d]{0,20}(\d{3,4})\s*[°º]?\s*C', 'onset_temp'),
    (r'flash\s+(?:occurred|initiated|started)[^\d]{0,20}(\d{3,4})\s*[°º]?\s*C', 'flash_occurred'),
    (r'(\d{3,4})\s*[°º]?\s*C[^\d]{0,20}(?:onset|flash)', 'temp_before_onset'),
]

# Electric field patterns  
FIELD_PATTERNS = [
    (r'(\d{2,4})\s*V\s*[/·]?\s*cm', 'V_per_cm'),
    (r'E\s*[=:≈]\s*(\d{2,4})\s*V', 'E_equals'),
    (r'(?:electric\s+)?field[^\d]{0,15}(\d{2,4})\s*V', 'field_of'),
    (r'(\d{2,4})\s*V\s*/\s*cm', 'V_slash_cm'),
]

# Arrhenius/conductivity indicators
ARRHENIUS_PATTERNS = [
    r'arrhenius',
    r'ln\s*\(?σ',
    r'log\s*\(?σ',
    r'activation\s+energy',
    r'1000\s*/\s*T',
    r'impedance\s+spectroscop',
    r'conductivity.*temperature',
    r'E[_a]?\s*[=:]\s*[\d.]+\s*eV',
]


def extract_doi(text: str) -> str:
    """Extract DOI from text."""
    matches = re.findall(DOI_PATTERN, text)
    if matches:
        cleaned = matches[0].rstrip('.,;:')
        return cleaned
    return ""


def extract_year_from_doi(doi: str) -> str:
    """Try to extract publication year from DOI patterns."""
    # Some DOIs contain year info
    year_match = re.search(r'20[12]\d', doi)
    if year_match:
        return year_match.group()
    return ""


def extract_title(text: str) -> str:
    """Extract title from first page text (usually first major text block)."""
    lines = text.split('\n')
    # Title is usually in the first few lines, before "Abstract"
    title_lines = []
    for line in lines[:15]:
        line = line.strip()
        if len(line) > 20 and not any(kw in line.lower() for kw in ['abstract', 'doi:', 'received', 'accepted', 'http']):
            title_lines.append(line)
            if len(' '.join(title_lines)) > 100:
                break
    return ' '.join(title_lines)[:200] if title_lines else ""


def extract_temperatures(text: str) -> Tuple[List[str], str]:
    """Extract temperature values and determine confidence."""
    temps = []
    confidence = "L"  # Default low
    
    for pattern, pattern_type in TEMP_PATTERNS:
        matches = re.findall(pattern, text, re.IGNORECASE)
        if matches:
            temps.extend(matches)
            if pattern_type in ['onset_text', 'T_onset_var', 'onset_temp']:
                confidence = "H"  # High confidence if explicitly stated
            elif confidence != "H":
                confidence = "M"  # Medium if found but not explicit
    
    return list(set(temps)), confidence


def extract_fields(text: str) -> Tuple[List[str], str]:
    """Extract electric field values and determine confidence."""
    fields = []
    confidence = "L"
    
    for pattern, pattern_type in FIELD_PATTERNS:
        matches = re.findall(pattern, text, re.IGNORECASE)
        if matches:
            fields.extend(matches)
            if pattern_type in ['E_equals', 'field_of']:
                confidence = "H"
            elif confidence != "H":
                confidence = "M"
    
    # Filter out unrealistic values (likely page numbers, etc.)
    valid_fields = [f for f in fields if 10 <= int(f) <= 5000]
    return list(set(valid_fields)), confidence


def has_arrhenius_data(text: str) -> bool:
    """Check if paper contains Arrhenius plot or conductivity data."""
    text_lower = text.lower()
    for pattern in ARRHENIUS_PATTERNS:
        if re.search(pattern, text_lower):
            return True
    return False


def determine_category(text: str, has_T: bool, has_E: bool, has_arrhenius: bool) -> Tuple[str, str]:
    """Determine paper category and technique."""
    text_lower = text.lower()
    
    # Check for review/theory papers
    if any(kw in text_lower for kw in REVIEW_KEYWORDS):
        return "REVIEW", "REVIEW"
    
    if any(kw in text_lower for kw in THEORY_KEYWORDS):
        return "THEORY", "SIMULATION"
    
    # Check for different techniques
    technique = "DC_FLASH"
    if any(kw in text_lower for kw in SPS_KEYWORDS):
        technique = "SPS" if "spark plasma" in text_lower or "SPS" in text else "OTHER"
        if not has_T and not has_E:
            return "DIFFERENT_TECHNIQUE", technique
    
    # Determine category based on data availability
    if has_T and has_E and has_arrhenius:
        return "FULL_DATA", technique
    elif has_T and has_E:
        return "ONSET_ONLY", technique
    elif has_arrhenius:
        return "CONDUCTIVITY_ONLY", technique
    else:
        return "INSUFFICIENT", technique


def calculate_completeness(data: PaperData) -> int:
    """Calculate completeness score (0-100)."""
    score = 0
    
    # Critical parameters (25 points each)
    if data.T_onset_C or data.T_onset_K:
        score += 25
    if data.E_field_Vcm:
        score += 25
    
    # Important parameters (15 points each)
    if data.has_arrhenius_plot:
        score += 15  # Ea available
        score += 15  # σ₀ available
    
    # Supplementary (10 points each)
    if data.has_thermodynamics:
        score += 10  # ΔH
        score += 10  # ΔS
    
    return score


def guess_material_from_filename(filename: str) -> str:
    """Extract material name from filename."""
    name = filename.replace('.pdf', '')
    parts = name.split()
    
    # Known author surnames to skip
    authors = ['Raj', 'Todd', 'Wang', 'Sglavo', 'Grasso', 'Kim', 'Jia', 'Yamamoto', 
               'Janssen', 'Bran', 'Jiang', 'Luo', 'Cheng', 'Perejon', 'Maqueda',
               'Tsur', 'Nishino', 'McClellan', 'Lupascu', 'Garcia', 'Steil',
               'Mccartney', 'An', 'Jo', 'Akdogan', 'Senos', 'Huang', 'Xu',
               'Yoshida', 'Muccillo', 'Mazo', 'Valdez', 'Badimele', 'Jalali']
    
    material_parts = []
    for part in parts:
        if part not in authors and not part.isdigit():
            material_parts.append(part)
        if len(material_parts) >= 2:
            break
    
    return ' '.join(material_parts) if material_parts else name


def extract_text_from_pdf(pdf_path: str, max_pages: int = 10) -> str:
    """Extract text content from PDF."""
    try:
        with pdfplumber.open(pdf_path) as pdf:
            text = ""
            for i, page in enumerate(pdf.pages[:max_pages]):
                page_text = page.extract_text()
                if page_text:
                    text += page_text + "\n"
            return text
    except Exception as e:
        return ""


def extract_paper_data(pdf_path: str) -> PaperData:
    """Extract all data from a single paper."""
    filename = os.path.basename(pdf_path)
    data = PaperData(filename=filename)
    
    # Extract text
    text = extract_text_from_pdf(pdf_path)
    if not text:
        data.screening_notes = "Failed to extract text from PDF"
        data.paper_category = "INSUFFICIENT"
        return data
    
    # Extract DOI
    data.doi = extract_doi(text)
    data.year = extract_year_from_doi(data.doi)
    
    # Extract title
    data.title = extract_title(text)[:150]
    
    # Guess material from filename
    data.material = guess_material_from_filename(filename)
    
    # Extract temperatures
    temps, t_conf = extract_temperatures(text)
    if temps:
        # Convert to sorted list
        temp_values = sorted([int(t) for t in temps if t.isdigit()])
        # Filter reasonable values (200-2000°C for flash sintering)
        valid_temps = [t for t in temp_values if 200 <= t <= 2000]
        if valid_temps:
            data.T_onset_C = ", ".join(str(t) for t in valid_temps[:5])  # Top 5
            data.T_confidence = t_conf
            data.has_onset_data = True
            data.extraction_notes = f"Found temps: {valid_temps}"
    
    # Extract electric fields
    fields, e_conf = extract_fields(text)
    if fields:
        field_values = sorted([int(f) for f in fields], key=int)
        data.E_field_Vcm = ", ".join(str(f) for f in field_values[:5])
        data.E_confidence = e_conf
        data.has_onset_data = data.has_onset_data or bool(fields)
    
    # Check for Arrhenius data
    data.has_arrhenius_plot = has_arrhenius_data(text)
    data.has_conductivity = data.has_arrhenius_plot
    
    # Determine category
    has_T = bool(data.T_onset_C)
    has_E = bool(data.E_field_Vcm)
    data.paper_category, data.technique_used = determine_category(
        text, has_T, has_E, data.has_arrhenius_plot
    )
    
    # Calculate completeness
    data.completeness_score = calculate_completeness(data)
    
    # Determine if applicable for solver
    data.is_applicable = (
        data.paper_category in ["FULL_DATA", "ONSET_ONLY"] and
        data.completeness_score >= 50
    )
    
    # Flag if needs digitization (has graphs but values not extracted)
    data.needs_digitization = (
        data.paper_category in ["FULL_DATA", "ONSET_ONLY", "CONDUCTIVITY_ONLY"] and
        (not has_T or not has_E or data.has_arrhenius_plot)
    )
    
    # Count experimental conditions
    if has_T and has_E:
        data.num_conditions = min(len(data.T_onset_C.split(',')), len(data.E_field_Vcm.split(',')))
    
    return data


def process_all_papers(papers_dir: str) -> List[PaperData]:
    """Process all PDFs in the papers directory."""
    results = []
    
    pdf_files = sorted([f for f in os.listdir(papers_dir) if f.endswith('.pdf')])
    total = len(pdf_files)
    
    print(f"\nProcessing {total} PDF files...")
    print("=" * 80)
    
    for i, filename in enumerate(pdf_files, 1):
        pdf_path = os.path.join(papers_dir, filename)
        print(f"[{i}/{total}] {filename[:55]}...")
        
        data = extract_paper_data(pdf_path)
        results.append(data)
        
        # Show status
        doi_status = "✓DOI" if data.doi else "✗DOI"
        cat_status = data.paper_category[:8]
        score_status = f"Score:{data.completeness_score}"
        print(f"        {doi_status} | {cat_status:<10} | {score_status}")
    
    return results


def save_results(results: List[PaperData], output_dir: str):
    """Save extraction results to CSV files."""
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Save full extraction results
    full_path = os.path.join(output_dir, "pdf_extraction_results.csv")
    with open(full_path, 'w', newline='', encoding='utf-8') as f:
        fieldnames = [
            'filename', 'doi', 'material', 'paper_category', 'completeness_score',
            'T_onset_C', 'T_confidence', 'E_field_Vcm', 'E_confidence',
            'has_arrhenius_plot', 'is_applicable', 'needs_digitization',
            'technique_used', 'extraction_notes'
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for data in results:
            row = {k: getattr(data, k, '') for k in fieldnames}
            writer.writerow(row)
    print(f"\nFull results: {full_path}")
    
    # 2. Update paper_screening.csv with extracted data
    screening_path = os.path.join(output_dir, "paper_screening.csv")
    with open(screening_path, 'w', newline='', encoding='utf-8') as f:
        fieldnames = [
            'filename', 'doi', 'year', 'first_author', 'title', 'paper_category',
            'technique_used', 'has_onset_data', 'has_conductivity', 'has_thermodynamics',
            'materials_list', 'num_conditions', 'is_applicable', 'screening_notes'
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for data in results:
            row = {
                'filename': data.filename,
                'doi': data.doi,
                'year': data.year,
                'first_author': '',  # Would need NLP to extract
                'title': data.title,
                'paper_category': data.paper_category,
                'technique_used': data.technique_used,
                'has_onset_data': 'TRUE' if data.has_onset_data else 'FALSE',
                'has_conductivity': 'TRUE' if data.has_conductivity else 'FALSE',
                'has_thermodynamics': 'FALSE',  # Rarely in papers
                'materials_list': data.material,
                'num_conditions': data.num_conditions,
                'is_applicable': 'TRUE' if data.is_applicable else 'FALSE',
                'screening_notes': data.screening_notes or data.extraction_notes,
            }
            writer.writerow(row)
    print(f"Paper screening: {screening_path}")
    
    # 3. Save DOI list
    doi_path = os.path.join(output_dir, "extracted_dois.csv")
    with open(doi_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['filename', 'doi', 'material', 'paper_category'])
        for data in results:
            writer.writerow([data.filename, data.doi, data.material, data.paper_category])
    print(f"DOI list: {doi_path}")
    
    # 4. Save papers needing digitization
    digitize_path = os.path.join(output_dir, "papers_needing_digitization.csv")
    needs_digitization = [d for d in results if d.needs_digitization]
    with open(digitize_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['filename', 'material', 'doi', 'missing_data', 'has_arrhenius', 'priority'])
        for data in needs_digitization:
            missing = []
            if not data.T_onset_C:
                missing.append('T_onset')
            if not data.E_field_Vcm:
                missing.append('E_field')
            if data.has_arrhenius_plot:
                missing.append('Ea_from_plot')
            
            # Priority: FULL_DATA > ONSET_ONLY > CONDUCTIVITY_ONLY
            priority_map = {'FULL_DATA': 1, 'ONSET_ONLY': 2, 'CONDUCTIVITY_ONLY': 3}
            priority = priority_map.get(data.paper_category, 4)
            
            writer.writerow([
                data.filename, data.material, data.doi,
                ', '.join(missing), 'YES' if data.has_arrhenius_plot else 'NO',
                priority
            ])
    print(f"Digitization queue: {digitize_path}")
    
    # Print summary
    print_summary(results)


def print_summary(results: List[PaperData]):
    """Print extraction summary."""
    print("\n" + "=" * 80)
    print("EXTRACTION SUMMARY")
    print("=" * 80)
    
    total = len(results)
    
    # DOI stats
    with_doi = sum(1 for r in results if r.doi)
    print(f"\nDOI Extraction:")
    print(f"  With DOI:     {with_doi}/{total} ({100*with_doi/total:.0f}%)")
    
    # Category breakdown
    print(f"\nPaper Categories:")
    categories = {}
    for r in results:
        cat = r.paper_category
        categories[cat] = categories.get(cat, 0) + 1
    
    for cat in ['FULL_DATA', 'ONSET_ONLY', 'CONDUCTIVITY_ONLY', 'REVIEW', 'THEORY', 'DIFFERENT_TECHNIQUE', 'INSUFFICIENT']:
        count = categories.get(cat, 0)
        if count > 0:
            print(f"  {cat:<20} {count:>3} ({100*count/total:>4.0f}%)")
    
    # Applicable papers
    applicable = sum(1 for r in results if r.is_applicable)
    print(f"\nApplicable for Solver:")
    print(f"  Ready to use:      {applicable}/{total} ({100*applicable/total:.0f}%)")
    
    # Digitization needed
    needs_dig = sum(1 for r in results if r.needs_digitization)
    print(f"  Needs digitization: {needs_dig}/{total} ({100*needs_dig/total:.0f}%)")
    
    # Completeness score distribution
    print(f"\nCompleteness Scores:")
    score_buckets = {'0-25': 0, '26-50': 0, '51-75': 0, '76-100': 0}
    for r in results:
        s = r.completeness_score
        if s <= 25:
            score_buckets['0-25'] += 1
        elif s <= 50:
            score_buckets['26-50'] += 1
        elif s <= 75:
            score_buckets['51-75'] += 1
        else:
            score_buckets['76-100'] += 1
    
    for bucket, count in score_buckets.items():
        bar = '█' * (count // 2)
        print(f"  {bucket:>6}: {bar} {count}")
    
    # Top materials
    print(f"\nMost Common Materials:")
    materials = {}
    for r in results:
        mat = r.material.split()[0] if r.material else 'Unknown'
        materials[mat] = materials.get(mat, 0) + 1
    
    for mat, count in sorted(materials.items(), key=lambda x: -x[1])[:10]:
        print(f"  {mat:<15} {count}")


def main():
    """Main entry point."""
    script_dir = Path(__file__).parent
    project_dir = script_dir.parent
    papers_dir = project_dir / "Expanded Solver Papers"
    output_dir = project_dir / "data"
    
    print("=" * 80)
    print("FLASH BALANCE PDF EXTRACTION PIPELINE v2.0")
    print("=" * 80)
    print(f"Papers directory: {papers_dir}")
    print(f"Output directory: {output_dir}")
    print("\nCategories: FULL_DATA | ONSET_ONLY | CONDUCTIVITY_ONLY | REVIEW | THEORY | INSUFFICIENT")
    print("Confidence: H (High) | M (Medium) | L (Low)")
    
    if not papers_dir.exists():
        print(f"\nERROR: Papers directory not found: {papers_dir}")
        return 1
    
    # Process all papers
    results = process_all_papers(str(papers_dir))
    
    # Save results
    save_results(results, str(output_dir))
    
    print("\n" + "=" * 80)
    print("Extraction complete!")
    print("=" * 80)
    print("\nNext steps:")
    print("1. Review papers_needing_digitization.csv")
    print("2. Use WebPlotDigitizer (https://automeris.io/WebPlotDigitizer/)")
    print("3. Save digitized data to data/digitized/[Material]_[Author]/")
    print("4. Run scripts/analyze_arrhenius.py for Ea extraction")
    print("5. Run scripts/calibrate_materials.py with complete data")
    
    return 0


if __name__ == "__main__":
    exit(main())
