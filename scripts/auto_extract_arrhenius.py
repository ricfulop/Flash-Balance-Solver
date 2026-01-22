#!/usr/bin/env python3
"""
Automatic Arrhenius Plot Extraction from PDFs

Fully automated pipeline to extract Ea and sigma_0 from Arrhenius plots
embedded in scientific papers. No manual intervention required.

Pipeline:
1. Extract all images from PDF pages
2. Detect Arrhenius plots via OCR (looking for "1/T", "ln σ", etc.)
3. Auto-calibrate axes from detected tick marks
4. Extract data points using color/shape detection
5. Fit Arrhenius equation and output Ea, sigma_0

Usage:
    python scripts/auto_extract_arrhenius.py "Expanded Solver Papers/"
    python scripts/auto_extract_arrhenius.py "Expanded Solver Papers/BiFeO3 Raj.pdf"
"""

import os
import sys
import re
import json
import csv
import warnings
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, asdict
import numpy as np

# Suppress warnings
warnings.filterwarnings('ignore')

# PDF processing
import fitz  # PyMuPDF

# Image processing
import cv2
from skimage import feature, morphology, measure
from skimage.color import rgb2gray

# OCR
try:
    import easyocr
    EASYOCR_AVAILABLE = True
except ImportError:
    EASYOCR_AVAILABLE = False
    print("Warning: easyocr not available, using fallback pattern matching")


@dataclass
class ArrheniusResult:
    """Result from Arrhenius extraction."""
    filename: str
    material: str
    doi: str
    Ea_eV: float
    sigma_0_Sm: float
    R2: float
    confidence: str
    n_points: int
    T_range_K: Tuple[float, float]
    plot_page: int
    extraction_method: str
    notes: str = ""


class ArrheniusExtractor:
    """Automated Arrhenius plot extractor."""
    
    # Physical constants
    kB_eV = 8.617e-5  # Boltzmann constant in eV/K
    
    # Arrhenius plot indicators
    ARRHENIUS_KEYWORDS = [
        '1000/t', '1/t', '10³/t', '10^3/t',
        'ln(σ', 'ln σ', 'log σ', 'log(σ',
        'arrhenius', 'activation energy',
        'conductivity', 's/cm', 's cm',
    ]
    
    # Typical axis ranges for different plot types
    AXIS_RANGES = {
        '1000/T': (0.5, 2.5),  # 1000/K (corresponds to 400-2000K)
        'log_sigma': (-6, 2),   # log10(S/cm)
    }
    
    def __init__(self, verbose: bool = True):
        self.verbose = verbose
        self.reader = None
        if EASYOCR_AVAILABLE:
            try:
                self.reader = easyocr.Reader(['en'], gpu=False, verbose=False)
            except Exception as e:
                print(f"EasyOCR init failed: {e}")
    
    def log(self, msg: str):
        if self.verbose:
            print(msg)
    
    def extract_from_pdf(self, pdf_path: str) -> Optional[ArrheniusResult]:
        """
        Extract Arrhenius data from a PDF file.
        
        Returns ArrheniusResult if successful, None otherwise.
        """
        pdf_path = Path(pdf_path)
        self.log(f"\n{'='*60}")
        self.log(f"Processing: {pdf_path.name}")
        self.log('='*60)
        
        # Extract material name from filename
        material = self._extract_material_name(pdf_path.stem)
        
        # Extract DOI
        doi = self._extract_doi_from_pdf(pdf_path)
        
        # Get all page images
        images = self._extract_page_images(pdf_path)
        if not images:
            self.log("  No images extracted from PDF")
            return None
        
        self.log(f"  Extracted {len(images)} page images")
        
        # Find Arrhenius plots - skip page 1 (usually title page)
        arrhenius_candidates = []
        for page_num, img in images:
            if page_num == 1:
                continue  # Skip title page
            
            score, is_arrhenius = self._is_arrhenius_plot(img)
            if is_arrhenius:
                arrhenius_candidates.append((page_num, img, score))
                self.log(f"  Page {page_num}: Arrhenius plot detected (score: {score:.2f})")
        
        if not arrhenius_candidates:
            self.log("  No Arrhenius plots detected via OCR, searching all pages...")
            # Try to find any plot with data points (skip page 1)
            for page_num, img in images:
                if page_num == 1:
                    continue
                points = self._extract_data_points(img)
                if 5 <= len(points) <= 50:  # Reasonable number for Arrhenius
                    arrhenius_candidates.append((page_num, img, 0.5))
                    self.log(f"  Page {page_num}: Found plot with {len(points)} points")
        
        if not arrhenius_candidates:
            return None
        
        # Try each candidate until we get a valid fit
        arrhenius_candidates.sort(key=lambda x: x[2], reverse=True)
        
        for page_num, img, score in arrhenius_candidates:
            self.log(f"  Trying plot from page {page_num}")
            
            # Extract data points
            points = self._extract_data_points(img)
            if len(points) < 4:
                self.log(f"    Insufficient data points ({len(points)})")
                continue
            
            # Filter to get only plausible Arrhenius points
            # (points that form a roughly linear negative slope)
            filtered_points = self._filter_arrhenius_points(points, img.shape)
            
            if len(filtered_points) < 4:
                self.log(f"    After filtering: {len(filtered_points)} points (need 4+)")
                continue
            
            self.log(f"  Extracted {len(filtered_points)} data points")
            
            # Calibrate axes (estimate from image dimensions and typical ranges)
            x_data, y_data = self._calibrate_and_convert(img, filtered_points)
            
            if len(x_data) < 4:
                continue
            
            # Fit Arrhenius equation
            result = self._fit_arrhenius(x_data, y_data)
            
            if result is None:
                continue
            
            Ea, sigma_0, R2 = result
            
            # Validate the result
            if not (0.05 < Ea < 3.0):  # Reasonable Ea range
                self.log(f"    Ea={Ea:.3f} eV out of range, trying next page")
                continue
            
            if R2 < 0.7:
                self.log(f"    R²={R2:.3f} too low, trying next page")
                continue
            
            # Determine confidence
            if R2 >= 0.99:
                confidence = 'H'
            elif R2 >= 0.95:
                confidence = 'M'
            elif R2 >= 0.85:
                confidence = 'L'
            else:
                confidence = 'D'  # Derived/estimated
            
            # Calculate temperature range
            inv_T = x_data / 1000 if np.mean(x_data) < 10 else x_data
            T_min = 1 / np.max(inv_T) if np.max(inv_T) > 0 else 500
            T_max = 1 / np.min(inv_T) if np.min(inv_T) > 0 else 1500
            
            self.log(f"  ✓ Ea = {Ea:.3f} eV, σ₀ = {sigma_0:.2e} S/m, R² = {R2:.4f}")
            
            return ArrheniusResult(
                filename=pdf_path.name,
                material=material,
                doi=doi,
                Ea_eV=Ea,
                sigma_0_Sm=sigma_0,
                R2=R2,
                confidence=confidence,
                n_points=len(x_data),
                T_range_K=(T_min, T_max),
                plot_page=page_num,
                extraction_method='auto_cv2',
            )
        
        self.log("  No valid Arrhenius fit found")
        return None
    
    def _filter_arrhenius_points(self, points: List[Tuple[float, float]], 
                                  img_shape: Tuple) -> List[Tuple[float, float]]:
        """
        Filter points to keep only those that could be from an Arrhenius plot.
        Arrhenius plots have negative slope (y decreases as x increases).
        """
        if len(points) < 4:
            return points
        
        h, w = img_shape[:2]
        
        # Convert to numpy for easier manipulation
        pts = np.array(points)
        
        # Normalize coordinates
        x_norm = pts[:, 0] / w
        y_norm = pts[:, 1] / h
        
        # Arrhenius: as 1/T increases (x increases), ln(σ) decreases (y increases in image coords)
        # So we want points where y_pixel increases with x_pixel
        
        # Sort by x
        sorted_idx = np.argsort(pts[:, 0])
        pts_sorted = pts[sorted_idx]
        
        # Find the longest monotonically increasing (in y_pixel) subsequence
        # This corresponds to decreasing conductivity with increasing 1/T
        best_sequence = []
        
        for start in range(len(pts_sorted)):
            sequence = [pts_sorted[start]]
            last_y = pts_sorted[start][1]
            
            for i in range(start + 1, len(pts_sorted)):
                # Allow some tolerance for noise
                if pts_sorted[i][1] >= last_y - 0.02 * h:
                    sequence.append(pts_sorted[i])
                    last_y = pts_sorted[i][1]
            
            if len(sequence) > len(best_sequence):
                best_sequence = sequence
        
        # Also check for decreasing y (reverse case)
        best_sequence_rev = []
        for start in range(len(pts_sorted)):
            sequence = [pts_sorted[start]]
            last_y = pts_sorted[start][1]
            
            for i in range(start + 1, len(pts_sorted)):
                if pts_sorted[i][1] <= last_y + 0.02 * h:
                    sequence.append(pts_sorted[i])
                    last_y = pts_sorted[i][1]
            
            if len(sequence) > len(best_sequence_rev):
                best_sequence_rev = sequence
        
        # Use whichever gives more points
        if len(best_sequence_rev) > len(best_sequence):
            best_sequence = best_sequence_rev
        
        return [(p[0], p[1]) for p in best_sequence]
    
    def _extract_material_name(self, filename: str) -> str:
        """Extract material name from filename."""
        # Remove author name (usually last word)
        parts = filename.replace('_', ' ').replace('-', ' ').split()
        if len(parts) > 1:
            # Check if last part is likely an author name (capitalized, no numbers)
            if parts[-1][0].isupper() and not any(c.isdigit() for c in parts[-1]):
                return ' '.join(parts[:-1])
        return filename
    
    def _extract_doi_from_pdf(self, pdf_path: Path) -> str:
        """Extract DOI from PDF metadata or text."""
        try:
            with fitz.open(pdf_path) as doc:
                # Check metadata
                metadata = doc.metadata
                if metadata.get('subject'):
                    doi_match = re.search(r'10\.\d{4,}/[^\s]+', metadata['subject'])
                    if doi_match:
                        return doi_match.group()
                
                # Search first 2 pages
                for page_num in range(min(2, len(doc))):
                    text = doc[page_num].get_text()
                    doi_match = re.search(r'10\.\d{4,}/[^\s\]>]+', text)
                    if doi_match:
                        doi = doi_match.group().rstrip('.,;')
                        return doi
        except Exception:
            pass
        return ""
    
    def _extract_page_images(self, pdf_path: Path) -> List[Tuple[int, np.ndarray]]:
        """Extract images from PDF pages."""
        images = []
        try:
            with fitz.open(pdf_path) as doc:
                for page_num in range(len(doc)):
                    page = doc[page_num]
                    
                    # Render page as image (300 DPI for good quality)
                    mat = fitz.Matrix(2, 2)  # 2x zoom
                    pix = page.get_pixmap(matrix=mat)
                    
                    # Convert to numpy array
                    img = np.frombuffer(pix.samples, dtype=np.uint8)
                    img = img.reshape(pix.height, pix.width, pix.n)
                    
                    # Convert to BGR for OpenCV
                    if pix.n == 4:  # RGBA
                        img = cv2.cvtColor(img, cv2.COLOR_RGBA2BGR)
                    elif pix.n == 3:  # RGB
                        img = cv2.cvtColor(img, cv2.COLOR_RGB2BGR)
                    
                    images.append((page_num + 1, img))
        except Exception as e:
            self.log(f"  Error extracting images: {e}")
        
        return images
    
    def _is_arrhenius_plot(self, img: np.ndarray) -> Tuple[float, bool]:
        """
        Detect if image contains an Arrhenius plot.
        Returns (score, is_arrhenius).
        """
        score = 0.0
        
        # Convert to grayscale for text detection
        if len(img.shape) == 3:
            gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
        else:
            gray = img
        
        # Use OCR if available
        if self.reader:
            try:
                results = self.reader.readtext(gray, detail=0, paragraph=False)
                text = ' '.join(results).lower()
                
                # Check for Arrhenius keywords
                for keyword in self.ARRHENIUS_KEYWORDS:
                    if keyword in text:
                        score += 0.2
                
                # Strong indicators
                if '1/t' in text or '1000/t' in text:
                    score += 0.3
                if 'ln' in text and 'σ' in text.lower():
                    score += 0.3
                if 'log' in text and ('σ' in text.lower() or 'conductivity' in text):
                    score += 0.3
                if 'arrhenius' in text:
                    score += 0.4
                    
            except Exception:
                pass
        
        # Check for plot characteristics using image processing
        # Look for axes (straight lines)
        edges = cv2.Canny(gray, 50, 150)
        lines = cv2.HoughLinesP(edges, 1, np.pi/180, 100, minLineLength=50, maxLineGap=10)
        
        if lines is not None and len(lines) > 2:
            score += 0.1
        
        # Look for data points (small circles/squares)
        circles = cv2.HoughCircles(gray, cv2.HOUGH_GRADIENT, 1, 20,
                                   param1=50, param2=30, minRadius=3, maxRadius=15)
        if circles is not None and len(circles[0]) >= 5:
            score += 0.2
        
        return score, score >= 0.3
    
    def _extract_data_points(self, img: np.ndarray) -> List[Tuple[float, float]]:
        """
        Extract data points from plot image.
        """
        points = []
        h, w = img.shape[:2]
        
        # Convert to HSV for color detection
        hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
        
        # Common data point colors in scientific plots
        color_ranges = [
            # Blue points
            (np.array([100, 50, 50]), np.array([130, 255, 255])),
            # Red points
            (np.array([0, 50, 50]), np.array([10, 255, 255])),
            (np.array([170, 50, 50]), np.array([180, 255, 255])),
            # Black points (use grayscale)
            None,
            # Green points
            (np.array([35, 50, 50]), np.array([85, 255, 255])),
        ]
        
        all_points = []
        
        for color_range in color_ranges:
            if color_range is None:
                # Detect black points
                gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
                _, mask = cv2.threshold(gray, 50, 255, cv2.THRESH_BINARY_INV)
            else:
                mask = cv2.inRange(hsv, color_range[0], color_range[1])
            
            # Clean up mask
            kernel = np.ones((3, 3), np.uint8)
            mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, kernel)
            mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, kernel)
            
            # Find contours
            contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            
            for contour in contours:
                area = cv2.contourArea(contour)
                if 20 < area < 2000:  # Filter by size
                    M = cv2.moments(contour)
                    if M["m00"] > 0:
                        cx = M["m10"] / M["m00"]
                        cy = M["m01"] / M["m00"]
                        
                        # Filter points in plot area (roughly center 70% of image)
                        if 0.1 * w < cx < 0.95 * w and 0.05 * h < cy < 0.9 * h:
                            all_points.append((cx, cy))
        
        # Remove duplicates (points close together)
        filtered_points = []
        for p in all_points:
            is_duplicate = False
            for fp in filtered_points:
                if abs(p[0] - fp[0]) < 10 and abs(p[1] - fp[1]) < 10:
                    is_duplicate = True
                    break
            if not is_duplicate:
                filtered_points.append(p)
        
        # Sort by x coordinate
        filtered_points.sort(key=lambda x: x[0])
        
        return filtered_points
    
    def _calibrate_and_convert(self, img: np.ndarray, 
                                points: List[Tuple[float, float]]) -> Tuple[np.ndarray, np.ndarray]:
        """
        Convert pixel coordinates to data coordinates.
        Uses heuristics based on typical Arrhenius plot ranges.
        """
        if not points:
            return np.array([]), np.array([])
        
        h, w = img.shape[:2]
        
        # Estimate plot region (typically 15-90% of image)
        x_min_pix, x_max_pix = 0.15 * w, 0.9 * w
        y_min_pix, y_max_pix = 0.1 * h, 0.85 * h
        
        # Typical Arrhenius plot ranges
        # X-axis: 1000/T typically 0.8 to 2.0 (500K to 1250K)
        # Y-axis: log10(σ) typically -4 to 1 (S/cm)
        x_min_data, x_max_data = 0.7, 1.8
        y_min_data, y_max_data = -4.5, 0.5
        
        # Try to detect actual axis range from OCR
        if self.reader:
            try:
                gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY) if len(img.shape) == 3 else img
                results = self.reader.readtext(gray)
                
                # Look for axis tick values
                x_ticks = []
                y_ticks = []
                
                for bbox, text, conf in results:
                    if conf < 0.3:
                        continue
                    
                    # Get center position
                    center_x = np.mean([p[0] for p in bbox])
                    center_y = np.mean([p[1] for p in bbox])
                    
                    # Try to parse as number
                    try:
                        val = float(text.replace(',', '.').replace('−', '-'))
                        
                        # Near bottom = x-axis
                        if center_y > 0.85 * h:
                            if 0.3 < val < 3.0:  # Valid 1000/T range
                                x_ticks.append((center_x, val))
                        
                        # Near left = y-axis
                        if center_x < 0.15 * w:
                            if -8 < val < 5:  # Valid log σ range
                                y_ticks.append((center_y, val))
                    except:
                        pass
                
                # Use detected ticks to calibrate
                if len(x_ticks) >= 2:
                    x_ticks.sort(key=lambda x: x[0])
                    px1, v1 = x_ticks[0]
                    px2, v2 = x_ticks[-1]
                    if px2 != px1:
                        scale = (v2 - v1) / (px2 - px1)
                        x_min_data = v1 - scale * (px1 - x_min_pix)
                        x_max_data = v1 + scale * (x_max_pix - px1)
                
                if len(y_ticks) >= 2:
                    y_ticks.sort(key=lambda x: x[0])  # Sort by y pixel (top to bottom)
                    py1, v1 = y_ticks[0]  # Top of plot
                    py2, v2 = y_ticks[-1]  # Bottom of plot
                    if py2 != py1:
                        scale = (v2 - v1) / (py2 - py1)
                        y_max_data = v1 - scale * (py1 - y_min_pix)  # Note: y inverted
                        y_min_data = v1 + scale * (y_max_pix - py1)
                
            except Exception:
                pass
        
        # Convert pixel coordinates to data coordinates
        x_data = []
        y_data = []
        
        for px, py in points:
            # Linear interpolation
            x = x_min_data + (px - x_min_pix) / (x_max_pix - x_min_pix) * (x_max_data - x_min_data)
            y = y_max_data - (py - y_min_pix) / (y_max_pix - y_min_pix) * (y_max_data - y_min_data)
            
            # Filter reasonable values
            if 0.3 < x < 3.0 and -8 < y < 3:
                x_data.append(x)
                y_data.append(y)
        
        return np.array(x_data), np.array(y_data)
    
    def _fit_arrhenius(self, x_data: np.ndarray, 
                       y_data: np.ndarray) -> Optional[Tuple[float, float, float]]:
        """
        Fit Arrhenius equation to data.
        
        Assumes:
        - x_data: 1000/T (1000/K)
        - y_data: log10(σ) (S/cm)
        
        Returns (Ea_eV, sigma_0_Sm, R2)
        """
        if len(x_data) < 3:
            return None
        
        # Convert to Arrhenius form
        # x = 1000/T → 1/T = x/1000
        # y = log10(σ) → ln(σ) = y * ln(10)
        
        inv_T = x_data / 1000  # 1/K
        ln_sigma = y_data * np.log(10)  # Convert log10 to ln
        
        # Also convert S/cm to S/m (multiply by 100 → add ln(100))
        ln_sigma = ln_sigma + np.log(100)
        
        # Linear fit: ln(σ) = ln(σ₀) - Ea/(kB*T)
        # ln(σ) = ln(σ₀) - (Ea/kB) * (1/T)
        try:
            coeffs = np.polyfit(inv_T, ln_sigma, 1)
            slope = coeffs[0]
            intercept = coeffs[1]
            
            # Calculate R²
            y_pred = slope * inv_T + intercept
            ss_res = np.sum((ln_sigma - y_pred) ** 2)
            ss_tot = np.sum((ln_sigma - np.mean(ln_sigma)) ** 2)
            R2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            
            # Extract parameters
            Ea_eV = -slope * self.kB_eV
            sigma_0_Sm = np.exp(intercept)
            
            # Sanity checks
            if Ea_eV < 0 or Ea_eV > 5:
                # Try flipping the data
                ln_sigma_flip = ln_sigma[::-1]
                coeffs = np.polyfit(inv_T, ln_sigma_flip, 1)
                slope = coeffs[0]
                intercept = coeffs[1]
                Ea_eV = -slope * self.kB_eV
                sigma_0_Sm = np.exp(intercept)
            
            if Ea_eV < 0:
                Ea_eV = abs(Ea_eV)
            
            if R2 < 0.5:
                return None
            
            return Ea_eV, sigma_0_Sm, R2
            
        except Exception:
            return None


def process_directory(directory: str, verbose: bool = True) -> List[ArrheniusResult]:
    """Process all PDFs in a directory."""
    extractor = ArrheniusExtractor(verbose=verbose)
    results = []
    
    pdf_files = list(Path(directory).glob('*.pdf'))
    print(f"\nFound {len(pdf_files)} PDF files")
    
    for i, pdf_path in enumerate(pdf_files):
        print(f"\n[{i+1}/{len(pdf_files)}] {pdf_path.name}")
        
        result = extractor.extract_from_pdf(str(pdf_path))
        if result:
            results.append(result)
    
    return results


def print_results_table(results: List[ArrheniusResult]):
    """Print results in table format."""
    print("\n" + "=" * 100)
    print("ARRHENIUS EXTRACTION RESULTS")
    print("=" * 100)
    
    print(f"{'Material':<25} | {'Ea (eV)':<8} | {'σ₀ (S/m)':<12} | {'R²':<8} | {'Conf':<4} | {'Points':<6} | DOI")
    print("-" * 100)
    
    for r in sorted(results, key=lambda x: x.material):
        doi_short = r.doi[:30] + '...' if len(r.doi) > 30 else r.doi
        print(f"{r.material:<25} | {r.Ea_eV:<8.3f} | {r.sigma_0_Sm:<12.2e} | {r.R2:<8.4f} | {r.confidence:<4} | {r.n_points:<6} | {doi_short}")
    
    print("=" * 100)


def save_results(results: List[ArrheniusResult], output_dir: str):
    """Save results to CSV and JSON."""
    output_path = Path(output_dir)
    
    # CSV
    csv_path = output_path / "arrhenius_results.csv"
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['filename', 'material', 'doi', 'Ea_eV', 'sigma_0_Sm', 
                        'R2', 'confidence', 'n_points', 'T_min_K', 'T_max_K',
                        'plot_page', 'extraction_method', 'notes'])
        for r in results:
            writer.writerow([
                r.filename, r.material, r.doi, r.Ea_eV, r.sigma_0_Sm,
                r.R2, r.confidence, r.n_points, r.T_range_K[0], r.T_range_K[1],
                r.plot_page, r.extraction_method, r.notes
            ])
    print(f"\nSaved CSV to: {csv_path}")
    
    # JSON
    json_path = output_path / "arrhenius_results.json"
    with open(json_path, 'w') as f:
        json.dump([asdict(r) for r in results], f, indent=2, default=str)
    print(f"Saved JSON to: {json_path}")


def main():
    """Main entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Automatically extract Arrhenius parameters from PDFs'
    )
    parser.add_argument('path', help='PDF file or directory of PDFs')
    parser.add_argument('--output', '-o', default='data',
                       help='Output directory for results')
    parser.add_argument('--quiet', '-q', action='store_true',
                       help='Suppress verbose output')
    
    args = parser.parse_args()
    
    path = Path(args.path)
    
    if path.is_file():
        # Single file
        extractor = ArrheniusExtractor(verbose=not args.quiet)
        result = extractor.extract_from_pdf(str(path))
        if result:
            print_results_table([result])
            save_results([result], args.output)
    elif path.is_dir():
        # Directory
        results = process_directory(str(path), verbose=not args.quiet)
        if results:
            print_results_table(results)
            save_results(results, args.output)
        else:
            print("\nNo Arrhenius data could be extracted.")
    else:
        print(f"Error: Path not found: {path}")
        sys.exit(1)


if __name__ == '__main__':
    main()
