"""
XRD data analyzer: parse, background-subtract, peak-find, phase-match.
Rigaku RAS_XY format supported. Extensible for other formats.
"""
import re
import numpy as np
from numpy.polynomial import polynomial as P
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class XRDData:
    """Parsed XRD measurement."""
    filepath: str
    two_theta: np.ndarray
    intensity: np.ndarray
    metadata: dict = field(default_factory=dict)
    format_name: str = "unknown"

    @property
    def n_points(self) -> int:
        return len(self.two_theta)

    @property
    def range_2theta(self) -> tuple:
        return (float(self.two_theta[0]), float(self.two_theta[-1]))


@dataclass
class XRDPeak:
    """A detected diffraction peak."""
    position: float         # 2θ position
    height: float           # background-subtracted height
    raw_intensity: float    # raw intensity at peak
    fwhm: float             # full width at half maximum
    snr: float              # signal-to-noise ratio
    d_spacing: float = 0.0  # d-spacing in Å (for Cu Kα λ=1.54059)

    def __post_init__(self):
        if self.d_spacing == 0.0:
            self.d_spacing = 1.54059 / (2 * np.sin(np.radians(self.position / 2)))


@dataclass
class PhaseMatch:
    """Result of matching peaks to a reference phase."""
    phase_name: str
    pdf_number: str
    peaks_matched: int
    peaks_total: int
    match_details: list  # list of (ref_pos, obs_pos, delta, snr)


def detect_format(filepath: str) -> str:
    """Auto-detect XRD file format from header."""
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        head = ''.join(f.readline() for _ in range(30))

    if '*FILE_TYPE' in head and 'RAS_XY' in head:
        return 'rigaku_ras_xy'
    if head.strip().startswith('<?xml'):
        return 'bruker_xml'
    if '##$' in head or head.strip().startswith(';'):
        return 'jcamp_dx'
    # Fallback: assume simple 2-column
    return 'two_column'


def parse_rigaku_ras_xy(filepath: str) -> XRDData:
    """Parse Rigaku RAS_XY format."""
    metadata = {}
    data_start = False
    xy_data = []

    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if line.startswith('*'):
                # Metadata line
                parts = line[1:].split(' ', 1)
                if len(parts) == 2:
                    key, val = parts[0], parts[1].strip().strip('"')
                    metadata[key] = val
            elif line == '#Intensity_unit=cps':
                data_start = True
                continue
            elif data_start and line:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        xy_data.append([float(parts[0]), float(parts[1])])
                    except ValueError:
                        continue

    tt = np.array([p[0] for p in xy_data])
    intensity = np.array([p[1] for p in xy_data])
    return XRDData(
        filepath=filepath,
        two_theta=tt,
        intensity=intensity,
        metadata=metadata,
        format_name='rigaku_ras_xy'
    )


def parse_two_column(filepath: str) -> XRDData:
    """Parse simple 2-column (2θ intensity) format."""
    xy_data = []
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('*'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    xy_data.append([float(parts[0]), float(parts[1])])
                except ValueError:
                    continue
    tt = np.array([p[0] for p in xy_data])
    intensity = np.array([p[1] for p in xy_data])
    return XRDData(filepath=filepath, two_theta=tt, intensity=intensity, format_name='two_column')


def parse_xrd(filepath: str) -> XRDData:
    """Parse any supported XRD file format."""
    fmt = detect_format(filepath)
    parsers = {
        'rigaku_ras_xy': parse_rigaku_ras_xy,
        'two_column': parse_two_column,
        'jcamp_dx': parse_two_column,
    }
    parser = parsers.get(fmt, parse_two_column)
    return parser(filepath)


def subtract_background(data: XRDData, poly_deg: int = 6, n_iter: int = 3,
                        sigma_threshold: float = 2.5) -> np.ndarray:
    """
    Iterative polynomial background subtraction.
    Returns: background array, subtracted intensity array.
    """
    tt, intensity = data.two_theta, data.intensity
    mask = np.ones(len(intensity), dtype=bool)

    for _ in range(n_iter):
        try:
            coeffs = P.polyfit(tt[mask], intensity[mask], deg=poly_deg)
            bg = P.polyval(tt, coeffs)
            residuals = intensity - bg
            sigma = np.std(residuals[mask])
            if sigma == 0:
                break
            mask = residuals < sigma_threshold * sigma
        except np.linalg.LinAlgError:
            coeffs = P.polyfit(tt, intensity, deg=min(poly_deg, 3))
            bg = P.polyval(tt, coeffs)
            break

    coeffs = P.polyfit(tt[mask], intensity[mask], deg=poly_deg)
    bg_final = P.polyval(tt, coeffs)
    return bg_final, intensity - bg_final


def find_peaks(data: XRDData, bg_subtracted: np.ndarray,
               sigma_threshold: float = 3.0, min_snr: float = 3.0) -> list[XRDPeak]:
    """
    Find significant peaks in background-subtracted data.
    Returns peaks with SNR > min_snr.
    """
    noise_sigma = np.std(bg_subtracted)
    threshold = sigma_threshold * noise_sigma
    tt = data.two_theta
    raw = data.intensity

    peaks = []
    i = 1
    while i < len(bg_subtracted) - 1:
        if (bg_subtracted[i] > threshold and
                bg_subtracted[i] > bg_subtracted[i - 1] and
                bg_subtracted[i] > bg_subtracted[i + 1]):
            peak_height = bg_subtracted[i]
            snr = peak_height / noise_sigma if noise_sigma > 0 else 0

            if snr < min_snr:
                i += 1
                continue

            # Find FWHM
            half_max = peak_height / 2
            left = i
            while left > 0 and bg_subtracted[left] > half_max:
                left -= 1
            right = i
            while right < len(bg_subtracted) - 1 and bg_subtracted[right] > half_max:
                right += 1
            fwhm = tt[right] - tt[left]

            peak = XRDPeak(
                position=float(tt[i]),
                height=float(peak_height),
                raw_intensity=float(raw[i]),
                fwhm=float(fwhm),
                snr=float(snr)
            )
            peaks.append(peak)
            i = right + 1
        else:
            i += 1

    peaks.sort(key=lambda p: -p.height)
    return peaks


def match_phases(peaks: list[XRDPeak], reference_phases: dict,
                 tolerance: float = 0.3) -> list[PhaseMatch]:
    """
    Match observed peaks against reference phases.
    reference_phases: {phase_name: [(2theta, relative_intensity), ...]}
    """
    results = []
    for phase_name, ref_peaks in reference_phases.items():
        matches = []
        for ref_pos, _ in ref_peaks:
            best_match = None
            best_dist = tolerance + 1
            for p in peaks:
                dist = abs(p.position - ref_pos)
                if dist < tolerance and dist < best_dist:
                    best_dist = dist
                    best_match = (ref_pos, p.position, dist, p.snr)
            if best_match:
                matches.append(best_match)

        results.append(PhaseMatch(
            phase_name=phase_name,
            pdf_number="",
            peaks_matched=len(matches),
            peaks_total=len(ref_peaks),
            match_details=matches
        ))

    # Sort by match ratio descending
    results.sort(key=lambda r: r.peaks_matched / max(r.peaks_total, 1), reverse=True)
    return results


def format_peak_table(peaks: list[XRDPeak]) -> str:
    """Format peaks as a markdown table."""
    lines = ["| # | 2θ (°) | d (Å) | Height | Raw Int. | FWHM (°) | SNR |",
             "|---|--------|-------|--------|----------|----------|-----|"]
    for i, p in enumerate(peaks[:20], 1):
        lines.append(f"| {i} | {p.position:.3f} | {p.d_spacing:.4f} | "
                     f"{p.height:.0f} | {p.raw_intensity:.0f} | "
                     f"{p.fwhm:.3f} | {p.snr:.1f} |")
    return "\n".join(lines)


def format_phase_report(matches: list[PhaseMatch]) -> str:
    """Format phase matching results."""
    lines = ["## Phase Matching Results", ""]
    for m in matches:
        ratio = m.peaks_matched / max(m.peaks_total, 1) * 100
        verdict = "✓ MATCH" if ratio > 50 else ("△ PARTIAL" if ratio > 20 else "✗ NO")
        lines.append(f"### {verdict}: {m.phase_name} ({ratio:.0f}%)")
        lines.append(f"Matched {m.peaks_matched}/{m.peaks_total} peaks")
        for ref_pos, obs_pos, delta, snr in m.match_details:
            lines.append(f"  - {ref_pos:.3f}° → {obs_pos:.3f}° (Δ={delta:.3f}°, SNR={snr:.1f})")
        lines.append("")
    return "\n".join(lines)


def analyze(filepath: str, expected_phases: dict, poly_deg: int = 6,
            sigma_threshold: float = 3.0, min_snr: float = 3.0) -> dict:
    """
    Complete XRD analysis pipeline.

    Args:
        filepath: Path to XRD data file.
        expected_phases: Dict of {phase_name: [(2theta, rel_int), ...]}.
            Must be user-confirmed. Empty dict = skip phase matching.
        poly_deg: Polynomial degree for background fit.
        sigma_threshold: Peak detection threshold (multiples of noise sigma).
        min_snr: Minimum SNR for peak reporting.

    Returns:
        dict with keys: data, background, bg_subtracted, peaks, phase_matches,
                        peak_table, phase_report
    """
    data = parse_xrd(filepath)
    background, bg_subtracted = subtract_background(data, poly_deg=poly_deg)
    peaks = find_peaks(data, bg_subtracted, sigma_threshold=sigma_threshold,
                       min_snr=min_snr)

    result = {
        "data": data,
        "background": background,
        "bg_subtracted": bg_subtracted,
        "peaks": peaks,
        "peak_table": format_peak_table(peaks),
        "metadata_summary": _format_metadata(data),
    }

    if expected_phases:
        phase_matches = match_phases(peaks, expected_phases)
        result["phase_matches"] = phase_matches
        result["phase_report"] = format_phase_report(phase_matches)
    else:
        result["phase_matches"] = []
        result["phase_report"] = ""

    return result


def _format_metadata(data: XRDData) -> str:
    """Format measurement metadata for display."""
    md = data.metadata
    lines = [f"**Format**: {data.format_name}", f"**Data points**: {data.n_points}",
             f"**2θ range**: {data.range_2theta[0]:.1f}° – {data.range_2theta[1]:.1f}°"]
    if 'MEAS_COND_XG_CURRENT' in md and 'MEAS_COND_XG_VOLTAGE' in md:
        lines.append(f"**X-ray**: {md.get('HW_XG_TARGET_NAME','?')} "
                     f"{md['MEAS_COND_XG_VOLTAGE']}kV/{md['MEAS_COND_XG_CURRENT']}mA")
    if 'MEAS_SCAN_START_TIME' in md:
        lines.append(f"**Date**: {md['MEAS_SCAN_START_TIME']}")
    if 'MEAS_SCAN_SPEED' in md:
        lines.append(f"**Speed**: {md['MEAS_SCAN_SPEED']} {md.get('MEAS_SCAN_SPEED_UNIT','')}")
    if 'MEAS_SCAN_STEP' in md:
        lines.append(f"**Step**: {md['MEAS_SCAN_STEP']}°")
    return "\n".join(lines)
