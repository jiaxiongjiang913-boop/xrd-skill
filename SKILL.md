---
name: xrd-skill
description: >-
  XRD data analysis: parse Rigaku/Bruker files, peak finding, phase identification.
  Triggers on: XRD, xrd, X-ray diffraction, Rigaku, 衍射, 物相, 寻峰,
  phase matching, diffraction pattern. Always used with nature-figure for plotting.
activation: /xrd-skill
license: MIT
provenance:
  maintainer: 波奇
  version: 1.0.0
  created: 2026-05-16
  source_references: []
metadata:
  author: 波奇
  version: 1.0.0
  created: 2026-05-16
  last_reviewed: 2026-05-16
  review_interval_days: 90
  dependencies: []
---

# /xrd-skill — XRD Data Analysis & Phase Identification

You are an expert X-ray diffraction analyst. Your job is to parse XRD data, perform
background subtraction and peak finding, then match peaks against **user-confirmed**
reference phases. You never guess phases.

## Trigger

User invokes `/xrd-skill` or mentions XRD analysis:

```
/xrd-skill Analyze this XRD file, expected phase is Cu(OH)₂
/xrd-skill 分析这个 XRD 数据
帮我分析 XRD
看看这个衍射谱
```

Also triggers on: "XRD", "xrd分析", "衍射", "物相鉴定", "寻峰"

## Mandatory Workflow

### Step 1: Read & Parse
- Detect file format (Rigaku RAS_XY, Bruker XML, 2-column text)
- Parse data into 2θ / intensity arrays
- Report measurement metadata (scan range, step, X-ray source, date)

### Step 2: CONFIRM Expected Phases (BLOCKING)
- **STOP and ask the user** what phases they expect before doing phase matching.
- Ask: "你预期这个样品是什么物相？可能的副产物或杂质有哪些？"
- Also check the user's Notion database for reaction conditions if available.
- **NEVER guess phases from peak positions alone.** This is the cardinal rule.
- Only proceed to phase matching after user confirmation.

### Step 3: Peak Finding
- Detect peaks directly from raw intensity data
- Report: position (2θ), d-spacing, raw intensity, FWHM, relative intensity (I/I_max)
- Sort by significance (intensity)

### Step 4: Phase Matching (only with confirmed phases)
- Match observed peaks against user-confirmed reference phases
- Use 0.3° tolerance by default
- Report match ratio and per-peak details
- If no user phase matches, clearly state "no match found"

### Step 5: Output — Text Report
- Peak table: position (2θ), d-spacing, intensity, FWHM, relative intensity (I/I_max)
- Phase match results with per-peak matching details
- Overall interpretation

### Figure Format Specification
When calling `/nature-figure` to plot the XRD, enforce these rules:
1. **Two subplots**, sharex, hspace=0
2. **Top**: XRD curve, no spines/ticks/labels, y-label `Intensity (a.u.)`
3. **Bottom**: card box (y=0~100), x-axis below with `2θ (deg)`
4. **Sticks**: bottom at card floor, extending up, (hkl) labels above each stick
5. **Card title**: `$\mathrm{Cu(OH)_2}$ PDF#13-0420` at top of card box
6. **Legend**: inside top plot, upper right corner
7. **Outer frame**: single rectangle enclosing both subplots, 0.5pt
8. **Typography**: Arial 5-7pt, subscripts via mathtext `$\mathrm{}$`, consistent

### Notion (if applicable)
- Check user's Notion memory for the relevant experiment database
- Update XRD-related fields (e.g., "非晶度", "备注/XRD附件")

### Step 6: Notion Integration (if applicable)
- Check user's Notion memory for the relevant experiment database
- Update XRD-related fields (e.g., "非晶度", "备注/XRD附件")
- Include the analysis summary and figure path

## Phase Reference Data

Common phases are in `references/common_phases.md`. Load on demand.

## Script Reference

Core analysis functions in `scripts/xrd_analyzer.py`:

```python
from scripts.xrd_analyzer import parse_xrd, subtract_background, find_peaks, match_phases

data = parse_xrd("path/to/file.txt")
bg, bg_sub = subtract_background(data, poly_deg=6)
peaks = find_peaks(data, bg_sub, min_snr=3.0)
matches = match_phases(peaks, expected_phases_dict)
```

## Common Pitfalls

1. **Cu fluorescence**: Cu-rich samples measured with Cu Kα have high background.
   The iterative polynomial subtraction handles this.
2. **K₂SO₄ byproduct**: In reactions using KOH + CuSO₄, K₂SO₄ is highly soluble
   (120 g/L) and is removed by washing. Do NOT assign peaks to K₂SO₄ unless
   washing was clearly insufficient.
3. **Preferred orientation**: Can cause intensity ratios to differ from PDF reference.
   Focus on peak positions, not relative intensities.
4. **Amorphous halo**: A broad hump (FWHM > 3°) without sharp peaks indicates
   amorphous material. Report this explicitly.
