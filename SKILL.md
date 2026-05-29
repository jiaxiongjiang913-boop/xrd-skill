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

---

## IRON RULE — Figure Format Specification (READ BEFORE WRITING ANY PLOTTING CODE)

**Every XRD figure you produce MUST follow this exact spec. No deviation. No "improvisation."
No "Nature-style" from a different skill. This is the single source of truth.**

### Layout
1. **Two subplots**, sharex, hspace=0
2. Height ratio: top ≈0.62, bottom ≈0.38 (card box is compact)
3. **Axis box frames** on both subplots — all 4 spines visible per panel, 0.4pt `#333333`
   - The axes themselves form the "方框". Do NOT add a separate outer rectangle

### Top panel — XRD Curve
4. No y-axis ticks or tick labels (`labelleft=False, left=False`)
5. All 4 spines visible (box frame), 0.4pt `#333333`
6. y-label: `Intensity (a.u.)`, fontsize=7
7. Sample label: text at upper-left of plot (not legend box), e.g.
   `r'$\mathregular{Sn_{2}S_{3}}$ bulk piece'`

### Bottom panel — Standard Card Box
8. All 4 spines visible (box frame), 0.4pt `#333333` — same as top panel
9. Card bottom edge = x-axis spine (ylim bottom = 0), no gap
10. No y-axis ticks/labels on card panel
11. **Card box height** = 54 data units (compact); adjust `figsize` height accordingly

### Stick patterns
12. **All stick heights / 2** (relative intensity values halved for compact display)
13. Sticks rooted at card floor or divider line, no (hkl) labels
14. **Multi-phase**: divide card vertically with a solid gray line (`#bbbbbb`, 0.5pt)
15. Phase labels at top of each section with PDF#, e.g.:
    - `r"$\mathregular{Sn_{2}S_{3}}$ PDF#27-0899"`

### Typography
16. Arial / DejaVu Sans, 5-7pt throughout
17. Chemical formulas via mathtext: `$\mathregular{Sn_{2}S_{3}}$` (handles Unicode subscripts)

### rcParams (MUST set before any plotting)
```python
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'svg.fonttype': 'none',
    'font.size': 7,
    'axes.linewidth': 0.4,
    'xtick.major.width': 0.4,
    'ytick.major.width': 0.4,
})
```

### Save
```python
fig.savefig("name.svg", bbox_inches='tight', dpi=300, facecolor='white')
fig.savefig("name.png", bbox_inches='tight', dpi=300, facecolor='white')
plt.close(fig)
```

---

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

### Step 6: Plotting — BLOCKING re-read + compliance gate

**Before writing ANY plotting code, you MUST:**

1. **Re-read the IRON RULE** section at the top of this file (lines 31–84). Do this every time,
   even if you think you remember it. Read it fresh.

2. **Write the plotting script**, matching every rule exactly — spine visibility, frame color,
   font size, label positions, stick heights, card box height, everything.

3. **Self-audit against the 17 rules** before running the script. For each rule, confirm:
   - Rule 1 (two subplots sharex hspace=0)? Yes/No
   - Rule 2 (height ratio 0.62/0.38)? Yes/No
   - Rule 3 (axis box frames, NO outer rectangle)? Yes/No
   - Rules 4-7 (top panel)? Yes/No
   - Rules 8-11 (bottom panel)? Yes/No
   - Rules 12-15 (stick patterns)? Yes/No
   - Rules 16-17 (typography)? Yes/No

4. **Only after all 17 rules pass**, execute the script and save SVG + PNG.

5. If any rule fails, fix the code and re-audit. Do not proceed with a failing audit.

### Step 7: Notion Integration (if applicable)
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
