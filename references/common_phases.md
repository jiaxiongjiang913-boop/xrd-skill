# Common XRD Reference Phases

Load this file when doing phase matching. Each phase lists 2θ positions
for Cu Kα radiation (λ = 1.54059 Å) and relative intensities (100 = strongest).

## Cu Compounds

### Cu(OH)₂ — Orthorhombic (PDF #13-0420)
```
2θ (°)     d (Å)     hkl    I/I₀
16.72      5.298     020     40
23.84      3.729     021     60
33.97      2.637     130     70
35.90      2.499     111     50
39.80      2.264     131     30
53.29      1.718     150     40
55.15      1.664     220     15
```
peaks: [(16.72, 40), (23.84, 60), (33.97, 70), (35.90, 50), (39.80, 30), (53.29, 40)]

### CuO — Tenorite, Monoclinic (PDF #48-1548)
```
2θ (°)     d (Å)     hkl    I/I₀
32.51      2.752     110     15
35.54      2.524     -111+002  100
38.71      2.324     111+200  96
48.72      1.867     -202     30
53.49      1.712     020     15
58.27      1.583     202     25
61.53      1.507     -113     25
66.22      1.410     022+311  15
68.12      1.375     220     15
```
peaks: [(32.51, 15), (35.54, 100), (38.71, 96), (48.72, 30), (53.49, 15), (58.27, 25), (61.53, 25), (66.22, 15), (68.12, 15)]

### Cu₂O — Cuprite, Cubic (PDF #05-0667)
```
2θ (°)     d (Å)     hkl    I/I₀
29.55      3.020     110     10
36.42      2.465     111     100
42.30      2.135     200     40
61.37      1.510     220     30
73.53      1.287     311     20
```
peaks: [(29.55, 10), (36.42, 100), (42.30, 40), (61.37, 30)]

### Cu — FCC Metal (PDF #04-0836)
```
2θ (°)     d (Å)     hkl    I/I₀
43.30      2.088     111     100
50.43      1.808     200     46
74.13      1.278     220     20
```
peaks: [(43.30, 100), (50.43, 46), (74.13, 20)]

## Common Byproducts / Impurities

### K₂SO₄ — Arcanite, Orthorhombic (PDF #05-0613)
```
2θ (°)     d (Å)     hkl    I/I₀
21.36      4.157     200     25
22.08      4.022     021     20
29.77      2.999     130     90
30.80      2.901     131     100
31.05      2.879     112     70
36.16      2.482     221     30
43.08      2.098     330     20
```
peaks: [(21.36, 25), (22.08, 20), (29.77, 90), (30.80, 100), (31.05, 70), (36.16, 30), (43.08, 20)]
**Note**: K₂SO₄ solubility = 120 g/L at 25°C. In reactions using KOH + sulfate salts,
it is typically removed by washing. Only assign K₂SO₄ if washing was confirmed insufficient.

### Na₂SO₄ — Thenardite (PDF #37-1465)
```
2θ (°)     d (Å)     hkl    I/I₀
19.01      4.664     111     35
28.09      3.174     131     20
32.10      2.786     040     30
33.46      2.676     220     25
```
peaks: [(19.01, 35), (28.09, 20), (32.10, 30), (33.46, 25)]

### KCl — Sylvite (PDF #41-1476)
```
2θ (°)     d (Å)     hkl    I/I₀
28.35      3.146     200     100
40.51      2.224     220     60
50.22      1.815     222     25
```
peaks: [(28.35, 100), (40.51, 60), (50.22, 25)]

### KNO₃ — Niter (PDF #05-0377)
```
2θ (°)     d (Å)     hkl    I/I₀
23.56      3.773     111     100
30.81      2.899     112     25
33.86      2.645     211     20
```
peaks: [(23.56, 100), (30.81, 25), (33.86, 20)]

## Usage

In Python, use the peaks lists directly:

```python
from scripts.xrd_analyzer import match_phases

cuoh2_peaks = [(16.72, 40), (23.84, 60), (33.97, 70), (35.90, 50), (39.80, 30), (53.29, 40)]
cuo_peaks = [(32.51, 15), (35.54, 100), (38.71, 96), ...]

phases = {"Cu(OH)₂ PDF#13-0420": cuoh2_peaks, "CuO PDF#48-1548": cuo_peaks}
results = match_phases(observed_peaks, phases)
```

**IMPORTANT**: Always confirm expected phases with the user. Never assign phases
from this list without user confirmation.
