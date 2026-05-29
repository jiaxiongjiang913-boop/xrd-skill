"""
XRD Web Analyzer — Flask + matplotlib backend.
Renders exact Nature-style figure per xrd-skill spec.
"""
import re, io, os, base64, uuid
import numpy as np
from numpy.polynomial import polynomial as P
from flask import Flask, request, jsonify, render_template, send_file

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

app = Flask(__name__)
UPLOAD_DIR = os.path.join(os.path.dirname(__file__), 'uploads')
os.makedirs(UPLOAD_DIR, exist_ok=True)
data_store = {}

# ============================================================
# FONT
# ============================================================
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Arial'
rcParams['mathtext.it'] = 'Arial:italic'
rcParams['font.size'] = 7

# ============================================================
# PHASES
# ============================================================
PHASES = {
    "Cu(OH)₂ PDF#13-0420": [(16.72,40),(23.84,60),(33.97,70),(35.90,50),(39.80,30),(53.29,40)],
    "CuO PDF#48-1548": [(32.51,15),(35.54,100),(38.71,96),(48.72,30),(53.49,15),(58.27,25),(61.53,25)],
    "Cu₂O PDF#05-0667": [(29.55,10),(36.42,100),(42.30,40),(61.37,30)],
    "Cu PDF#04-0836": [(43.30,100),(50.43,46),(74.13,20)],
    "K₂SO₄ PDF#05-0613": [(21.36,25),(22.08,20),(29.77,90),(30.80,100),(31.05,70),(36.16,30)],
    "KCl PDF#41-1476": [(28.35,100),(40.51,60),(50.22,25)],
    "α-MoO₃ PDF#05-0508": [(12.76,31),(23.33,43),(25.70,100),(27.33,90),(33.72,20),(35.40,14),(38.97,27),(39.63,14),(42.90,5),(45.75,8),(46.28,14),(49.33,6),(52.80,7),(55.22,10),(56.33,8),(58.84,6)],
    "Sn₂S₃ PDF#27-0899": [(12.635,90),(14.802,60),(21.136,25),(21.498,40),(23.771,30),(25.427,18),(29.756,12),(30.484,100),(30.807,15),(32.902,70),(33.536,20),(34.466,12),(35.891,10),(37.12,8),(38.1,17),(39.311,11)],
    "SnS PDF#39-0354": [(21.997,30),(22.298,10),(26.001,25),(27.46,20),(30.469,35),(31.523,100),(31.959,50),(32.863,8),(34.512,12),(39.046,16)],
}
PCOLORS = ['#457b9d','#e63946','#2a9d8f','#e76f51','#264653','#f4a261','#a8dadc','#6d597a']

# ============================================================
# PARSING
# ============================================================
def parse_rigaku_raw(filepath):
    metadata = {}
    start, stop, step, count = None, None, None, None
    intensities = []
    data_section = False
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line == '*BEGIN': continue
            if line.startswith('*') and not data_section:
                m = re.match(r'\*(.+?)\s*=\s*(.+)', line)
                if m:
                    k, v = m.group(1).strip(), m.group(2).strip()
                    metadata[k] = v
                    if k == 'START': start = float(v)
                    elif k == 'STOP': stop = float(v)
                    elif k == 'STEP': step = float(v)
                    elif k == 'COUNT': count = int(v); data_section = True
            elif data_section:
                for p in line.replace(',', ' ').split():
                    try: intensities.append(float(p))
                    except ValueError: pass
    tt = np.linspace(start, stop, count)
    raw = np.array(intensities[:count])
    return metadata, tt, raw

def parse_two_column(filepath):
    xy = []
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#') or s.startswith('*'): continue
            parts = s.split()
            if len(parts) >= 2:
                try: xy.append([float(parts[0]), float(parts[1])])
                except ValueError: continue
    return {}, np.array([p[0] for p in xy]), np.array([p[1] for p in xy])

def auto_parse(filepath):
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        head = ''.join(f.readline() for _ in range(3))
    if '*TYPE' in head: return parse_rigaku_raw(filepath)
    return parse_two_column(filepath)

# ============================================================
# ANALYSIS
# ============================================================
def subtract_background(tt, raw, poly_deg=4, n_iter=3, sigma_thresh=2.5):
    mask = np.ones(len(raw), dtype=bool)
    for _ in range(n_iter):
        try:
            c = P.polyfit(tt[mask], raw[mask], deg=poly_deg)
            bg = P.polyval(tt, c)
            res = raw - bg
            sigma = np.std(res[mask])
            if sigma == 0: break
            mask = res < sigma_thresh * sigma
        except np.linalg.LinAlgError: break
    c = P.polyfit(tt[mask], raw[mask], deg=poly_deg)
    bg_final = P.polyval(tt, c)
    return bg_final, raw - bg_final

def find_peaks(tt, signal, raw_intensity, sigma_thresh=2.0, min_snr=1.5):
    noise = np.std(signal)
    threshold = sigma_thresh * noise
    peaks = []
    i = 1
    n = len(signal)
    while i < n - 1:
        if signal[i] > threshold and signal[i] > signal[i-1] and signal[i] >= signal[i+1]:
            snr = signal[i] / noise if noise > 0 else 0
            if snr >= min_snr:
                half = signal[i] / 2
                left = i
                while left > 0 and signal[left] > half: left -= 1
                right = i
                while right < n - 1 and signal[right] > half: right += 1
                fwhm = tt[right] - tt[left]
                d = 1.54059 / (2 * np.sin(np.radians(tt[i] / 2)))
                peaks.append(dict(pos=round(float(tt[i]),3), d=round(float(d),4),
                    intensity=float(raw_intensity[i]), snr=round(float(snr),1), fwhm=round(float(fwhm),3)))
            i = right + 1
        else: i += 1
    peaks.sort(key=lambda p: -p['snr'])
    return peaks

def match_phase(obs, ref, tol=0.3):
    matched, details = 0, []
    for rp, _ in ref:
        best, best_d = None, tol + 1
        for p in obs:
            d = abs(p['pos'] - rp)
            if d < tol and d < best_d: best_d = d; best = p
        if best: matched += 1; details.append(dict(rp=rp, op=best['pos'], delta=round(best_d,3)))
    return matched, len(ref), details

# ============================================================
# FIGURE — IRON RULE spec from SKILL.md (17 rules, enforced)
# ============================================================
def build_figure(tt, raw_norm, filename, selected_phases, peaks, params):
    """Build XRD figure strictly per IRON RULE.
    Rules 1-3: 2 subplots sharex hspace=0, 0.62/0.38 ratio, axis box frames (NO outer rect)
    Rules 4-7: top panel — no y ticks, 4 spines, y-label, sample label
    Rules 8-11: bottom panel — 4 spines, ylim(0,card_h), no y ticks, card_h=54
    Rules 12-15: stick height/2, rooted at floor, divider for multi-phase, phase labels
    Rules 16-17: Arial 5-7pt, mathtext formulas
    """
    card_h = params.get('cardHeight', 54)
    stick_s = params.get('stickScale', 0.5)

    # ── rcParams per IRON RULE ──
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'svg.fonttype': 'none',
        'font.size': 7,
        'axes.linewidth': 0.4,
        'xtick.major.width': 0.4,
        'ytick.major.width': 0.4,
    })

    fig = plt.figure(figsize=(6.3, 2.8))
    gs = fig.add_gridspec(2, 1, hspace=0, height_ratios=[0.62, 0.38])
    ax_top = fig.add_subplot(gs[0])
    ax_bot = fig.add_subplot(gs[1], sharex=ax_top)

    xmin, xmax = tt[0], tt[-1]

    # ===== TOP PANEL — XRD Curve (rules 4-7) =====
    ax_top.plot(tt, raw_norm, color='#1a1a2e', lw=1.0, clip_on=True)
    ax_top.set_xlim(xmin, xmax)
    ax_top.set_ylim(-2, 108)
    # Rule 4: no y ticks/labels
    ax_top.tick_params(axis='y', left=False, labelleft=False)
    ax_top.tick_params(axis='x', bottom=False, labelbottom=False)
    # Rule 5: all 4 spines visible = box frame
    for sp in ax_top.spines.values():
        sp.set_visible(True)
        sp.set_linewidth(0.4)
        sp.set_color('#333333')
    # Rule 6: y-label
    ax_top.set_ylabel('Intensity (a.u.)', fontsize=7, color='#333333', labelpad=1)
    # Rule 7: sample label upper-left
    ax_top.text(0.02, 0.92, filename.replace('_', ' '),
                transform=ax_top.transAxes, fontsize=7, va='top', color='#1a1a2e')

    # ===== BOTTOM PANEL — Card Box (rules 8-11) =====
    # Rule 9+11: ylim bottom=0, top=card_h
    ax_bot.set_ylim(0, card_h)
    ax_bot.set_xlim(xmin, xmax)
    # Rule 8: all 4 spines visible = box frame
    for sp in ax_bot.spines.values():
        sp.set_visible(True)
        sp.set_linewidth(0.4)
        sp.set_color('#333333')
    # Rule 10: no y ticks/labels
    ax_bot.tick_params(axis='y', left=False, labelleft=False)
    ax_bot.set_xlabel(r'2$\theta$ (deg)', fontsize=7, color='#333333', labelpad=1)
    ax_bot.tick_params(axis='x', which='major', labelsize=6, pad=1)

    # ===== Sticks (rules 12-15) =====
    np_h = len(selected_phases)
    if np_h > 0:
        sec_h = card_h / np_h
        for idx, (pname, ref) in enumerate(selected_phases):
            clr = PCOLORS[idx % len(PCOLORS)]
            base_y = idx * sec_h  # rooted at floor or divider
            for rp, ri in ref:
                stick_h = (ri / 100.0) * sec_h * stick_s  # Rule 12: heights / 2 via stick_s
                ax_bot.plot([rp, rp], [base_y, base_y + stick_h],
                            color=clr, lw=0.8, solid_capstyle='butt', clip_on=True)
            # Rule 15: phase label at top of section
            ax_bot.text(xmin + 0.5, base_y + sec_h - 2, pname,
                        fontsize=6, va='top', ha='left', color=clr)
            # Rule 14: divider line between phases
            if idx > 0:
                ax_bot.plot([xmin, xmax], [base_y, base_y],
                            color='#bbbbbb', lw=0.5, zorder=1)

    # Rule 3: NO outer rectangle — axis spines ARE the frame

    # ── Save PNG for preview ──
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight',
                facecolor='white', pad_inches=0.04)
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()


def build_figure_svg_pdf(tt, raw_norm, filename, selected_phases, peaks, params):
    """Same as build_figure but returns dict with SVG and PDF base64."""
    card_h = params.get('cardHeight', 54)
    stick_s = params.get('stickScale', 0.5)

    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'svg.fonttype': 'none',
        'font.size': 7,
        'axes.linewidth': 0.4,
        'xtick.major.width': 0.4,
        'ytick.major.width': 0.4,
    })

    fig = plt.figure(figsize=(6.3, 2.8))
    gs = fig.add_gridspec(2, 1, hspace=0, height_ratios=[0.62, 0.38])
    ax_top = fig.add_subplot(gs[0])
    ax_bot = fig.add_subplot(gs[1], sharex=ax_top)

    xmin, xmax = tt[0], tt[-1]

    ax_top.plot(tt, raw_norm, color='#1a1a2e', lw=1.0, clip_on=True)
    ax_top.set_xlim(xmin, xmax); ax_top.set_ylim(-2, 108)
    ax_top.tick_params(axis='y', left=False, labelleft=False)
    ax_top.tick_params(axis='x', bottom=False, labelbottom=False)
    for sp in ax_top.spines.values():
        sp.set_visible(True); sp.set_linewidth(0.4); sp.set_color('#333333')
    ax_top.set_ylabel('Intensity (a.u.)', fontsize=7, color='#333333', labelpad=1)
    ax_top.text(0.02, 0.92, filename.replace('_', ' '),
                transform=ax_top.transAxes, fontsize=7, va='top', color='#1a1a2e')

    ax_bot.set_ylim(0, card_h); ax_bot.set_xlim(xmin, xmax)
    for sp in ax_bot.spines.values():
        sp.set_visible(True); sp.set_linewidth(0.4); sp.set_color('#333333')
    ax_bot.tick_params(axis='y', left=False, labelleft=False)
    ax_bot.set_xlabel(r'2$\theta$ (deg)', fontsize=7, color='#333333', labelpad=1)
    ax_bot.tick_params(axis='x', which='major', labelsize=6, pad=1)

    np_h = len(selected_phases)
    if np_h > 0:
        sec_h = card_h / np_h
        for idx, (pname, ref) in enumerate(selected_phases):
            clr = PCOLORS[idx % len(PCOLORS)]
            base_y = idx * sec_h
            for rp, ri in ref:
                stick_h = (ri / 100.0) * sec_h * stick_s
                ax_bot.plot([rp, rp], [base_y, base_y + stick_h],
                            color=clr, lw=0.8, solid_capstyle='butt', clip_on=True)
            ax_bot.text(xmin + 0.5, base_y + sec_h - 2, pname,
                        fontsize=6, va='top', ha='left', color=clr)
            if idx > 0:
                ax_bot.plot([xmin, xmax], [base_y, base_y],
                            color='#bbbbbb', lw=0.5, zorder=1)

    fig.tight_layout(pad=0.5)
    buf_png = io.BytesIO()
    fig.savefig(buf_png, format='png', dpi=300, bbox_inches='tight', facecolor='white')
    buf_png.seek(0)
    buf_svg = io.BytesIO()
    fig.savefig(buf_svg, format='svg', bbox_inches='tight', facecolor='white')
    buf_svg.seek(0)
    buf_pdf = io.BytesIO()
    fig.savefig(buf_pdf, format='pdf', bbox_inches='tight', facecolor='white')
    buf_pdf.seek(0)
    plt.close(fig)
    return {
        'png': base64.b64encode(buf_png.read()).decode(),
        'svg': base64.b64encode(buf_svg.read()).decode(),
        'pdf': base64.b64encode(buf_pdf.read()).decode(),
    }

# ============================================================
# ROUTES
# ============================================================
@app.route('/')
def index():
    return render_template('index.html', phases=list(PHASES.keys()))

@app.route('/api/phases')
def api_phases():
    return jsonify(list(PHASES.keys()))

@app.route('/api/upload', methods=['POST'])
def api_upload():
    if 'file' not in request.files:
        return jsonify({'error': 'No file'}), 400
    file = request.files['file']
    if file.filename == '': return jsonify({'error': 'Empty'}), 400
    fpath = os.path.join(UPLOAD_DIR, file.filename)
    file.save(fpath)
    try:
        meta, tt, raw = auto_parse(fpath)
    except Exception as e:
        return jsonify({'error': str(e)}), 400
    sid = str(uuid.uuid4())
    data_store[sid] = dict(filepath=fpath, meta=meta, tt=tt.tolist(), raw=raw.tolist(), fname=file.filename)
    return jsonify(dict(session_id=sid, fname=file.filename,
        scan_range=f"{tt[0]:.1f}°–{tt[-1]:.1f}°",
        n_points=len(tt), i_range=f"{raw.min():.0f}–{raw.max():.0f}"))

@app.route('/api/analyze', methods=['POST'])
def api_analyze():
    data = request.get_json()
    sid = data.get('session_id')
    params = data.get('params', {})
    if sid not in data_store: return jsonify({'error': 'Session lost'}), 400
    stored = data_store[sid]
    tt = np.array(stored['tt'])
    raw = np.array(stored['raw'])

    bg = params.get('bgSubtract', False)
    pd_ = int(params.get('polyDeg', 4))
    st = float(params.get('sigmaThresh', 2.0))
    ms = float(params.get('minSnr', 1.5))

    if bg:
        _, sig = subtract_background(tt, raw, poly_deg=pd_)
    else:
        sig = raw.copy()

    peaks = find_peaks(tt, sig, raw, sigma_thresh=st, min_snr=ms)

    sel_phases = []
    pm = []
    for pn in params.get('phases', []):
        if pn in PHASES:
            ref = PHASES[pn]
            sel_phases.append((pn, ref))
            m, t, d = match_phase(peaks, ref)
            pm.append(dict(name=pn, matched=m, total=t, ratio=round(m/max(t,1)*100), details=d))

    raw_norm = np.array(sig, dtype=float)
    if raw_norm.max() > 0: raw_norm = raw_norm / raw_norm.max()

    img_b64 = build_figure(tt, raw_norm, stored['fname'].rsplit('.',1)[0],
                           sel_phases, peaks, params)

    # Preview PNG for the interactive UI
    return jsonify(dict(image=img_b64, peaks=peaks[:15], n_peaks=len(peaks), phase_matches=pm))


@app.route('/api/download/<fmt>', methods=['POST'])
def api_download(fmt):
    """Download SVG or PDF of the current figure."""
    data = request.get_json()
    sid = data.get('session_id')
    params = data.get('params', {})
    if sid not in data_store: return jsonify({'error': 'Session lost'}), 400
    stored = data_store[sid]
    tt = np.array(stored['tt'])
    raw = np.array(stored['raw'])

    bg = params.get('bgSubtract', False)
    pd_ = int(params.get('polyDeg', 4))
    if bg:
        _, sig = subtract_background(tt, raw, poly_deg=pd_)
    else:
        sig = raw.copy()

    peaks = find_peaks(tt, sig, raw,
        sigma_thresh=float(params.get('sigmaThresh', 2.0)),
        min_snr=float(params.get('minSnr', 1.5)))

    sel_phases = []
    for pn in params.get('phases', []):
        if pn in PHASES:
            sel_phases.append((pn, PHASES[pn]))

    raw_norm = np.array(sig, dtype=float)
    if raw_norm.max() > 0: raw_norm = raw_norm / raw_norm.max()

    result = build_figure_svg_pdf(tt, raw_norm,
        stored['fname'].rsplit('.', 1)[0], sel_phases, peaks, params)

    if fmt not in result:
        return jsonify({'error': f'Unknown format: {fmt}'}), 400

    buf = io.BytesIO(base64.b64decode(result[fmt]))
    buf.seek(0)
    mime = {'svg': 'image/svg+xml', 'pdf': 'application/pdf', 'png': 'image/png'}[fmt]
    return send_file(buf, mimetype=mime,
        download_name=f"xrd_{stored['fname'].rsplit('.',1)[0]}.{fmt}")


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
