#!/usr/bin/env python3
"""Repair saved Figure 6/7 display geometry without recomputing scientific fits.

The existing raster panels 6a (apart from its title), 6b, 7a, and 7b are
preserved at native resolution. The modified panels use the frozen CSV rows.
All files are written to a new output directory. No scientific inputs mutate.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import sys
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import PIL
from PIL import Image, ImageDraw, ImageFont

FONT_DIR = Path('/usr/share/fonts/truetype/dejavu')
COLORS = {'F': '#3B6FB6', 'C': '#D8902F', 'S': '#6B4C9A', 'P': '#4C956C', 'O': '#8A8A8A'}
TEXT_BOXES = []
PANEL_DATA = {}


def apply_figure_style():
    """Apply the figure-style role ladder for the existing Pillow renderer.

    Matplotlib is absent on the execution host; this is a Pillow adaptation,
    retaining the established colors and native 300-dpi raster geometry.
    """
    return {'base': 27, 'legend': 24, 'tick': 22, 'letter': 36, 'dpi': 300}


STYLE = apply_figure_style()


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path):
    with Path(path).open(newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle))


def font(role='base', bold=False):
    return ImageFont.truetype(str(FONT_DIR / ('DejaVuSans-Bold.ttf' if bold else 'DejaVuSans.ttf')), STYLE[role])


def text(image, xy, value, role='base', anchor='la', bold=False, panel='', fill='#222222'):
    draw = ImageDraw.Draw(image)
    box = draw.textbbox(xy, str(value), font=font(role, bold), anchor=anchor)
    draw.text(xy, str(value), font=font(role, bold), anchor=anchor, fill=fill)
    TEXT_BOXES.append({'figure': image.info['figure'], 'panel': panel, 'text': str(value), 'bbox': box, 'role': role})


def rotated_text(image, xy, value, panel=''):
    f = font()
    box = f.getbbox(value)
    tile = Image.new('RGBA', (box[2] - box[0] + 4, box[3] - box[1] + 4))
    ImageDraw.Draw(tile).text((2 - box[0], 2 - box[1]), value, font=f, fill='#222222')
    tile = tile.rotate(90, expand=True)
    left, top = round(xy[0] - tile.width / 2), round(xy[1] - tile.height / 2)
    image.alpha_composite(tile, (left, top))
    TEXT_BOXES.append({'figure': image.info['figure'], 'panel': panel, 'text': value,
                       'bbox': (left, top, left + tile.width, top + tile.height), 'role': 'base'})


class Panel:
    def __init__(self, image, rect, title, xlabel, ylabel, xlim, ylim, letter,
                 xticks, yticks, xlabels=None, ylabels=None):
        self.image = image
        self.draw = ImageDraw.Draw(image)
        self.rect = rect
        self.letter = letter
        l, t, r, b = rect
        self.plot = (l + 150, t + 92, r - 35, b - 120)
        self.xlim, self.ylim = xlim, ylim
        text(image, (l + 16, t + 12), letter, 'letter', bold=True, panel=letter)
        text(image, ((l + r) / 2 + 22, t + 35), title, anchor='ma', bold=True, panel=letter)
        text(image, ((self.plot[0] + self.plot[2]) / 2, b - 35), xlabel, anchor='ma', panel=letter)
        rotated_text(image, (l + 37, (self.plot[1] + self.plot[3]) / 2), ylabel, letter)
        x0, y0, x1, y1 = self.plot
        self.draw.line((x0, y1, x1, y1), fill='#222222', width=3)
        self.draw.line((x0, y0, x0, y1), fill='#222222', width=3)
        for i, v in enumerate(xticks):
            px = self.tx(v)
            self.draw.line((px, y1, px, y1 + 9), fill='#222222', width=2)
            text(image, (px, y1 + 18), xlabels[i] if xlabels else f'{v:g}', 'tick', anchor='ma', panel=letter)
        for i, v in enumerate(yticks):
            py = self.ty(v)
            self.draw.line((x0 - 9, py, x0, py), fill='#222222', width=2)
            text(image, (x0 - 15, py), ylabels[i] if ylabels else f'{v:g}', 'tick', anchor='rm', panel=letter)

    def tx(self, x):
        return self.plot[0] + (x - self.xlim[0]) / (self.xlim[1] - self.xlim[0]) * (self.plot[2] - self.plot[0])

    def ty(self, y):
        return self.plot[3] - (y - self.ylim[0]) / (self.ylim[1] - self.ylim[0]) * (self.plot[3] - self.plot[1])

    def line(self, xs, ys, color, width=4):
        self.draw.line([(self.tx(x), self.ty(y)) for x, y in zip(xs, ys)], fill=color, width=width)

    def hline(self, y, color='#BBBBBB', width=2):
        self.draw.line((self.plot[0], self.ty(y), self.plot[2], self.ty(y)), fill=color, width=width)

    def vline(self, x, color='#BBBBBB', width=2):
        self.draw.line((self.tx(x), self.plot[1], self.tx(x), self.plot[3]), fill=color, width=width)

    def marker(self, x, y, color, shape='o', radius=7):
        self.pixel_marker(self.tx(x), self.ty(y), color, shape, radius)

    def pixel_marker(self, px, py, color, shape='o', radius=7):
        if shape == 's':
            self.draw.rectangle((px-radius, py-radius, px+radius, py+radius), fill=color)
        elif shape == '^':
            self.draw.polygon([(px, py-radius), (px-radius, py+radius), (px+radius, py+radius)], fill=color)
        else:
            self.draw.ellipse((px-radius, py-radius, px+radius, py+radius), fill=color)


def figure6(inputs):
    src = Image.open(inputs / 'Figure6_original.png').convert('RGBA')
    assert src.size == (2160, 1560)
    im = src.copy()
    im.info['figure'] = 'Figure6'
    # Preserve the curves, windows, identities, and information-share marks.
    ImageDraw.Draw(im).rectangle((70, 15, 1079, 79), fill='white')
    text(im, (562, 35), 'Predeclared windows separate critical and SPR features', anchor='ma', bold=True, panel='a')
    ImageDraw.Draw(im).rectangle((0, 780, 2159, 1559), fill='white')
    equal = rows(inputs / 'V24_EQUAL_BUDGET_COMPARISON.csv')
    selected = [r for r in equal if r['selected_within_strategy'] == '1' and r['strategy'] in ('F', 'C', 'S', 'P')]
    order = ['F', 'C', 'S', 'P']
    vals = [float(next(r for r in selected if r['strategy'] == s)['local_sd_riu']) * 1e6 for s in order]
    p = Panel(im, (0, 780, 1080, 1560), 'Equal budget favors the focused SPR window',
              'Strategy', 'Local SD (micro-RIU)', (-0.5, 3.5), (0, max(vals)*1.12), 'c',
              list(range(4)), [0, 50, 100, 150], ['Full F', 'Critical C', 'SPR S', 'Bare prism P'])
    for i, (s, val) in enumerate(zip(order, vals)):
        p.line([i, i], [0, val], '#CCCCCC', 4)
        p.marker(i, val, COLORS[s], radius=10)
    checks = rows(inputs / 'V24_DESIGN_NUMERICAL_CHECKS.csv')
    phase = [r for r in checks if r['check_type'] == 'grid_phase' and r['strategy'] in order]
    series = {}
    for s in order:
        g = sorted([r for r in phase if r['strategy'] == s], key=lambda r: float(r['phase_fraction_of_spacing']))
        xs = [float(r['phase_fraction_of_spacing']) for r in g]
        ys = [float(r['local_sd_riu']) for r in g]
        assert xs == [-0.5, 0.0, 0.5]
        series[s] = {'phase': xs, 'ratio_to_central_phase': [v/ys[1] for v in ys]}
    ratios = [v for a in series.values() for v in a['ratio_to_central_phase']]
    p = Panel(im, (1080, 780, 2160, 1560), 'Critical-angle designs retain grid-phase sensitivity',
              'Grid phase (fraction of spacing)', 'SD / centre-phase SD', (-0.54, 0.54),
              (min(ratios)*.98, max(ratios)*1.02), 'd', [-0.5, -0.25, 0, 0.25, 0.5], [1.0, 1.05, 1.10, 1.14])
    p.hline(1)
    for s in order:
        xx, yy = series[s]['phase'], series[s]['ratio_to_central_phase']
        p.line(xx, yy, COLORS[s])
        for x, y in zip(xx, yy):
            p.marker(x, y, COLORS[s], radius=6)
    for i, s in enumerate(order):
        col, row = i % 2, i // 2
        x, y = p.plot[0] + 12 + col*215, p.plot[1] + 18 + row*38
        p.draw.line((x, y+10, x+35, y+10), fill=COLORS[s], width=5)
        text(im, (x+45, y+10), s, 'legend', anchor='lm', panel='d')
    # Pixel-level preservation checks apply outside the one title-edit box.
    assert src.crop((0, 80, 1080, 780)).tobytes() == im.crop((0, 80, 1080, 780)).tobytes()
    assert src.crop((1080, 0, 2160, 780)).tobytes() == im.crop((1080, 0, 2160, 780)).tobytes()
    PANEL_DATA['Figure6'] = {'panel_count': 4, 'discretization_bins': 201, 'curve_recomputed': False,
                            'panel_a_plot_preserved_pixel_exact': True, 'panel_b_preserved_pixel_exact': True,
                            'panel_c_local_sd_micro_riu': dict(zip(order, vals)), 'panel_d': series}
    return im, {'a': (0, 0, 1080, 780), 'b': (1080, 0, 2160, 780),
                'c': (0, 780, 1080, 1560), 'd': (1080, 780, 2160, 1560)}


def figure7(inputs):
    src = Image.open(inputs / 'Figure7_original.png').convert('RGBA')
    assert src.size == (2160, 900)
    im = src.copy()
    im.info['figure'] = 'Figure7'
    ImageDraw.Draw(im).rectangle((1440, 0, 2159, 899), fill='white')
    predictions = {r['anchor_id']: r for r in rows(inputs / 'V24_THRESHOLD_PREDICTIONS.csv')}
    paired = [r for r in rows(inputs / 'V24_CONFIRMATION_SUMMARY.csv') if r['summary_type'] == 'paired_D_T0_minus_T1']
    assert len(predictions) == 9 and len(paired) == 63
    points = []
    for r in paired:
        pr = predictions[r['anchor_id']]
        scale = float(pr['v1_minus_v0_riu2'])
        x, y = float(r['local_predicted_D_riu2'])/scale, float(r['mean_D_riu2'])/scale
        err = 1.96*float(r['D_mcse_riu2'])/scale
        lo, hi = y-err, y+err
        assert math.isclose(lo, float(r['D_mc_interval_low_riu2'])/scale, abs_tol=1e-12)
        assert math.isclose(hi, float(r['D_mc_interval_high_riu2'])/scale, abs_tol=1e-12)
        assert int(r['n']) == 200
        points.append({'anchor_id': r['anchor_id'], 'rho': float(r['rho_delta_over_pred']),
                       'x': x, 'y': y, 'low': lo, 'high': hi, 'n': int(r['n']),
                       'thickness_nm': int(float(pr['thickness_nm'])), 'photons': float(pr['total_incident_photons'])})
    p = Panel(im, (1440, 0, 2160, 900), 'Paired inversion tests local risk',
              'Local predicted D / (v1-v0)', 'Observed paired D / (v1-v0)',
              (-1.3, 3.5), (-4.5, 9.6), 'c', [-1, 0, 1, 2, 3], [-4, -2, 0, 2, 4, 6, 8])
    p.line([-1.3, 3.5], [-1.3, 3.5], '#777777', 3)
    p.hline(0)
    p.vline(0)
    shapes = {45: 'o', 50: 's', 55: '^'}
    colors = {1e5: '#9EC1E6', 1e6: '#3B6FB6', 1e7: '#163A63'}
    for r in points:
        assert p.ylim[0] < r['low'] <= r['y'] <= r['high'] < p.ylim[1]
        x, y0, y1 = p.tx(r['x']), p.ty(r['low']), p.ty(r['high'])
        assert p.plot[0] + 7 < x < p.plot[2] - 7
        assert p.plot[1] + 7 < y1 < y0 < p.plot[3] - 7
        p.draw.line((x, y0, x, y1), fill=colors[r['photons']], width=2)
    # Draw all 63 means after the intervals so an interval cannot cover a mean.
    for r in points:
        p.marker(r['x'], r['y'], colors[r['photons']], shapes[r['thickness_nm']], 6)
    # Two keys are adjacent to panel c and occupy data-free upper whitespace.
    text(im, (p.plot[0]+16, p.plot[1]+17), 'Thickness', 'legend', panel='c')
    text(im, (p.plot[0]+221, p.plot[1]+17), 'Photons', 'legend', panel='c')
    for i, (d, n) in enumerate(zip([45, 50, 55], [1e5, 1e6, 1e7])):
        yy = p.plot[1] + 65 + i*37
        p.pixel_marker(p.plot[0]+27, yy+9, '#444444', shapes[d], 7)
        text(im, (p.plot[0]+45, yy+9), f'{d} nm', 'legend', anchor='lm', panel='c')
        p.draw.line((p.plot[0]+224, yy+9, p.plot[0]+253, yy+9), fill=colors[n], width=5)
        text(im, (p.plot[0]+264, yy+9), {1e5: '10\u2075', 1e6: '10\u2076', 1e7: '10\u2077'}[n], 'legend', anchor='lm', panel='c')
    assert src.crop((0, 0, 1440, 900)).tobytes() == im.crop((0, 0, 1440, 900)).tobytes()
    PANEL_DATA['Figure7'] = {'panel_count': 3, 'panel_a_b_preserved_pixel_exact': True,
        'panel_c_x_quantity': 'local_predicted_D_riu2 / v1_minus_v0_riu2', 'x_jitter_applied': False,
        'xlim': p.xlim, 'ylim': p.ylim, 'condition_count': len(points), 'replicates_per_condition': 200,
        'mean_markers_drawn': len(points), 'complete_intervals_drawn': len(points),
        'clipped_intervals': 0, 'omitted_mean_markers': 0,
        'original_clipped_intervals': sum(r['low'] < -2.2 or r['high'] > 5.2 for r in points),
        'original_omitted_mean_markers': sum(not -2.2 <= r['y'] <= 5.2 for r in points),
        'interval_range': [min(r['low'] for r in points), max(r['high'] for r in points)], 'points': points}
    return im, {'a': (0, 0, 720, 900), 'b': (720, 0, 1440, 900), 'c': (1440, 0, 2160, 900)}


def bbox_audit(images):
    overlaps, outside = [], []
    for t in TEXT_BOXES:
        w, h = images[t['figure']].size
        x0, y0, x1, y1 = t['bbox']
        if x0 < 0 or y0 < 0 or x1 > w or y1 > h:
            outside.append(t)
    for a, b in combinations(TEXT_BOXES, 2):
        if a['figure'] != b['figure']:
            continue
        aa, bb = a['bbox'], b['bbox']
        if min(aa[2], bb[2]) > max(aa[0], bb[0]) and min(aa[3], bb[3]) > max(aa[1], bb[1]):
            overlaps.append({'a': a, 'b': b})
    return {'scope': 'all newly rendered text boxes; unchanged native panels reviewed visually',
            'text_count': len(TEXT_BOXES), 'outside_canvas': outside, 'overlapping_text_boxes': overlaps,
            'boxes': TEXT_BOXES, 'pass': not overlaps and not outside}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input-dir', type=Path, required=True)
    ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args()
    start = now()
    if args.output_dir.exists():
        raise FileExistsError(f'Refusing to overwrite {args.output_dir}')
    args.output_dir.mkdir(parents=True)
    input_hashes = {str(p.resolve()): sha(p) for p in sorted(args.input_dir.iterdir()) if p.is_file()}
    rendered = {}
    for name, builder in [('Figure6', figure6), ('Figure7', figure7)]:
        im, boxes = builder(args.input_dir)
        rendered[name] = im
        im.convert('RGB').save(args.output_dir/f'{name}_display_fixed.png', dpi=(300, 300), optimize=True)
        im.convert('RGB').save(args.output_dir/f'{name}_display_fixed.pdf', 'PDF', resolution=300.0)
        cropdir = args.output_dir/'panel_crops'
        cropdir.mkdir(exist_ok=True)
        for letter, box in boxes.items():
            im.crop(box).convert('RGB').save(cropdir/f'{name}_{letter}.png', dpi=(300, 300))
    audit = bbox_audit(rendered)
    (args.output_dir/'text_bbox_audit.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
    (args.output_dir/'plot_data_audit.json').write_text(json.dumps(PANEL_DATA, indent=2), encoding='utf-8')
    assert input_hashes == {str(p.resolve()): sha(p) for p in sorted(args.input_dir.iterdir()) if p.is_file()}
    manifest = {'start_utc': start, 'end_utc': now(), 'actual_command': ' '.join(sys.argv),
                'platform': platform.platform(), 'python': sys.version, 'pillow': PIL.__version__,
                'renderer': 'Pillow; unchanged source raster panels retained at native resolution',
                'style': STYLE, 'random_seed': 'NA: no stochastic computation',
                'cpu_threads': {k: os.environ.get(k) for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']},
                'input_hashes': input_hashes, 'input_hashes_unchanged': True,
                'outputs': {str(p.resolve()): sha(p) for p in sorted(args.output_dir.rglob('*')) if p.is_file()},
                'scientific_fit_recomputed': False, 'csv_values_modified': False,
                'bbox_check_pass': audit['pass'], 'exit_status': 0 if audit['pass'] else 2}
    (args.output_dir/'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in manifest.items() if k not in ['outputs', 'input_hashes']}, indent=2))
    print('INPUT_COUNT', len(input_hashes), 'OUTPUT_COUNT', len(manifest['outputs']))
    if not audit['pass']:
        print(json.dumps({'overlaps': audit['overlapping_text_boxes'], 'outside': audit['outside_canvas']}, indent=2))
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
