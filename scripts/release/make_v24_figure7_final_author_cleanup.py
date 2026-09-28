#!/usr/bin/env python3
"""Render the two CalSPR V2.4 figures on Kylin without matplotlib."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import shutil
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from calspr.v24_design import ResponseScenario, angle_centres, gold_reflectance, physical_windows

COLORS = {"F": "#3B6FB6", "C": "#D8902F", "S": "#6B4C9A", "P": "#4C956C", "O": "#8A8A8A"}
FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"


def rows(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def font(size, bold=False):
    path = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else FONT
    return ImageFont.truetype(path, size=size)


def text(draw, xy, value, size=26, fill="#222222", anchor="la", bold=False):
    draw.text(xy, str(value), font=font(size, bold), fill=fill, anchor=anchor)


def rotated_text(image, xy, value, size=28):
    f = font(size)
    box = f.getbbox(value)
    tile = Image.new("RGBA", (box[2] - box[0] + 10, box[3] - box[1] + 10), (255, 255, 255, 0))
    ImageDraw.Draw(tile).text((5, 5), value, font=f, fill="#222222")
    tile = tile.rotate(90, expand=True)
    image.alpha_composite(tile, (int(xy[0] - tile.width / 2), int(xy[1] - tile.height / 2)))


class Panel:
    def __init__(self, image, rect, title, xlabel, ylabel, xlim, ylim, xlog=False, ylog=False, letter="a", show_xticks=True, show_yticks=True):
        self.image = image
        self.draw = ImageDraw.Draw(image)
        self.rect = rect
        l, t, r, b = rect
        self.plot = (l + 150, t + 92, r - 35, b - 120)
        self.xlim, self.ylim, self.xlog, self.ylog = xlim, ylim, xlog, ylog
        self.show_xticks, self.show_yticks = show_xticks, show_yticks
        text(self.draw, (l + 16, t + 12), letter, 36, bold=True)
        text(self.draw, ((l + r) / 2 + 22, t + 35), title, 26 if len(title) < 48 else 23, anchor="ma", bold=True)
        text(self.draw, ((self.plot[0] + self.plot[2]) / 2, b - 35), xlabel, 27, anchor="ma")
        rotated_text(image, (l + 37, (self.plot[1] + self.plot[3]) / 2), ylabel, 27)
        self.axes()

    def tx(self, x):
        lo, hi = self.xlim
        if self.xlog:
            x, lo, hi = math.log10(x), math.log10(lo), math.log10(hi)
        return self.plot[0] + (x - lo) / (hi - lo) * (self.plot[2] - self.plot[0])

    def ty(self, y):
        lo, hi = self.ylim
        if self.ylog:
            y, lo, hi = math.log10(y), math.log10(lo), math.log10(hi)
        return self.plot[3] - (y - lo) / (hi - lo) * (self.plot[3] - self.plot[1])

    def axes(self):
        d = self.draw
        x0, y0, x1, y1 = self.plot
        d.line((x0, y1, x1, y1), fill="#222222", width=3)
        d.line((x0, y0, x0, y1), fill="#222222", width=3)
        if self.xlog:
            ticks = [self.xlim[0], *[10 ** k for k in range(math.ceil(math.log10(self.xlim[0])), math.floor(math.log10(self.xlim[1])) + 1)], self.xlim[1]]
            ticks = sorted(set(ticks))
        else:
            ticks = np.linspace(self.xlim[0], self.xlim[1], 5)
        for v in ticks if self.show_xticks else []:
            px = self.tx(float(v)); d.line((px, y1, px, y1 + 9), fill="#222222", width=2)
            lab = f"1e{int(round(math.log10(v)))}" if self.xlog and abs(math.log10(v)-round(math.log10(v))) < 1e-8 else f"{v:.3g}"
            text(d, (px, y1 + 18), lab, 22, anchor="ma")
        if self.ylog:
            p0, p1 = math.floor(math.log10(self.ylim[0])), math.ceil(math.log10(self.ylim[1]))
            ticks = [self.ylim[0], *[10 ** k for k in range(p0, p1 + 1) if self.ylim[0] <= 10 ** k <= self.ylim[1]], self.ylim[1]]
            ticks = sorted(set(ticks))
        else:
            ticks = np.linspace(self.ylim[0], self.ylim[1], 5)
        for v in ticks if self.show_yticks else []:
            py = self.ty(float(v)); d.line((x0 - 9, py, x0, py), fill="#222222", width=2)
            lab = f"{v:.2g}" if abs(v) < 1 or abs(v) >= 1000 else f"{v:g}"
            text(d, (x0 - 15, py), lab, 22, anchor="rm")

    def line(self, xs, ys, color, width=5):
        pts = [(self.tx(float(x)), self.ty(float(y))) for x, y in zip(xs, ys)]
        self.draw.line(pts, fill=color, width=width, joint="curve")

    def marker(self, x, y, color, shape="o", radius=7):
        px, py = self.tx(float(x)), self.ty(float(y))
        if shape == "s":
            self.draw.rectangle((px-radius, py-radius, px+radius, py+radius), fill=color)
        elif shape == "^":
            self.draw.polygon([(px, py-radius), (px-radius, py+radius), (px+radius, py+radius)], fill=color)
        else:
            self.draw.ellipse((px-radius, py-radius, px+radius, py+radius), fill=color)

    def points(self, xs, ys, color, shape="o", radius=7):
        for x, y in zip(xs, ys):
            self.marker(x, y, color, shape, radius)

    def hline(self, y, color="#BBBBBB", width=2):
        self.draw.line((self.plot[0], self.ty(y), self.plot[2], self.ty(y)), fill=color, width=width)

    def vline(self, x, color="#BBBBBB", width=2):
        self.draw.line((self.tx(x), self.plot[1], self.tx(x), self.plot[3]), fill=color, width=width)

    def shade_x(self, lo, hi, color):
        self.draw.rectangle((self.tx(lo), self.plot[1], self.tx(hi), self.plot[3]), fill=color)

    def errorbar(self, x, y, err, color, shape):
        px, py, p0, p1 = self.tx(x), self.ty(y), self.ty(y-err), self.ty(y+err)
        p0 = min(self.plot[3], max(self.plot[1], p0)); p1 = min(self.plot[3], max(self.plot[1], p1))
        self.draw.line((px, p0, px, p1), fill=color, width=2)
        if self.ylim[0] <= y <= self.ylim[1]: self.marker(x, y, color, shape, 5)

    def legend(self, items, x=None, y=None, columns=1):
        x = self.plot[0] + 12 if x is None else x
        y = self.plot[1] + 12 if y is None else y
        for i, (label, color) in enumerate(items):
            col, row = i % columns, i // columns
            xx, yy = x + col * 215, y + row * 34
            self.draw.line((xx, yy+11, xx+35, yy+11), fill=color, width=5)
            text(self.draw, (xx+45, yy+11), label, 21, anchor="lm")


def save(image, out, stem):
    png, pdf = out / f"{stem}.png", out / f"{stem}.pdf"
    image.convert("RGB").save(png, dpi=(300, 300), optimize=True)
    image.convert("RGB").save(pdf, "PDF", resolution=300.0)
    return png, pdf


def copy_source(path, target):
    if target.exists():
        raise FileExistsError(target)
    shutil.copy2(path, target)


def figure1(cfg, scan, out, src):
    attrib, equal = rows(scan/"V24_INFORMATION_ATTRIBUTION.csv"), rows(scan/"V24_EQUAL_BUDGET_COMPARISON.csv")
    checks = rows(scan/"V24_DESIGN_NUMERICAL_CHECKS.csv")
    model_cfg = {"angle_min_deg":60.0,"angle_max_deg":80.0,"angle_bins":201,"baseline_index":1.333,"sample_delta_index":0.005,"wavelength_nm":633.0,"window_definition":{"nominal_thickness_nm":50.0,"critical_half_width_deg":1.0,"spr_half_width_deg":2.0}}
    angles = angle_centres(model_cfg); scenario = ResponseScenario("top_hat_0p10deg_main", "top_hat", 0.1, 16)
    gold = ROOT/cfg["material"]; physical_windows(model_cfg, angles, gold)
    rr = gold_reflectance(angles,633.0,50.0,1.333,gold,scenario); rs = gold_reflectance(angles,633.0,50.0,1.338,gold,scenario)
    im = Image.new("RGBA", (2160,1560), "white")
    rects=[(0,0,1080,780),(1080,0,2160,780),(0,780,1080,1560),(1080,780,2160,1560)]
    p=Panel(im,rects[0],"Frozen windows separate critical and SPR features","Incident angle (degrees)","Reflectance",(60,80),(0,1.05),letter="a")
    p.shade_x(60.7,62.6,"#F8E9D2"); p.shade_x(71,75,"#EDE7F4"); p.line(angles,rr,"#777777",4); p.line(angles,rs,COLORS["F"],5); p.legend([("Reference q=1.333","#777777"),("Sample q=1.338",COLORS["F"])])
    full=float(next(x for x in attrib if x["design"]=="F")["effective_information"]); contrib=[]
    for name in ("C","S","O"):
        row=next(x for x in attrib if x["design"]==f"F_projected_score_contribution_{name}"); contrib.append(float(row["effective_information"])/full*100)
    p=Panel(im,rects[1],"SPR carries most projected full-design information","Share of projected F information (%)","Window",(0,max(contrib)*1.15),(-0.6,2.6),letter="b",show_yticks=False)
    for i,(lab,val,col) in enumerate(zip(("Critical C","SPR S","Other O"),contrib,(COLORS["C"],COLORS["S"],COLORS["O"]))):
        p.draw.rectangle((p.tx(0),p.ty(i+0.27),p.tx(val),p.ty(i-0.27)),fill=col); text(p.draw,(p.plot[0]-15,p.ty(i)),lab,22,anchor="rm")
    selected=[x for x in equal if x["selected_within_strategy"]=="1" and x["strategy"] in ("F","C","S","P")]; order=["F","C","S","P"]
    vals=[float(next(x for x in selected if x["strategy"]==s)["local_sd_riu"])*1e6 for s in order]
    p=Panel(im,rects[2],"Equal budget favors the focused SPR window","Strategy","Local SD (micro-RIU)",(-0.5,3.5),(0,max(vals)*1.12),letter="c",show_xticks=False)
    for i,(s,v) in enumerate(zip(order,vals)):
        p.draw.line((p.tx(i),p.ty(0),p.tx(i),p.ty(v)),fill="#CCCCCC",width=4); p.marker(i,v,COLORS[s],"o",10); text(p.draw,(p.tx(i),p.plot[3]+18),("Full F","Critical C","SPR S","Bare prism P")[i],20,anchor="ma")
    phase=[x for x in checks if x["check_type"]=="grid_phase" and x["strategy"] in order]
    ratios=[]
    for s in order:
        g=sorted([x for x in phase if x["strategy"]==s],key=lambda x:float(x["phase_fraction_of_spacing"])); yy=np.array([float(x["local_sd_riu"]) for x in g]); ratios.extend((yy/yy[1]).tolist())
    p=Panel(im,rects[3],"Critical-angle designs retain grid-phase sensitivity","Grid phase (fraction of spacing)","SD / centre-phase SD",(-0.5,0.5),(min(ratios)*.98,max(ratios)*1.02),letter="d"); p.hline(1)
    for s in order:
        g=sorted([x for x in phase if x["strategy"]==s],key=lambda x:float(x["phase_fraction_of_spacing"])); xx=[float(x["phase_fraction_of_spacing"]) for x in g]; yy=np.array([float(x["local_sd_riu"]) for x in g]); yy=yy/yy[1]; p.line(xx,yy,COLORS[s],4); p.points(xx,yy,COLORS[s],radius=6)
    p.legend([(s,COLORS[s]) for s in order],columns=2)
    paths=save(im,out,"CalSPR_V24_Figure1_information_sources")
    for source in (scan/"V24_INFORMATION_ATTRIBUTION.csv",scan/"V24_EQUAL_BUDGET_COMPARISON.csv",scan/"V24_DESIGN_NUMERICAL_CHECKS.csv"): copy_source(source,src/source.name)
    return paths


def figure2(cfg,scan,confirm,out,src):
    pred=rows(scan/"V24_THRESHOLD_PREDICTIONS.csv"); summaries=rows(confirm/"V24_CONFIRMATION_SUMMARY.csv"); paired=[x for x in summaries if x["summary_type"]=="paired_D_T0_minus_T1"]; pmap={x["anchor_id"]:x for x in pred}
    im=Image.new("RGBA",(2160,900),"white"); rects=[(0,0,720,900),(720,0,1440,900),(1440,0,2160,900)]
    ys=[float(x["delta_star_schur_nm"]) for x in pred]
    p=Panel(im,rects[0],"Precision narrows the reuse region","Total incident photons","Predicted mismatch scale (nm)",(1e5,1e7),(min(ys)*.85,max(ys)*1.15),True,True,"a")
    for d,c in zip((45,50,55),("#3B6FB6","#D8902F","#6B4C9A")):
        g=sorted([x for x in pred if float(x["thickness_nm"])==d],key=lambda x:float(x["total_incident_photons"])); xx=[float(x["total_incident_photons"]) for x in g]; yy=[float(x["delta_star_schur_nm"]) for x in g]; p.line(xx,yy,c,4); p.points(xx,yy,c,radius=7)
    p.legend([("45 nm","#3B6FB6"),("50 nm","#D8902F"),("55 nm","#6B4C9A")])
    x=np.array([float(r["delta_star_schur_nm"]) for r in pred]); y=np.array([float(r["delta_star_general_nm"]) for r in pred]); lo=min(x.min(),y.min())*.9; hi=max(x.max(),y.max())*1.1
    p=Panel(im,rects[1],"Equivalent local-scale formulations","1 / sqrt(Sdelta) (nm)","General local scale (nm)",(lo,hi),(lo,hi),False,False,"b"); p.line([lo,hi],[lo,hi],"#888888",3); p.points(x,y,COLORS["F"],radius=8)
    p=Panel(im,rects[2],"Paired inversion tests local risk","Local predicted D / (v1-v0)","Observed paired D / (v1-v0)",(-1.2,4.0),(-2.2,5.2),False,False,"c"); p.line([-1.2,4],[-1.2,4],"#777777",3); p.hline(0); p.vline(0)
    marker={45:"o",50:"s",55:"^"}; color={1e5:"#9EC1E6",1e6:"#3B6FB6",1e7:"#163A63"}
    for r in paired:
        pr=pmap[r["anchor_id"]]; scale=float(pr["v1_minus_v0_riu2"]); xp=float(r["local_predicted_D_riu2"])/scale; yp=float(r["mean_D_riu2"])/scale; err=1.96*float(r["D_mcse_riu2"])/scale; d=int(float(pr["thickness_nm"])); n=float(pr["total_incident_photons"]); p.errorbar(xp,yp,err,color[n],marker[d])
    paths=save(im,out,"CalSPR_V24_Figure2_threshold_prediction")
    for source in (scan/"V24_THRESHOLD_PREDICTIONS.csv",confirm/"V24_CONFIRMATION_SUMMARY.csv"): copy_source(source,src/source.name)
    return paths


def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""): h.update(chunk)
    return h.hexdigest()


def main(config_path):
    cfg=json.loads(Path(config_path).read_text()); scan=ROOT/cfg["scan_dir"]; confirm=ROOT/cfg["confirmation_dir"]; out=ROOT/cfg["output_dir"]; src=ROOT/cfg["source_data_dir"]
    if out.exists(): raise FileExistsError(f"refusing to overwrite {out}")
    src.mkdir(parents=True); f2=figure2(cfg,scan,confirm,out,src)
    files=[*f2,*sorted(src.iterdir())]; manifest={"status":"V24_FIGURE7_TEXT_CLEANUP_COMPLETE","run_id":cfg["run_id"],"renderer":cfg["renderer"],"files":[{"path":str(p),"size":p.stat().st_size,"sha256":sha(p)} for p in files]}
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n"); print(json.dumps(manifest,indent=2))


if __name__=="__main__":
    if len(sys.argv)!=2: raise SystemExit("usage: make_v24_figures_pillow.py CONFIG")
    main(sys.argv[1])
