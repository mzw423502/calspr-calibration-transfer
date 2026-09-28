#!/usr/bin/env python3
"""Create the two bounded publication figures for CalSPR V2.4."""
from __future__ import annotations
import csv, json, math, shutil, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from calspr.figure_style_kernel import apply_figure_style, panel_letter
from calspr.v24_design import ResponseScenario, angle_centres, gold_reflectance, physical_windows

COLORS={"F":"#3B6FB6","C":"#D8902F","S":"#6B4C9A","P":"#4C956C","O":"#8A8A8A"}

def rows(path):
    with open(path,newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def save(fig,out,stem):
    png=out/f"{stem}.png"; pdf=out/f"{stem}.pdf"; fig.savefig(png); fig.savefig(pdf); plt.close(fig); return png,pdf
def copy_source(path,target):
    if target.exists(): raise FileExistsError(target)
    shutil.copy2(path,target)

def figure1(cfg,scan,out,src):
    attrib=rows(scan/"V24_INFORMATION_ATTRIBUTION.csv"); equal=rows(scan/"V24_EQUAL_BUDGET_COMPARISON.csv"); checks=rows(scan/"V24_DESIGN_NUMERICAL_CHECKS.csv")
    base=json.loads((scan/"summary.json").read_text()); model_cfg={"angle_min_deg":60.0,"angle_max_deg":80.0,"angle_bins":201,"baseline_index":1.333,"sample_delta_index":0.005,"wavelength_nm":633.0,"window_definition":{"nominal_thickness_nm":50.0,"critical_half_width_deg":1.0,"spr_half_width_deg":2.0}}
    angles=angle_centres(model_cfg); scenario=ResponseScenario("top_hat_0p10deg_main","top_hat",0.1,16); gold=ROOT/cfg["material"]; win=physical_windows(model_cfg,angles,gold)
    rr=gold_reflectance(angles,633.0,50.0,1.333,gold,scenario); rs=gold_reflectance(angles,633.0,50.0,1.338,gold,scenario)
    apply_figure_style(sizes=(8,7,6)); fig,axs=plt.subplots(2,2,figsize=(7.2,5.2)); ax=axs[0,0]
    ax.axvspan(60.7,62.6,color=COLORS["C"],alpha=.16,label="Critical window (C)"); ax.axvspan(71,75,color=COLORS["S"],alpha=.14,label="SPR window (S)")
    ax.plot(angles,rr,color="#777777",label="Reference, q=1.333"); ax.plot(angles,rs,color=COLORS["F"],label="Sample, q=1.338"); ax.set(xlabel="Incident angle (degrees)",ylabel="Reflectance",title="Frozen windows separate the critical and SPR features"); ax.legend(loc="lower left"); panel_letter(ax,"a")
    ax=axs[0,1]; contrib=[]
    for name in ("C","S","O"):
        r=next(x for x in attrib if x["design"]==f"F_projected_score_contribution_{name}"); contrib.append(float(r["effective_information"])/float(next(x for x in attrib if x["design"]=="F")["effective_information"])*100)
    ax.barh(["Critical C","SPR S","Other O"],contrib,color=[COLORS["C"],COLORS["S"],COLORS["O"]]); ax.set(xlabel="Share of jointly projected F information (%)",title="The SPR window carries most full-design target information"); ax.set_xlim(0,max(contrib)*1.15); panel_letter(ax,"b")
    ax=axs[1,0]; selected=[x for x in equal if x["selected_within_strategy"]=="1" and x["strategy"] in ("F","C","S","P")]; order=["F","C","S","P"]; y=[float(next(x for x in selected if x["strategy"]==s)["local_sd_riu"])*1e6 for s in order]
    ax.scatter(range(4),y,s=34,color=[COLORS[s] for s in order],zorder=3); [ax.plot([i,i],[0,v],color="#CCCCCC",lw=1,zorder=1) for i,v in enumerate(y)]; ax.set_xticks(range(4),["Full F","Critical C","SPR S","Bare prism P"]); ax.set(ylabel="Local SD (micro-RIU)",title="Equal budget favors the focused SPR window"); ax.set_ylim(0,max(y)*1.12); panel_letter(ax,"c")
    ax=axs[1,1]
    phase=[x for x in checks if x["check_type"]=="grid_phase" and x["strategy"] in order]
    for s in order:
        g=sorted([x for x in phase if x["strategy"]==s],key=lambda x:float(x["phase_fraction_of_spacing"])); yy=np.array([float(x["local_sd_riu"]) for x in g]); yy=yy/yy[1]; xx=[float(x["phase_fraction_of_spacing"]) for x in g]; ax.plot(xx,yy,marker="o",color=COLORS[s],label=s)
    ax.axhline(1,color="#BBBBBB",lw=.8); ax.set(xlabel="Grid phase (fraction of 0.1-degree spacing)",ylabel="Local SD / centre-phase SD",title="Critical-angle designs retain grid-phase sensitivity"); ax.legend(ncol=2); panel_letter(ax,"d")
    fig.tight_layout(); paths=save(fig,out,"CalSPR_V24_Figure1_information_sources")
    for p in (scan/"V24_INFORMATION_ATTRIBUTION.csv",scan/"V24_EQUAL_BUDGET_COMPARISON.csv",scan/"V24_DESIGN_NUMERICAL_CHECKS.csv"):copy_source(p,src/p.name)
    return paths

def figure2(cfg,scan,confirm,out,src):
    pred=rows(scan/"V24_THRESHOLD_PREDICTIONS.csv"); summ=rows(confirm/"V24_CONFIRMATION_SUMMARY.csv"); paired=[x for x in summ if x["summary_type"]=="paired_D_T0_minus_T1"]
    pmap={x["anchor_id"]:x for x in pred}; apply_figure_style(sizes=(8,7,6)); fig=plt.figure(figsize=(7.2,3.0)); gs=fig.add_gridspec(1,3,wspace=.42); ax=fig.add_subplot(gs[0,0])
    for d,c in zip((45,50,55),("#3B6FB6","#D8902F","#6B4C9A")):
        g=sorted([x for x in pred if float(x["thickness_nm"])==d],key=lambda x:float(x["total_incident_photons"])); ax.plot([float(x["total_incident_photons"]) for x in g],[float(x["delta_star_schur_nm"]) for x in g],marker="o",color=c,label=f"{d} nm")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set(xlabel="Total incident photons",ylabel="Predicted threshold (nm)",title="Higher precision narrows the reuse region"); ax.legend(); panel_letter(ax,"a",dx=-.28)
    ax=fig.add_subplot(gs[0,1]); x=np.array([float(r["delta_star_schur_nm"]) for r in pred]); y=np.array([float(r["delta_star_general_nm"]) for r in pred]); lo=min(x.min(),y.min())*.9; hi=max(x.max(),y.max())*1.1; ax.plot([lo,hi],[lo,hi],color="#888888",lw=1); ax.scatter(x,y,color=COLORS["F"],s=26); ax.set(xlabel="1 / sqrt(Sdelta) (nm)",ylabel="sqrt(v1-v0) / |Kdelta| (nm)",title="Both threshold forms coincide",xlim=(lo,hi),ylim=(lo,hi)); panel_letter(ax,"b",dx=-.28)
    ax=fig.add_subplot(gs[0,2]); marker={45:"o",50:"s",55:"^"}; color={1e5:"#9EC1E6",1e6:"#3B6FB6",1e7:"#163A63"}
    for r in paired:
        p=pmap[r["anchor_id"]]; scale=float(p["v1_minus_v0_riu2"]); xp=float(r["local_predicted_D_riu2"])/scale; yp=float(r["mean_D_riu2"])/scale; err=1.96*float(r["D_mcse_riu2"])/scale; d=int(float(p["thickness_nm"])); n=float(p["total_incident_photons"]); ax.errorbar(xp,yp,yerr=err,fmt=marker[d],ms=3,color=color[n],alpha=.48,lw=.55,capsize=0)
    lim=4.0; ax.plot([-1.2,lim],[-1.2,lim],color="#777777",lw=1); ax.axhline(0,color="#BBBBBB",lw=.7); ax.axvline(0,color="#BBBBBB",lw=.7); ax.set(xlim=(-1.2,lim),ylim=(-2.2,5.2),xlabel="Local predicted D / (v1-v0)",ylabel="Observed paired D / (v1-v0)",title="Paired inversion tests the local risk prediction"); panel_letter(ax,"c",dx=-.28)
    fig.tight_layout(); paths=save(fig,out,"CalSPR_V24_Figure2_threshold_prediction")
    for p in (scan/"V24_THRESHOLD_PREDICTIONS.csv",confirm/"V24_CONFIRMATION_SUMMARY.csv"):copy_source(p,src/p.name)
    return paths

def main(cp):
    cfg=json.loads(Path(cp).read_text()); scan=ROOT/cfg["scan_dir"]; confirm=ROOT/cfg["confirmation_dir"]; out=ROOT/cfg["output_dir"]; src=ROOT/cfg["source_data_dir"]
    if out.exists() and any(out.iterdir()):raise FileExistsError(f"refusing to overwrite nonempty {out}")
    out.mkdir(parents=True,exist_ok=True); src.mkdir(parents=True,exist_ok=True); f1=figure1(cfg,scan,out,src); f2=figure2(cfg,scan,confirm,out,src)
    manifest={"status":"V24_FIGURES_COMPLETE","figure_files":[str(x) for x in (*f1,*f2)],"source_data_files":[str(x) for x in sorted(src.iterdir())]}; (out/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n"); print(json.dumps(manifest,indent=2))
if __name__=="__main__":
    if len(sys.argv)!=2:raise SystemExit("usage: make_v24_figures.py CONFIG")
    main(sys.argv[1])
