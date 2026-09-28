#!/usr/bin/env python3
"""Render CalSPR V2.2 manuscript figures from verified remote-result snapshots.

This script is intentionally run on Windows because the remote Matplotlib
runtime is execution-blocked. It performs no stochastic or inverse calculation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_ROOT/"scripts"))
from vendor.figure_style import apply_figure_style, panel_letter

COLORS = {"T0":"#666666", "T1":"#0072B2", "T2":"#E69F00", "T_R":"#7A5195",
          "O":"#222222", "fixed":"#0072B2", "aware":"#CC79A7"}


def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest()


def rows(path):
    with open(path,newline="",encoding="utf-8") as f: return list(csv.DictReader(f))


def num(row,key): return float(row[key])


def save_all(fig, stem):
    fig.savefig(stem.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    plt.close(fig)


def figure1(data, out, boundary):
    paired=rows(data["symmetric"]); cond=rows(data["boundary_conditions"])
    fig=plt.figure(figsize=(7.1,6.0),layout="constrained")
    gs=fig.add_gridspec(2,2,height_ratios=[0.72,1.0])
    ax=fig.add_subplot(gs[0,:]); ax.axis("off"); panel_letter(ax,"a")
    ax.set_title("Four actionable choices spend the same photon budget",loc="left")
    cards=[("T0\nrigid reuse","source ref 50%","target sample 50%","thickness difference fixed"),
           ("T1\nestimate difference","source ref 50%","target sample 50%","thickness difference fitted"),
           ("T2\ntransfer + target ref","source ref + target ref","target sample","three-way allocation"),
           ("T_R\nrecalibrate target","target ref","target sample","no source parameter")]
    for j,(title,a,b,c) in enumerate(cards):
        x=0.02+j*0.245
        ax.add_patch(plt.Rectangle((x,0.16),0.22,0.65,facecolor="#F7F7F7",edgecolor=COLORS[title.split()[0]],lw=1.4))
        ax.text(x+0.11,0.70,title,ha="center",va="center",weight="bold",fontsize=7,color=COLORS[title.split()[0]],linespacing=1.25)
        ax.text(x+0.11,0.51,a,ha="center",va="center")
        ax.text(x+0.11,0.36,b,ha="center",va="center")
        ax.text(x+0.11,0.21,c,ha="center",va="center",fontsize=6,color="#555555")
    ax.text(0.5,0.02,"Total incident photons = 10M in every strategy",ha="center",va="bottom",fontsize=7)

    ax=fig.add_subplot(gs[1,0]); panel_letter(ax,"b")
    x=np.array([num(r,"absolute_mismatch_nm") for r in paired])
    y=np.array([num(r,"mean_mse_difference_T0_minus_T1_riu2") for r in paired])*1e11
    lo=np.array([num(r,"mc_interval_low_riu2") for r in paired])*1e11
    hi=np.array([num(r,"mc_interval_high_riu2") for r in paired])*1e11
    ax.errorbar(x,y,yerr=[y-lo,hi-y],fmt="o-",color=COLORS["T1"],capsize=3,lw=1.5)
    ax.axhline(0,color="#888888",lw=0.8)
    ax.axvline(boundary,color="#333333",ls="--",lw=1)
    ax.axvspan(0.03,0.05,color="#E5E5E5",alpha=.6,zorder=0)
    ax.text(boundary+0.001,13.5,"local prediction\n0.031 nm",fontsize=7,va="top")
    ax.set(xlabel="Absolute thickness mismatch (nm)",ylabel="Mean paired MSE difference\nT0 - T1 (10$^{-11}$ RIU$^2$)")
    ax.set_title("The observed transition is a band, not a point",loc="left")
    ax.margins(x=.06,y=.12)

    ax=fig.add_subplot(gs[1,1]); panel_letter(ax,"c")
    methods=["T0","T1"]
    for method in methods:
        rr=[r for r in cond if r["scenario"]=="thickness_mismatch" and r["strategy"]==method]
        xs=np.array([abs(num(r,"true_mismatch_nm")) for r in rr])
        ys=np.array([num(r,"rmse_riu") for r in rr])*1e5
        ax.scatter(xs,ys,s=26,color=COLORS[method],alpha=.55,label=method)
        ux=sorted(set(xs))
        mean=[np.mean(ys[xs==z]) for z in ux]
        ax.plot(ux,mean,"-",color=COLORS[method],lw=1.5)
    ax.axvspan(0.03,0.05,color="#E5E5E5",alpha=.6,zorder=0)
    ax.set(xlabel="Absolute thickness mismatch (nm)",ylabel="RMSE (10$^{-5}$ RIU)")
    ax.set_title("Rigid-transfer error rises beyond the boundary",loc="left")
    ax.legend(frameon=False); ax.margins(.08)
    save_all(fig,out/"Figure1_transfer_boundary")


def figure2(data,out,boundary):
    local=rows(data["local"]); pseudo=rows(data["pseudo"])
    anchor=next(r for r in local if num(r,"budget_photons")==1e7 and num(r,"source_reference_fraction")==.5 and num(r,"sample_delta_index")==.005 and num(r,"nominal_thickness_nm")==50)
    fig,axs=plt.subplots(2,2,figsize=(7.1,5.8),layout="constrained")
    ax=axs[0,0]; panel_letter(ax,"a"); ax.set_title("Bias and variance trade places at the crossing",loc="left")
    ax.axis("off")
    ax.annotate("",(.75,.72),(.18,.28),arrowprops=dict(arrowstyle="->",lw=2,color=COLORS["T1"]))
    ax.annotate("",(.63,.43),(.18,.28),arrowprops=dict(arrowstyle="->",lw=2,color="#D55E00"))
    ax.text(.76,.73,"target derivative  $r_q$",color=COLORS["T1"],ha="right",va="bottom",fontsize=7)
    ax.text(.64,.44,"omitted mismatch  $r_\\delta$",color="#D55E00",ha="right",va="bottom",fontsize=7)
    ax.text(.06,.03,"$K_\\delta=(r_q^T r_\\delta)/(r_q^T r_q)$\n$MSE_{T0}\\approx v_0+(K_\\delta\\delta)^2$\n$MSE_{T1}\\approx v_1$",fontsize=8)

    ax=axs[0,1]; panel_letter(ax,"b")
    k=num(anchor,"K_delta_riu_per_nm"); s0=num(anchor,"T0_local_sd_riu"); s1=num(anchor,"T1_local_sd_riu")
    x=np.linspace(0,.10,200)
    ax.plot(x,np.sqrt(s0*s0+(k*x)**2)*1e5,color=COLORS["T0"],lw=2,label="T0 rigid reuse")
    ax.axhline(s1*1e5,color=COLORS["T1"],lw=2,label="T1 fit thickness")
    ax.axvline(boundary,color="#333333",ls="--",lw=1)
    ax.set(xlabel="Absolute thickness mismatch (nm)",ylabel="Local predicted RMSE (10$^{-5}$ RIU)")
    ax.set_title("T0 pays bias; T1 pays nuisance variance",loc="left"); ax.legend(frameon=False); ax.margins(.04)

    ax=axs[1,0]; panel_letter(ax,"c")
    rr=[r for r in pseudo if r["strategy"]=="T0" and abs(num(r,"true_mismatch_nm"))<=.2]
    x=np.array([num(r,"true_mismatch_nm") for r in rr]); obs=np.array([num(r,"q_bias_riu") for r in rr])*1e5
    pred=np.array([num(r,"linear_predicted_T0_bias_riu") for r in rr])*1e5
    order=np.argsort(x)
    ax.plot(x[order],pred[order],"--",color="#333333",lw=1.5,label="forward derivative")
    ax.scatter(x,obs,facecolor="white",edgecolor=COLORS["T0"],s=28,label="pseudo-true fit")
    ax.axhline(0,color="#AAAAAA",lw=.8)
    ax.set(xlabel="Signed thickness mismatch (nm)",ylabel="T0 index bias (10$^{-5}$ RIU)")
    ax.set_title("First-order bias matches finite fits locally",loc="left"); ax.legend(frameon=False); ax.margins(.07)

    ax=axs[1,1]; panel_letter(ax,"d")
    rr=[r for r in local if num(r,"source_reference_fraction")==.5 and num(r,"sample_delta_index")==.005 and num(r,"nominal_thickness_nm")==50]
    rr=sorted(rr,key=lambda r:num(r,"budget_photons")); x=np.array([num(r,"budget_photons") for r in rr]); y=np.array([num(r,"T0_vs_T1_boundary_abs_delta_nm") for r in rr])
    ax.loglog(x,y,"o-",color=COLORS["T1"],lw=1.5)
    guide=y[-1]*(x/x[-1])**-.5; ax.loglog(x,guide,"--",color="#666666",lw=1,label="$N^{-1/2}$")
    ax.set(xlabel="Total incident photons",ylabel="Predicted crossing |delta| (nm)")
    ax.set_title("Higher photon budgets narrow rigid reuse",loc="left"); ax.legend(frameon=False); ax.margins(.08)
    save_all(fig,out/"Figure2_geometry_and_pseudotrue")


def figure3(data,out):
    cond=rows(data["boundary_conditions"])
    fig,axs=plt.subplots(2,2,figsize=(7.1,5.8),layout="constrained")
    thick=[r for r in cond if r["scenario"]=="thickness_mismatch" and r["strategy"] in ("T0","T1","T2")]
    for ax,metric,label,title in ((axs[0,0],"bias_riu","Bias (10$^{-5}$ RIU)","Bias changes sign with mismatch"),
                                  (axs[0,1],"rmse_riu","RMSE (10$^{-5}$ RIU)","Risk separates only outside the crossing band"),
                                  (axs[1,0],"p95_absolute_error_riu","95th-percentile |error| (10$^{-5}$ RIU)","Tail error remains finite across the grid")):
        for method in ("T0","T1","T2"):
            rr=sorted([r for r in thick if r["strategy"]==method],key=lambda r:num(r,"true_mismatch_nm"))
            ax.plot([num(r,"true_mismatch_nm") for r in rr],[num(r,metric)*1e5 for r in rr],"o-",label=method,color=COLORS[method],lw=1.3,ms=4)
        ax.axhline(0,color="#AAAAAA",lw=.7); ax.set(xlabel="Signed thickness mismatch (nm)",ylabel=label); ax.set_title(title,loc="left"); ax.margins(.06)
    axs[0,0].legend(frameon=False,ncol=3)
    ax=axs[1,1]; panel_letter(ax,"d")
    ang=[r for r in cond if r["scenario"]=="target_angle_offset"]
    order=["T1","T1_angle","T2","T2_angle"]
    x=np.arange(len(order)); width=.34
    for j,tau in enumerate((-.002,.002)):
        vals=[num(next(r for r in ang if r["strategy"]==m and num(r,"true_target_angle_offset_deg")==tau),"rmse_riu")*1e5 for m in order]
        ax.bar(x+(j-.5)*width,vals,width,label=f"{tau:+.3f} deg",color="#56B4E9" if j==0 else "#D55E00",alpha=.85)
    ax.set_xticks(x, ["T1 fixed","T1 + angle","T2 fixed","T2 + angle"],rotation=18,ha="right")
    ax.set_ylabel("RMSE (10$^{-5}$ RIU)"); ax.set_title("Target reference prevents angle-fit variance inflation",loc="left"); ax.legend(frameon=False)
    for L,ax in zip("abc",axs.flat[:3]): panel_letter(ax,L)
    save_all(fig,out/"Figure3_frozen_confirmation")


def figure4(data,out):
    scan=rows(data["tr_allocation"]); confirm=rows(data["tr_conditions"]); diffs=rows(data["tr_differences"])
    fig,axs=plt.subplots(2,2,figsize=(7.1,5.8),layout="constrained")
    ax=axs[0,0]; panel_letter(ax,"a")
    for method,marker in (("T1","o"),("T2","s"),("T_R","D")):
        rr=[r for r in scan if r["strategy"]==method and r["angle_model"] in ("stable_angle","target_common_angle_fitted")]
        ax.scatter([num(r,"source_reference_fraction") for r in rr],[num(r,"local_sd_riu")*1e5 for r in rr],s=28,marker=marker,color=COLORS[method],alpha=.75,label=method)
    ax.set(xlabel="Source-reference photon fraction",ylabel="Local SD (10$^{-5}$ RIU)")
    ax.set_title("Stable-angle T2 moves toward the zero-source boundary",loc="left"); ax.legend(frameon=False); ax.margins(.08)

    ax=axs[0,1]; panel_letter(ax,"b")
    for method,marker in (("T1","o"),("T2","s"),("T_R","D")):
        rr=[r for r in scan if r["strategy"]==method and r["angle_model"] in ("angle_aware","target_common_angle_fitted")]
        ax.scatter([num(r,"source_reference_fraction") for r in rr],[num(r,"local_sd_riu")*1e5 for r in rr],s=28,marker=marker,color=COLORS[method],alpha=.75,label=method)
    ax.set(xlabel="Source-reference photon fraction",ylabel="Local SD (10$^{-5}$ RIU)")
    ax.set_title("Target-only calibration retains angle information",loc="left"); ax.set_yscale("log"); ax.margins(.08)

    ax=axs[1,0]; panel_letter(ax,"c")
    scenarios=sorted({r["scenario"] for r in confirm}); x=np.arange(len(scenarios)); width=.24
    for j,method in enumerate(("T1","T2","T_R")):
        vals=[num(next(r for r in confirm if r["scenario"]==s and r["strategy"]==method),"rmse_riu")*1e5 for s in scenarios]
        ax.bar(x+(j-1)*width,vals,width,label=method,color=COLORS[method])
    ax.set_xticks(x,["-0.03 nm","+0.03 nm","-0.002 deg","+0.002 deg"],rotation=18,ha="right")
    ax.set_ylabel("RMSE (10$^{-5}$ RIU)"); ax.set_title("Noisy confirmation compares complete protocols",loc="left"); ax.legend(frameon=False,ncol=3)

    ax=axs[1,1]; panel_letter(ax,"d")
    rr=[r for r in diffs if r["strategy_b"]=="T_R" and r["strategy_a"] in ("T1","T2")]
    labels=[]; y=[]; lo=[]; hi=[]
    for r in rr:
        labels.append(("T1" if r["strategy_a"]=="T1" else "T2")+" | "+r["scenario"].replace("stable_angle_near_boundary_","").replace("target_reinstallation_","angle "))
        v=num(r,"mean_mse_difference_a_minus_b_riu2")*1e10; y.append(v); lo.append(v-num(r,"mc_interval_low_riu2")*1e10); hi.append(num(r,"mc_interval_high_riu2")*1e10-v)
    yy=np.arange(len(y)); ax.errorbar(y,yy,xerr=[lo,hi],fmt="o",color="#333333",capsize=3)
    ax.axvline(0,color="#888888",lw=.8); ax.set_xscale("symlog",linthresh=1.0,linscale=1.0); ax.set_yticks(yy,labels); ax.set_xlabel("MSE difference versus T_R (10$^{-10}$ RIU$^2$; symlog)")
    ax.set_title("Source reference adds no resolved benefit",loc="left"); ax.invert_yaxis(); ax.margins(y=.06)
    save_all(fig,out/"Figure4_target_recalibration")


def figure5(data,out):
    mat=rows(data["material"])
    fig=plt.figure(figsize=(7.1,5.8),layout="constrained"); gs=fig.add_gridspec(2,2,height_ratios=[.72,1])
    ax=fig.add_subplot(gs[0,:]); panel_letter(ax,"a"); ax.axis("off"); ax.set_title("Material conditions and model mismatch answer different questions",loc="left")
    cards=[("Matched condition","same named model\ngenerates and fits","conditional comparison"),
           ("Public-model mismatch","Olmon or Rakic truth\nJC fit","forward-family error"),
           ("Interdevice change","JC source + non-JC target\nJC fit","cross-device material change")]
    for j,(t,b,c) in enumerate(cards):
        x=.04+j*.32; ax.add_patch(plt.Rectangle((x,.18),.27,.60,facecolor="#F7F7F7",edgecolor="#666666",lw=1.2)); ax.text(x+.135,.68,t,ha="center",weight="bold"); ax.text(x+.135,.45,b,ha="center"); ax.text(x+.135,.23,c,ha="center",fontsize=7,color="#555555")
    ax=fig.add_subplot(gs[1,0]); panel_letter(ax,"b")
    cases=[]
    for r in mat:
        label=("same non-JC" if r["case"]=="public_model_mismatch" else "JC to non-JC")+"\n"+r["truth_target_material"].replace("_"," ")
        if label not in cases: cases.append(label)
    x=np.arange(len(cases)); width=.18
    for j,method in enumerate(("T0","T1","T2","O")):
        vals=[]
        for case in cases:
            rr=next(r for r in mat if r["strategy"]==method and (("same non-JC" if r["case"]=="public_model_mismatch" else "JC to non-JC")+"\n"+r["truth_target_material"].replace("_"," "))==case)
            vals.append(num(rr,"q_bias_riu")*1e3)
        ax.bar(x+(j-1.5)*width,vals,width,label=method,color=COLORS[method])
    ax.axhline(0,color="#888888",lw=.8); ax.set_xticks(x,cases,rotation=16,ha="right"); ax.set_ylabel("Index bias (10$^{-3}$ RIU)"); ax.set_title("Optical error persists after fitting thickness",loc="left"); ax.legend(frameon=False,ncol=4)
    ax=fig.add_subplot(gs[1,1]); panel_letter(ax,"c")
    for method in ("T0","T1","T2","O"):
        rr=[r for r in mat if r["strategy"]==method]
        ax.scatter([abs(num(r,"q_bias_riu")) for r in rr],[num(r,"objective_half_deviance") for r in rr],s=34,color=COLORS[method],label=method,alpha=.8)
    ax.set(xlabel="Absolute index bias (RIU)",ylabel="Half deviance")
    ax.set_yscale("log"); ax.set_title("Large bias coincides with lack of fit",loc="left"); ax.legend(frameon=False); ax.margins(.08)
    save_all(fig,out/"Figure5_material_boundary")


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("project_root"); ap.add_argument("output_root")
    a=ap.parse_args(); root=Path(a.project_root).resolve(); out=Path(a.output_root).resolve()
    if out.exists(): raise FileExistsError(f"refusing to overwrite {out}")
    out.mkdir(parents=True); source=out.parent/f"{out.name}_source_data"
    if source.exists(): raise FileExistsError(f"refusing to overwrite {source}")
    source.mkdir(parents=True)
    data={
        "symmetric":root/"outputs/integrity_audit/boundary_confirmation_audit_20260913_v22/symmetric_magnitude_risk_differences.csv",
        "boundary_conditions":root/"outputs/integrity_audit/boundary_confirmation_audit_20260913_v22/condition_summary_verified.csv",
        "local":root/"outputs/boundary_scan/boundary_20260913_e1/local_geometry.csv",
        "pseudo":root/"outputs/boundary_scan/boundary_20260913_e1/pseudo_true_boundary.csv",
        "material":root/"outputs/boundary_scan/boundary_20260913_e1/material_mismatch.csv",
        "tr_allocation":root/"outputs/target_recalibration_scan/target_recalibration_scan_20260913_e1/allocation_scan.csv",
        "tr_conditions":root/"outputs/target_recalibration_confirmation/target_recalibration_confirm_20260913_c1/condition_summary.csv",
        "tr_differences":root/"outputs/target_recalibration_confirmation/target_recalibration_confirm_20260913_c1/independent_risk_differences.csv",
        "boundary_summary":root/"outputs/boundary_scan/boundary_20260913_e1/summary.json",
    }
    for key,path in data.items():
        if not path.is_file(): raise FileNotFoundError(f"missing {key}: {path}")
        shutil.copy2(path,source/f"{key}{path.suffix}")
    boundary=json.loads(data["boundary_summary"].read_text(encoding="utf-8"))["anchor_T0_vs_T1_boundary_abs_delta_nm"]
    apply_figure_style(frame="open",sizes=(8,7,6),grid=False)
    figure1(data,out,boundary); figure2(data,out,boundary); figure3(data,out); figure4(data,out); figure5(data,out)
    manifest=[]
    for path in sorted(list(out.glob("Figure*"))+list(source.glob("*"))):
        manifest.append({"path":str(path),"size_bytes":path.stat().st_size,"sha256":sha(path)})
    (out/"figure_manifest.json").write_text(json.dumps({"status":"RENDERED_PENDING_VISUAL_REVIEW","matplotlib":plt.matplotlib.__version__,"files":manifest},indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":"RENDERED_PENDING_VISUAL_REVIEW","output":str(out),"file_count":len(manifest)},indent=2))


if __name__=="__main__": main()
