#!/usr/bin/env python3
"""Replay original and approved display-repair commands into true vector files.

No fit, Monte Carlo, Fisher-information or optimizer execution occurs. The two
display-only reflectance curves use the original renderer and unchanged
parameters; their native Pillow pixels are checked against the frozen plot.
All other numeric marks use the saved CSVs. No raster enters SVG or PDF.
"""
from __future__ import annotations
import argparse, copy, csv, hashlib, importlib.util, json, math, os, platform, subprocess, sys
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import cairo
import PIL
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[2]
ORIGINAL = ROOT/'scripts/release/make_v24_figure7_final_author_cleanup.py'
FIX = ROOT/'scripts/figures/fix_frozen_figures.py'
INPUTS = ROOT/'figures/revision_20260928/inputs'
SCAN = ROOT/'figures/v24_final/source_data'
FROZEN = ROOT/'figures/revision_20260928/verified_render'
RAW_DRAW = ImageDraw.Draw

def now(): return datetime.now(timezone.utc).isoformat()
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod); return mod
def points(xy):
    xy=list(xy)
    return [list(map(float,p)) for p in xy] if isinstance(xy[0],(tuple,list)) else [list(map(float,xy[i:i+2])) for i in range(0,len(xy),2)]

class RecordDraw:
    def __init__(self,image): self.im=image; self.native=RAW_DRAW(image)
    @property
    def commands(self): return self.im.info.setdefault('vector_commands',[])
    def textbbox(self,*args,**kwargs): return self.native.textbbox(*args,**kwargs)
    def text(self,xy,value,font,fill='#222222',anchor=None,**kwargs):
        self.native.text(xy,value,font=font,fill=fill,anchor=anchor,**kwargs)
        actual=self.native.textbbox(xy,value,font=font,anchor=anchor)
        base=font.getbbox(value,anchor='ls')
        self.commands.append({'op':'text','text':value,'x':actual[0]-base[0],'y':actual[1]-base[1],
                              'size':font.size,'bold':'Bold' in str(font.path),'fill':fill,'bbox':actual})
    def line(self,xy,fill,width=1,**kwargs):
        self.native.line(xy,fill=fill,width=width,**kwargs)
        self.commands.append({'op':'line','points':points(xy),'fill':fill,'width':width})
    def rectangle(self,xy,fill=None,outline=None,width=1):
        self.native.rectangle(xy,fill=fill,outline=outline,width=width)
        self.commands.append({'op':'rectangle','box':list(xy),'fill':fill,'outline':outline,'width':width})
    def ellipse(self,xy,fill=None,outline=None,width=1):
        self.native.ellipse(xy,fill=fill,outline=outline,width=width)
        self.commands.append({'op':'ellipse','box':list(xy),'fill':fill,'outline':outline,'width':width})
    def polygon(self,xy,fill=None,outline=None):
        self.native.polygon(xy,fill=fill,outline=outline)
        self.commands.append({'op':'polygon','points':points(xy),'fill':fill,'outline':outline,'width':1})

def instrument(mod):
    mod.ImageDraw=SimpleNamespace(Draw=RecordDraw)
    old_rotate=mod.rotated_text
    if mod.__name__=='original':
        def rotate(image,xy,value,size=28):
            old_rotate(image,xy,value,size)
            image.info.setdefault('vector_commands',[]).append({'op':'rotated_text','xy':xy,'text':value,'size':size,'fill':'#222222'})
    else:
        def rotate(image,xy,value,panel=''):
            old_rotate(image,xy,value,panel)
            image.info.setdefault('vector_commands',[]).append({'op':'rotated_text','xy':xy,'text':value,'size':mod.STYLE['base'],'fill':'#222222'})
    mod.rotated_text=rotate

def color(ctx,c):
    if c=='white': c='#FFFFFF'
    if c=='black': c='#000000'
    ctx.set_source_rgb(*[int(c[k:k+2],16)/255 for k in (1,3,5)])
def paint_path(ctx,c):
    if c.get('fill'):
        color(ctx,c['fill']); ctx.fill_preserve() if c.get('outline') else ctx.fill()
    if c.get('outline'):
        color(ctx,c['outline']); ctx.set_line_width(c['width']); ctx.stroke()
    ctx.new_path()
def replay(ctx,commands):
    color(ctx,'white');ctx.paint()
    for c in commands:
        op=c['op'];ctx.new_path()
        if op=='line':
            for i,(x,y) in enumerate(c['points']): ctx.move_to(x,y) if i==0 else ctx.line_to(x,y)
            color(ctx,c['fill']);ctx.set_line_width(c['width']);ctx.set_line_join(cairo.LINE_JOIN_ROUND);ctx.stroke()
        elif op=='rectangle':
            x0,y0,x1,y1=c['box'];ctx.rectangle(x0,y0,x1-x0,y1-y0);paint_path(ctx,c)
        elif op=='ellipse':
            x0,y0,x1,y1=c['box'];ctx.save();ctx.translate((x0+x1)/2,(y0+y1)/2);ctx.scale((x1-x0)/2,(y1-y0)/2);ctx.arc(0,0,1,0,2*math.pi);ctx.restore();paint_path(ctx,c)
        elif op=='polygon':
            for i,(x,y) in enumerate(c['points']):ctx.move_to(x,y) if i==0 else ctx.line_to(x,y)
            ctx.close_path();paint_path(ctx,c)
        elif op in ('text','rotated_text'):
            ctx.save();ctx.select_font_face('DejaVu Sans',cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_BOLD if c.get('bold') else cairo.FONT_WEIGHT_NORMAL);ctx.set_font_size(c['size'])
            if op=='rotated_text':
                ctx.translate(*c['xy']);ctx.rotate(-math.pi/2);xb,yb,w,h,_,_=ctx.text_extents(c['text']);ctx.move_to(-w/2-xb,-h/2-yb)
            else:ctx.move_to(c['x'],c['y'])
            ctx.text_path(c['text']);color(ctx,c['fill']);ctx.fill();ctx.restore()
        else:raise ValueError(op)

def main():
    sys.dont_write_bytecode=True
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output-dir',type=Path,required=True);args=ap.parse_args();out=args.output_dir.expanduser().resolve()
    if out==ROOT or ROOT in out.parents:
        raise ValueError('Output must be outside the frozen repository tree')
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise FileExistsError(f'Output must be absent or empty: {out}')
    out.mkdir(parents=True,exist_ok=True);start=now()
    input_paths=[Path(__file__).resolve(),ORIGINAL,FIX,ROOT/'scripts/calspr/v24_design.py',ROOT/'scripts/calspr/optics.py',ROOT/'data/raw/materials/Au_Johnson.yml',*sorted(INPUTS.glob('*.csv')),*sorted(SCAN.glob('*.csv')),INPUTS/'Figure6_original.png',INPUTS/'Figure7_original.png',FROZEN/'Figure6_display_fixed.png',FROZEN/'Figure7_display_fixed.png',FROZEN/'plot_data_audit.json']
    for p in input_paths:
        if not p.is_file():raise FileNotFoundError(p)
    assert all(sha(SCAN/p.name)==sha(p) for p in INPUTS.glob('*.csv')), 'Archived source CSV copies differ'
    input_hashes={str(p):sha(p) for p in input_paths}
    original=load('original',ORIGINAL);instrument(original)
    images={}
    def store(image,ignored,stem):images[stem]=image;return []
    original.save=store;original.copy_source=lambda *args:None
    # Fixed deterministic drawing curves only; the original renderer supplies
    # identical 201 samples. All fitted summaries are read from frozen CSVs.
    original.figure1({'material':'data/raw/materials/Au_Johnson.yml'},SCAN,out,out)
    original.figure2({},SCAN,SCAN,out,out)
    im6=images['CalSPR_V24_Figure1_information_sources'];im7=images['CalSPR_V24_Figure2_threshold_prediction']
    comparison={}
    for name,im,box in [('Figure1',im6,(0,80,2160,780)),('Figure3',im7,(0,0,1440,900))]:
        ref=Image.open(INPUTS/('Figure6_original.png' if name=='Figure1' else 'Figure7_original.png')).convert('RGBA')
        comparison[name]={'original_pillow_region':box,'pixel_exact_before_vector_replay':im.crop(box).tobytes()==ref.crop(box).tobytes()}
        assert comparison[name]['pixel_exact_before_vector_replay'], comparison[name]
    # Store display-only curve vertices as an additional reproducibility input.
    curves=[c for c in im6.info['vector_commands'] if c['op']=='line' and len(c['points'])==201]
    assert len(curves)==2
    with (out/'Figure1a_curve_vertices.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['series','sample_index','angle_deg','reflectance'])
        for series,c in zip(['reference_q1.333','sample_q1.338'],curves):
            for i,(x,y) in enumerate(c['points']):w.writerow([series,i,60+(x-150)/895*20,(660-y)/568*1.05])
    fixed=load('fixed',FIX);instrument(fixed)
    def open_base(path):
        value=im6 if Path(path).name=='Figure6_original.png' else im7
        ret=value.copy();ret.info['vector_commands']=copy.deepcopy(value.info['vector_commands']);return ret
    fixed.Image=SimpleNamespace(open=open_base,new=Image.new)
    (out/'panel_crops').mkdir()
    metadata={}
    for old,name,fn in [('Figure6','Figure1',fixed.figure6),('Figure7','Figure3',fixed.figure7)]:
        im,boxes=fn(INPUTS);commands=im.info['vector_commands'];w,h=im.size
        # Assert complete native reconstruction matches approved repaired image.
        ref=Image.open(FROZEN/f'{old}_display_fixed.png').convert('RGBA')
        comparison[name]['complete_native_reconstruction_pixel_exact']=im.tobytes()==ref.tobytes()
        assert comparison[name]['complete_native_reconstruction_pixel_exact']
        for ext,surface_cls in [('svg',cairo.SVGSurface),('pdf',cairo.PDFSurface)]:
            surf=surface_cls(str(out/f'{name}.{ext}'),w*72/300,h*72/300);ctx=cairo.Context(surf);ctx.scale(72/300,72/300);replay(ctx,commands);surf.finish()
        surf=cairo.ImageSurface(cairo.FORMAT_ARGB32,w,h);replay(cairo.Context(surf),commands);surf.write_to_png(str(out/f'{name}.png'))
        png=Image.open(out/f'{name}.png');png.save(out/f'{name}_300dpi.png',dpi=(300,300))
        for letter,box in boxes.items():
            png.crop(box).save(out/'panel_crops'/f'{name}_{letter}.png')
            ref.crop(box).save(out/'panel_crops'/f'{name}_{letter}_approved_raster.png')
        (out/f'{name}_drawing_commands.json').write_text(json.dumps(commands,indent=2),encoding='utf-8')
        svg=(out/f'{name}.svg').read_text();pdfinfo=subprocess.run(['pdfimages','-list',str(out/f'{name}.pdf')],capture_output=True,text=True,check=True).stdout
        (out/f'{name}_pdfimages.txt').write_text(pdfinfo)
        assert '<image' not in svg
        assert len([line for line in pdfinfo.splitlines()[2:] if line.strip()])==0
        metadata[name]={'dimensions_px':im.size,'physical_width_inches':w/300,'panel_count':len(boxes),'vector_commands':len(commands),'svg_raster_images':0,'pdf_raster_images':0,'text_as_vector_outlines':True}
    (out/'plot_data_audit.json').write_text(json.dumps(fixed.PANEL_DATA,indent=2))
    (out/'native_reconstruction_audit.json').write_text(json.dumps(comparison,indent=2))
    assert input_hashes=={str(p):sha(p) for p in input_paths}
    manifest={'start_utc':start,'end_utc':now(),'actual_command':' '.join(sys.argv),'host':platform.node(),'python':sys.version,'pillow':PIL.__version__,'pycairo':cairo.version,'cairo':cairo.cairo_version_string(),'random_seed':'NA: deterministic display regeneration','cpu_threads':1,'scientific_fit_recomputed':False,'monte_carlo_rerun':False,'information_metrics_recomputed':False,'display_reflectance_curves_re_evaluated_with_original_renderer':True,'native_reconstruction':comparison,'input_hashes':input_hashes,'input_hashes_unchanged':True,'figures':metadata,'outputs':{str(p):sha(p) for p in sorted(out.rglob('*')) if p.is_file()},'exit_status':0}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2));print(json.dumps({k:v for k,v in manifest.items() if k not in ('input_hashes','outputs')},indent=2))

if __name__=='__main__':main()
