#!/usr/bin/env python3
"""Read-only geometry/data audit; write a separate vector QA report."""
import argparse, hashlib, itertools, json, subprocess
from pathlib import Path
import cairo

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('output_dir',type=Path);a=ap.parse_args();out=a.output_dir.resolve()
    root=Path(__file__).resolve().parents[2]
    if out==root or root in out.parents:raise ValueError('Output must be outside the frozen repository tree')
    target=out/'vector_geometry_and_data_audit.json'
    for p in [target,out/'Figure1_PDFproof.png',out/'Figure3_PDFproof.png']:
        if p.exists():raise FileExistsError(p)
    report={'figures':{}}
    for name,size in [('Figure1',(2160,1560)),('Figure3',(2160,900))]:
        cmds=json.loads((out/f'{name}_drawing_commands.json').read_text());texts=[]
        ctx=cairo.Context(cairo.ImageSurface(cairo.FORMAT_ARGB32,1,1))
        for i,c in enumerate(cmds):
            if c['op'] not in ('text','rotated_text'):continue
            ctx.select_font_face('DejaVu Sans',cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_BOLD if c.get('bold') else cairo.FONT_WEIGHT_NORMAL);ctx.set_font_size(c['size'])
            xb,yb,w,h,_,_=ctx.text_extents(c['text'])
            box=[c['x']+xb,c['y']+yb,c['x']+xb+w,c['y']+yb+h] if c['op']=='text' else [c['xy'][0]-h/2,c['xy'][1]-w/2,c['xy'][0]+h/2,c['xy'][1]+w/2]
            covered=False
            for later in cmds[i+1:]:
                if later['op']=='rectangle' and later.get('fill') in ('white','#FFFFFF','#ffffff'):
                    p,q,r,s=later['box']
                    if p<=box[0] and q<=box[1] and r>=box[2] and s>=box[3]:covered=True;break
            if not covered:texts.append({'text':c['text'],'bbox':box})
        outside=[t for t in texts if t['bbox'][0]<0 or t['bbox'][1]<0 or t['bbox'][2]>size[0] or t['bbox'][3]>size[1]]
        overlap=[]
        for p,q in itertools.combinations(texts,2):
            x,y=p['bbox'],q['bbox']
            if min(x[2],y[2])>max(x[0],y[0]) and min(x[3],y[3])>max(x[1],y[1]):overlap.append([p,q])
        subprocess.run(['pdftoppm','-png','-singlefile','-r','300',str(out/f'{name}.pdf'),str(out/f'{name}_PDFproof')],check=True)
        report['figures'][name]={'visible_text_count':len(texts),'outside':outside,'overlapping_text':overlap,'geometric_pass':not outside and not overlap,'pdf_proof':str(out/f'{name}_PDFproof.png')}
    old=Path(__file__).resolve().parents[2]/'figures/revision_20260928/verified_render/plot_data_audit.json';new=out/'plot_data_audit.json'
    report['plot_data_audit_identical_to_approved_repair']=json.loads(old.read_text())==json.loads(new.read_text())
    report['scientific_audit_sha256']={'approved':sha(old),'vector':sha(new)}
    report['all_pass']=all(f['geometric_pass'] for f in report['figures'].values()) and report['plot_data_audit_identical_to_approved_repair']
    target.write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
    if not report['all_pass']:raise SystemExit(2)
if __name__=='__main__':main()
