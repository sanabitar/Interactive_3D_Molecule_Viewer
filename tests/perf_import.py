#!/usr/bin/env python3
"""Import timing in real Chromium (Three.js SHIM: measures parsing, bond perception, bond orders, selection UI; NOT GPU rendering)."""
import sys,time,tempfile,json;sys.argv=['x'];__file__='/home/claude/laura/tests/browser/run_browser_tests.py';exec(open(__file__).read().split("@case('01")[0])
src=(FX/'alpha12.pdb').read_text().split('\n');atoms=[l for l in src if l.startswith('ATOM')];d=Path(tempfile.mkdtemp())
def make(n_copies,name):
    out=[];serial=1
    for c in range(n_copies):
        dx=(c%50)*30.0;dy=((c//50)%50)*30.0;dz=(c//2500)*30.0
        for l in atoms:
            x=float(l[30:38])+dx;y=float(l[38:46])+dy;z=float(l[46:54])+dz
            out.append(f"ATOM  {serial%100000:5d}"+l[11:30]+f"{x:8.3f}{y:8.3f}{z:8.3f}"+l[54:])
            serial+=1
    p=d/name;p.write_text('\n'.join(out)+'\nEND\n');return p
rows=[]
with sync_playwright() as pw:
    for copies,label in[(int(c),c+' copies') for c in __import__('os').environ.get('PERF','1,100,500,2000,4000').split(',')]:
        f=make(copies,f'big{copies}.pdb');mb=f.stat().st_size/1e6;pg=new_page(pw);t=time.time()
        pg.set_input_files('#file',str(f))
        try:pg.wait_for_function("(n)=>document.querySelector('#cmn').textContent===n",arg=f.name,timeout=170000);ms=(time.time()-t)*1000;st='loaded'
        except Exception as e:ms=(time.time()-t)*1000;st='TIMEOUT/not loaded: '+pg.inner_text('#ov')[:80]
        info=pg.inner_text('#cmm');rows.append((label,round(mb,1),st,int(ms),info,pg.errs[:1]))
        if copies==500:
            sec(pg,'Teaching, session');pg.click('#bl');pg.wait_for_timeout(300);rows.append(('  Blender export on 30k atoms',None,'toast: '+pg.inner_text('#toast')[:95],0,'',[]))
        pg._b.close()
    big=d/'over_limit.pdb'
    with open(big,'w') as fh:
        line=atoms[0]+'\n'
        for _ in range(int(61e6/len(line))+1000):fh.write(line)
    pg=new_page(pw);t=time.time();pg.set_input_files('#file',str(big));pg.wait_for_selector('#ov >> text=failed',timeout=120000);rows.append(('61 MB file (over the 60 MB limit)',round(big.stat().st_size/1e6,1),'rejected: '+pg.inner_text('#ov')[:90].replace('\n',' '),int((time.time()-t)*1000),'',pg.errs[:1]));pg._b.close()
for r in rows:print(f'{r[0]:34} {str(r[1]) if r[1] else "":>6} MB  {r[2][:78]:78} {r[3]:>7} ms  {r[4]}  {r[5] or ""}')
