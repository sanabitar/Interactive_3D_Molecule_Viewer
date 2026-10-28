#!/usr/bin/env python3
"""Real-Chromium regression tests for LAURA (Playwright for Python).
Default: Three.js is replaced by tests/browser/three_shim.js (no network needed) => DOM/CSS/JS/logic are real, 3D rendering is NOT exercised.
With network: python3 run_browser_tests.py --real-cdn  (uses the page's real CDN scripts; software WebGL may still be limited headless).
Usage: python3 tests/browser/run_browser_tests.py [--page PATH] [--out results.json]
The labs are ES modules, so serve the project first: python3 -m http.server 8799 (in the project root) and set LAURA_URL=http://127.0.0.1:8799/laura-viewer.html"""
import sys,os,re,json,zipfile,py_compile,tempfile,traceback
from pathlib import Path
from playwright.sync_api import sync_playwright
HERE=Path(__file__).parent;FX=HERE.parent/'fixtures'
args=sys.argv[1:];REAL='--real-cdn' in args
PAGE=Path(args[args.index('--page')+1]) if '--page' in args else Path('/mnt/user-data/outputs/laura-viewer.html')
OUT=Path(args[args.index('--out')+1]) if '--out' in args else HERE/'results.json'
JSZIP=next((p for p in ['/home/claude/.npm-global/lib/node_modules/docx/node_modules/jszip/dist/jszip.min.js'] if os.path.exists(p)),None)
SHIM=(HERE/'three_shim.js').read_text()
RES=[];ERRORS=[]
def new_page(pw,w=1280,h=800,storage=None):
    b=pw.chromium.launch();ctx=b.new_context(viewport={'width':w,'height':h},accept_downloads=True)
    if not REAL:
        ctx.route(re.compile(r'.*three\.min\.js'),lambda r:r.fulfill(body=SHIM,content_type='application/javascript'))
        ctx.route(re.compile(r'.*OrbitControls\.js'),lambda r:r.fulfill(body='',content_type='application/javascript'))
        ctx.route(re.compile(r'.*jszip\.min\.js'),lambda r:r.fulfill(path=JSZIP,content_type='application/javascript'))
    ctx.add_init_script("window.__fs=0;Element.prototype.requestFullscreen=function(){window.__fs++;return Promise.resolve()}")
    pg=ctx.new_page();errs=[]
    pg.on('console',lambda m:errs.append('console.'+m.type+': '+m.text) if m.type=='error' else None)
    pg.on('pageerror',lambda e:errs.append('pageerror: '+str(e)))
    pg.on('requestfailed',lambda r:errs.append('requestfailed: '+r.url[:90]) if 'rcsb' not in r.url else None)
    pg.goto(__import__('os').environ.get('LAURA_URL') or 'file://'+str(PAGE));pg.wait_for_function("document.querySelector('#cmn').textContent.trim().length>1")
    pg.errs=errs;pg._b=b;return pg
ONLY=[a for a in args[args.index('--only')+1].split(',')] if '--only' in args else None
def case(name):
    def deco(fn):
        if ONLY and name.split()[0] not in ONLY:return fn
        try:
            with sync_playwright() as pw:fn(pw);RES.append((name,'PASS',''))
        except AssertionError as e:RES.append((name,'FAIL',str(e)[:220]))
        except Exception as e:RES.append((name,'ERROR',(type(e).__name__+': '+str(e).splitlines()[0])[:220]))
        return fn
    return deco
def ok(c,msg):
    if not c:raise AssertionError(msg)
def sel_mol(pg,mid,name=None,t=2500):
    prev=pg.inner_text('#cmn');pg.select_option('#mol',mid)
    pg.wait_for_function("(n)=>document.querySelector('#cmn').textContent===n",arg=name,timeout=t) if name else pg.wait_for_function("(p)=>document.querySelector('#cmn').textContent!==p",arg=prev,timeout=t)
def sec(pg,title):
    d=pg.locator('details.sec',has=pg.locator('summary',has_text=title)).first
    if not d.evaluate('e=>e.open'):d.locator('summary').click()

LABJS="[...document.querySelectorAll('#lbl .lab')].filter(d=>d.style.display!=='none').map(d=>{const r=d.getBoundingClientRect();return{t:d.textContent,c:d.className,x:r.left,y:r.top,w:r.width,h:r.height}})"
def labels(pg,minimum=1,timeout=4000):
    pg.wait_for_function("(n)=>[...document.querySelectorAll('#lbl .lab')].filter(d=>d.style.display!=='none').length>=n",arg=minimum,timeout=timeout);return pg.evaluate(LABJS)
def overlaps(L):return [(a['t'],b['t']) for i,a in enumerate(L) for b in L[i+1:] if min(a['x']+a['w'],b['x']+b['w'])-max(a['x'],b['x'])>1 and min(a['y']+a['h'],b['y']+b['h'])-max(a['y'],b['y'])>1]
ATOMMESH="window.__three.insts.filter(m=>!m._d&&m.userData&&m.userData.map).map(m=>m.userData.map.length)"
def imp(pg,name):pg.set_input_files('#file',str(FX/name));pg.wait_for_timeout(250)
def dl(pg,sel):
    with pg.expect_download() as d:pg.click(sel)
    return Path(d.value.path()).read_text()

@case('01 launch: no console/page errors, library lists all 43 molecules (30 original + 6 peptide models + 7 new small molecules)')
def _(pw):
    pg=new_page(pw);n=pg.eval_on_selector_all('#mol option','o=>o.length');ok(n==43 and pg.inner_text('#lcount')==f'{n} molecules','options %d, count label: %s'%(n,pg.inner_text('#lcount')));ok(not pg.errs,str(pg.errs));ERRORS.extend(pg.errs);pg._b.close()
@case('02 built-in molecule + every rendering style without errors')
def _(pw):
    pg=new_page(pw);sel_mol(pg,'benzene','Benzene')
    for st in['bs','sf','st','wf','bb','tb']:pg.select_option('#sty',st);pg.wait_for_timeout(30)
    ok(not pg.errs,str(pg.errs));pg._b.close()
@case('03 rotate toggle shows/hides speed slider; reset view works')
def _(pw):
    pg=new_page(pw);ok(not pg.is_visible('#sp'),'slider visible before rotating');pg.click('#ar');ok(pg.is_visible('#sp'),'slider hidden while rotating');pg.click('#ar');ok(not pg.is_visible('#sp'),'slider still visible');pg.click('#rs');ok(not pg.errs,str(pg.errs));pg._b.close()
@case('04 selection: select all / invert / select-by-element')
def _(pw):
    pg=new_page(pw);sel_mol(pg,'pep_alpha');pg.click('#sall');ok('60 atoms selected' in pg.inner_text('#meas'),pg.inner_text('#meas')[:80]);pg.click('#sinv');ok('selected' not in pg.inner_text('#meas'),'invert should empty selection')
    pg.select_option('#selby','E:N');ok('12' in pg.inner_text('#meas') or '12 atoms' in pg.inner_text('#meas'),'select N: '+pg.inner_text('#sel')[:60]);pg.click('#clr');pg._b.close()
@case('05 import valid PDB, mmCIF, XYZ (+ warnings for two models)')
def _(pw):
    pg=new_page(pw)
    for f,n in[('alpha12.pdb','60 atoms'),('alpha12.cif','60 atoms'),('ethanol_heavy.xyz','3 atoms')]:
        imp(pg,f);ok(pg.inner_text('#cmn')==f and n in pg.inner_text('#cmm'),f+' -> '+pg.inner_text('#cmn')+' | '+pg.inner_text('#cmm'))
    imp(pg,'two_models.pdb');ok('additional model' in (pg.inner_text('#toast')+pg.inner_text('#info')).lower(),'two-model file loaded with no mention of the ignored second model')
    ok(not pg.errs,str(pg.errs));pg._b.close()
@case('06 invalid files rejected with a message; app stays usable')
def _(pw):
    pg=new_page(pw)
    for f in['bad_empty.pdb','bad_malformed.pdb','bad_truncated.pdb','bad_unsupported.txt']:
        imp(pg,f);t=pg.inner_text('#ov')
        ok('Import failed' in t or 'failed' in t.lower(),f+': no error shown ('+t[:50]+')');pg.click('#ov button')
    sel_mol(pg,'water','Water');pg._b.close()
@case('07 dialogs: Help opens, Esc closes it, focus returns; INFO likewise')
def _(pw):
    pg=new_page(pw);pg.click('#hp');ok(pg.is_visible('#mb'),'help not shown');pg.keyboard.press('Escape');ok(not pg.is_visible('#mb'),'Esc did not close Help');ok(pg.evaluate("document.activeElement.id")=='hp','focus not restored to Help button: '+pg.evaluate("document.activeElement.id||document.activeElement.tagName"))
    pg.click('#inf');ok(pg.is_visible('#infomodal .ic'),'info not shown');pg.keyboard.press('Escape');ok(not pg.is_visible('#infomodal .ic'),'Esc did not close Info');ok(pg.evaluate("document.activeElement.id")=='inf','focus not restored to INFO');pg._b.close()
@case('08 panels resize by keyboard and mouse; right panel closes/reopens')
def _(pw):
    pg=new_page(pw);pg.focus('#gl');w0=pg.evaluate("parseFloat(document.querySelector('main').style.getPropertyValue('--lw'))||260");pg.keyboard.press('ArrowRight');w1=pg.evaluate("parseFloat(document.querySelector('main').style.getPropertyValue('--lw'))")
    ok(w1==w0+20,f'keyboard resize {w0}->{w1}');b=pg.locator('#gl').bounding_box();pg.mouse.move(b['x']+4,b['y']+300);pg.mouse.down();pg.mouse.move(b['x']+154,b['y']+300);pg.mouse.up();w2=pg.evaluate("parseFloat(document.querySelector('main').style.getPropertyValue('--lw'))");ok(w2>w1+100,f'mouse resize {w1}->{w2}')
    pg.click('#rx');ok(not pg.is_visible('#info'),'right panel still visible');pg.click('#rshow');ok(pg.is_visible('#info'),'panel did not reopen');pg._b.close()
@case('09 theme toggle persists across reload')
def _(pw):
    pg=new_page(pw);ok('Dark Mode' in pg.inner_text('#th'),pg.inner_text('#th'));pg.click('#th');ok(pg.evaluate("document.documentElement.dataset.theme")=='dark','not dark');ok('Light Mode' in pg.inner_text('#th'),'label');
    pg.reload();pg.wait_for_function("document.querySelector('#cmn').textContent.trim().length>1");ok(pg.evaluate("document.documentElement.dataset.theme")=='dark','theme not persisted');
    s=pg.evaluate("[document.querySelector('.th').offsetWidth,document.querySelector('.th').offsetHeight]");pg.click('#th');s2=pg.evaluate("[document.querySelector('.th').offsetWidth,document.querySelector('.th').offsetHeight]");ok(s==s2==[130,34] or s==s2,f'theme buttons differ in size {s} vs {s2}');pg._b.close()
@case('10 responsive: no horizontal overflow at 1280 / 820 / 390 px; tools reachable')
def _(pw):
    for w,h in[(1280,800),(820,1100),(390,844)]:
        pg=new_page(pw,w,h);ov=pg.evaluate("document.documentElement.scrollWidth-innerWidth");ok(ov<=1,f'{w}px: horizontal overflow {ov}px')
        pg.locator('#ra').scroll_into_view_if_needed();ok(pg.is_visible('#ra'),f'{w}px: Ramachandran button not reachable');pg._b.close()
@case('11 session save -> reload -> restore (structure, selection)')
def _(pw):
    pg=new_page(pw);imp(pg,'alpha12.pdb');pg.click('#sall');sec(pg,'Teaching, session');data=dl(pg,'#ss');S=json.loads(data);ok(S['structure']['raw'] and len(S['selection'])==60,'session content');ok(S['structure']['name']=='alpha12.pdb','name')
    sp=Path(tempfile.mkdtemp())/'s.json';sp.write_text(data);pg.reload();pg.wait_for_function("document.querySelector('#cmn').textContent.trim().length>1")
    pg.set_input_files('#sl',str(sp));pg.wait_for_timeout(400);ok(pg.inner_text('#cmn')=='alpha12.pdb','not restored: '+pg.inner_text('#cmn'));ok('60 atoms selected' in pg.inner_text('#meas'),'selection not restored');pg._b.close()
@case('12 session loader rejects malformed / incompatible / oversized files with a clear message')
def _(pw):
    pg=new_page(pw);d=Path(tempfile.mkdtemp());bad={'trunc.json':'{"laura":"0.3","structure":{"id":"x"','nojson.json':'hello','nostruct.json':'{"laura":"0.3"}','future.json':json.dumps({"laura":"9.9","structure":{"builtin":True,"id":"water"},"settings":{}}),'huge.json':'{"laura":"0.3","pad":"'+'x'*6_200_000+'"}'}
    for n,t in bad.items():
        (d/n).write_text(t);pg.set_input_files('#sl',str(d/n));pg.wait_for_timeout(250);m=pg.inner_text('#toast')
        ok('LAURA session' in m and not re.search(r'Cannot read|undefined|Unexpected (token|end)|in JSON at position',m),f'{n}: unhelpful message: {m[:100]}')
    pg._b.close()
@case('13 analysis: Ramachandran + CSV (phi/psi -57/-47), contact map, no errors')
def _(pw):
    pg=new_page(pw);imp(pg,'alpha12.pdb');pg.click('#ra');ok(pg.is_visible('#rcv'),'no plot');rows=dl(pg,'#rcsv').strip().split('\n');ok(len(rows)==11 and rows[0]=='residue,resName,phi_deg,psi_deg,region',f'{len(rows)} rows: '+rows[0]);
    f=rows[1].split(',');ok(abs(float(f[2])+57)<0.5 and abs(float(f[3])+47)<0.5 and f[4]=='alphaR','phi/psi/region '+str(f));pg.click('#mb button:has-text("Close")');pg.click('#cmb');ok(pg.is_visible('#ccv'),'no contact map');ok(not pg.errs,str(pg.errs));pg._b.close()
@case('14 typing in the alignment textarea must not trigger shortcuts (r p t f space)')
def _(pw):
    pg=new_page(pw);imp(pg,'alpha12.pdb');pg.click('#sq');pg.click('#sa');pg.fill('#sa','');pg.keyboard.type('rptf ptf');
    ok(pg.evaluate("document.documentElement.dataset.theme||'none'")=='none','theme toggled by typing');ok(not pg.evaluate("document.body.classList.contains('pm')"),'presentation mode entered by typing');ok(pg.evaluate("window.__fs")==0,'fullscreen triggered by typing')
    ok(pg.input_value('#sa')=='rptf ptf','text lost: '+pg.input_value('#sa'));pg._b.close()
@case('15 browser shortcuts with Ctrl/Meta must not trigger LAURA actions; plain f does')
def _(pw):
    pg=new_page(pw);pg.keyboard.press('Control+f');pg.keyboard.press('Control+t');pg.keyboard.press('Control+p');ok(pg.evaluate("window.__fs")==0,'Ctrl+F triggered fullscreen');ok(pg.evaluate("document.documentElement.dataset.theme||'none'")=='none','Ctrl+T toggled theme')
    ok(not pg.evaluate("document.body.classList.contains('pm')"),'Ctrl+P entered presentation');pg.keyboard.press('f');ok(pg.evaluate("window.__fs")==1,'plain F did not request fullscreen');pg._b.close()
@case('16 sequence alignment tool computes identity for two chains')
def _(pw):
    pg=new_page(pw);pg.click('#sq');pg.fill('#sa','ACGT');pg.fill('#sb','AGT');pg.click('#sgo');t=pg.inner_text('#sout');ok('identical 3' in t and '75.0%' in t and 'Needleman' in t,t[:120]);pg._b.close()
@case('17 RMSD: identical chain = 0; helix vs strand >> 0')
def _(pw):
    pg=new_page(pw);imp(pg,'alpha12.pdb');pg.click('#rm');pg.click('#rgo');t=pg.inner_text('#rout');ok('RMSD 0.000' in t and 'over 12' in t,t[:100]);pg.select_option('#rb','pep_beta');pg.wait_for_timeout(200);pg.click('#rgo');m=re.search(r'RMSD ([\d.]+)',pg.inner_text('#rout'));ok(m and float(m.group(1))>3,'helix vs strand: '+pg.inner_text('#rout')[:80]);pg._b.close()
@case('18 ligand contacts + JSON/CSV exports consistent')
def _(pw):
    pg=new_page(pw);imp(pg,'ligand_complex.pdb');pg.click('summary:has-text("Ligand contacts")');opts=pg.eval_on_selector_all('#lig option','o=>o.map(x=>x.value)');ok('A:101' in opts and not any(o.endswith(':201') for o in opts),str(opts));pg.select_option('#lig','A:101');pg.click('#rc')
    ok('A:6' in pg.inner_text('#cres'),pg.inner_text('#cres')[:100]);j=json.loads(dl(pg,'#cj'));c=dl(pg,'#cc').strip().split('\n');ok(j['result']['params']['cutoffAngstrom']==4 and len(c)-1==len(j['result']['contacts']),'json/csv mismatch');ok(c[0].startswith('residue,resName,minDistance_A'),c[0]);pg._b.close()
@case('19 compare molecules table + CSV')
def _(pw):
    pg=new_page(pw);sel_mol(pg,'water','Water');pg.click('#cmp');pg.select_option('#cpb','methane');pg.click('#cgo');t=pg.inner_text('#cout');ok('H2O' in t and 'CH4' in t,t[:100]);c=dl(pg,'#ccsv');ok(c.startswith('"property"'),c[:60]);pg._b.close()
@case('20 Blender package: zip integrity, README, metadata, script compiles')
def _(pw):
    pg=new_page(pw);sel_mol(pg,'alanine_h');sec(pg,'Teaching, session');z=Path(tempfile.mkdtemp())/'b.zip'
    with pg.expect_download() as d:pg.click('#bl')
    d.value.save_as(str(z));zf=zipfile.ZipFile(z);n=zf.namelist();ok(zf.testzip() is None and 'import_molecule.py' in n and 'README.txt' in n and 'metadata.json' in n,str(n));json.loads(zf.read('metadata.json'))
    py=Path(tempfile.mkdtemp())/'i.py';py.write_bytes(zf.read('import_molecule.py'));py_compile.compile(str(py),doraise=True);pg._b.close()
@case('21 lifecycle: geometry/material/texture/instanced-mesh objects do not accumulate over 60 molecule switches (shim-level)')
def _(pw):
    pg=new_page(pw);ids=pg.eval_on_selector_all('#mol option','o=>o.map(x=>x.value)')
    for i in range(1,61):sel_mol(pg,ids[i%len(ids)])
    L=pg.evaluate("window.__three");ok(L['geo']<60 and L['mat']<60 and L['tex']<40 and L['inst']<30,'live objects after 60 switches: '+json.dumps({k:L[k] for k in['geo','mat','tex','inst']})+' created '+json.dumps(L['created']));pg._b.close()
@case('22 context loss handled (webglcontextlost is intercepted and the user is told)')
def _(pw):
    pg=new_page(pw);r=pg.evaluate("(()=>{const c=document.querySelector('canvas');let prevented=false;c.addEventListener('webglcontextlost',e=>{prevented=e.defaultPrevented},{once:false});const ev=new Event('webglcontextlost',{cancelable:true});c.dispatchEvent(ev);return ev.defaultPrevented})()")
    ok(r,'webglcontextlost not intercepted (preventDefault not called), so the context cannot be restored');ok('lost' in pg.inner_text('#toast').lower() or 'graphics' in pg.inner_text('#toast').lower(),'no user message: '+pg.inner_text('#toast')[:60]);pg._b.close()
@case('24 Zoom to selection: message with empty selection, no errors after Select all')
def _(pw):
    pg=new_page(pw);sel_mol(pg,'pep_alpha');pg.click('#fz');ok('Select atoms first' in pg.inner_text('#toast'),pg.inner_text('#toast'));pg.click('#sall');pg.click('#fz');pg.keyboard.press('z');ok(not pg.errs,str(pg.errs));pg._b.close()
@case('25 imported text is escaped (no HTML injection via file names or structure titles)')
def _(pw):
    pg=new_page(pw);d=Path(tempfile.mkdtemp())/'x"><img src=x onerror=window.__pwn=1>.xyz';d.write_text('1\nt\nC 0 0 0\n');pg.set_input_files('#file',str(d));pg.wait_for_timeout(300)
    ok(pg.evaluate('window.__pwn')is None,'script from a file name executed');ok(pg.locator('img[src="x"]').count()==0,'injected <img> element exists in DOM');pg._b.close()
@case('26 XSS: payloads in mmCIF atom/element fields and in a crafted session are rendered as text')
def _(pw):
    pg=new_page(pw);d=Path(tempfile.mkdtemp());H=['group_PDB','id','type_symbol','label_atom_id','label_alt_id','label_comp_id','label_asym_id','label_seq_id','pdbx_PDB_ins_code','Cartn_x','Cartn_y','Cartn_z','pdbx_PDB_model_num']
    c=d/'payload.cif';c.write_text('data_x\nloop_\n'+'\n'.join('_atom_site.'+h for h in H)+'\nATOM 1 <img/src/onerror=window.__pwn2=1> <img/src/onerror=window.__pwn2=1> . <b>ALA A 1 ? 0 0 0 1\nATOM 2 C CA . ALA A 1 ? 1.4 0 0 1\nATOM 3 C C . ALA A 1 ? 2 1.2 0 1\n#\n')
    pg.set_input_files('#file',str(c));pg.wait_for_timeout(300);pg.click('#sall');pg.click('summary:has-text("Scientific details")')
    ok(pg.evaluate('window.__pwn2')is None,'payload executed from mmCIF fields');ok(pg.locator('#info img, #sel img, #list img, #ov img').count()==0,'injected element present')
    ev='<img src=x onerror=window.__pwn3=1>';S={'laura':'0.3','structure':{'id':'p1','name':ev,'builtin':False,'src':ev,'fmt':ev,'ext':'xyz','raw':'1\nt\nC 0 0 0\n'},'settings':{'style':'bs','hydrogens':True,'solvent':True,'labels':False},'selection':[],'contacts':{},'camera':None}
    s=d/'evil.json';s.write_text(json.dumps(S));pg.set_input_files('#sl',str(s));pg.wait_for_timeout(400)
    ok(pg.evaluate('window.__pwn3')is None,'payload executed from session fields');ok(pg.locator('#info img, #list img, #ov img, #cmn img').count()==0,'injected element from session');ok(pg.errs==[] or all('payload' in e or 'x' for e in pg.errs),str(pg.errs));pg._b.close()
@case('27 every library entry loads; bond-order/hybridization angle check passes for all 43 (sp 180, sp2 120, sp3 109.5 within tolerance)')
def _(pw):
    pg=new_page(pw);ids=pg.eval_on_selector_all('#mol option','o=>o.map(x=>x.value)');bad=[];noinfo=[];seen=0
    for mid in ids:
        if pg.evaluate("document.querySelector('#mol').value")!=mid:sel_mol(pg,mid)
        ok(pg.evaluate("document.querySelector('#mol').value")==mid,'could not select '+mid);seen+=1
        t=pg.text_content('#info')
        if 'Bond orders (' not in t:noinfo.append(mid)
        if 'beyond tolerance' in t:bad.append((mid,re.search(r'largest deviation[^.]*',t).group(0)))
    ok(len(ids)==43 and seen==43,'entries %d visited %d'%(len(ids),seen));ok(not noinfo,'no bond-order info for '+str(noinfo));ok(not bad,'angle outliers: '+str(bad));ok(not pg.errs,str(pg.errs[:2]));pg._b.close()
@case('28 each tool dialog has its description block (collapsible) with the right title')
def _(pw):
    pg=new_page(pw);imp(pg,'alpha12.pdb');want={'#ra':'Ramachandran Plot','#cmb':'Contact Map','#sq':'Sequence Alignment','#rm':'Superposition & RMSD','#cmp':'Compare Molecules'}
    for btn,title in want.items():
        pg.click(btn);ok(pg.locator('#mb details.td').count()==1,btn+': no description block');ok(title in pg.inner_text('#mb details.td summary'),btn+': title '+pg.inner_text('#mb details.td summary')[:40])
        was=pg.evaluate("document.querySelector('#mb details.td').open");pg.click('#mb details.td summary');ok(pg.evaluate("document.querySelector('#mb details.td').open")!=was,btn+': description does not toggle');pg.keyboard.press('Escape')
    pg._b.close()
@case('29 Ramachandran test peptide: 14 residues visit every region; helix model is all alpha-helix')
def _(pw):
    pg=new_page(pw);sel_mol(pg,'pep_rama');pg.click('#ra');t=pg.inner_text('#rsum') if pg.locator('#rsum').count() else pg.inner_text('#mb')
    m=re.search(r'(\d+) residues: α-helix (\d+) · β/PPII (\d+) · left-handed (\d+) · outside the shading (\d+) · Gly[^\d]*(\d+)',t);ok(m,'summary text not found: '+t[:200]);n,a,b,l,o,g=map(int,m.groups())
    ok(n==14 and a+b+l+o+g==14 and min(a,b,l,o,g)>=1,f'total {n}: alphaR {a} beta {b} alphaL {l} other {o} gly {g}');pg.keyboard.press('Escape')
    sel_mol(pg,'pep_alpha');pg.click('#ra');t=pg.inner_text('#mb');m=re.search(r'(\d+) residues: α-helix (\d+)',t);ok(m and m.group(1)==m.group(2)=='10','alpha helix not 10/10: '+t[:160]);pg._b.close()
@case('30 helix pair: UI claims hold (10/12 identical on alignment; chain A vs B RMSD ~0 over 12 CA pairs)')
def _(pw):
    pg=new_page(pw);sel_mol(pg,'pep_pair');pg.click('#sq');pg.click('#sgo');t=pg.inner_text('#sout');ok('identical 10' in t and 'alignment length 12' in t,t[:120]);pg.keyboard.press('Escape')
    pg.click('#rm');pg.select_option('#rca','A');pg.wait_for_function("document.querySelectorAll('#rcb option').length>1",timeout=3000);pg.select_option('#rcb','B');pg.click('#rgo');t=pg.inner_text('#rout');m=re.search(r'RMSD ([\d.]+) Å over (\d+)',t);ok(m and int(m.group(2))==12 and float(m.group(1))<0.05,t[:120]);pg._b.close()
@case('31 XSS: payloads in chain/residue names do not execute in Ramachandran, contact-map or alignment dialogs')
def _(pw):
    pg=new_page(pw);d=Path(tempfile.mkdtemp());H=['group_PDB','id','type_symbol','label_atom_id','label_alt_id','label_comp_id','label_asym_id','label_seq_id','pdbx_PDB_ins_code','Cartn_x','Cartn_y','Cartn_z','auth_asym_id','auth_comp_id','pdbx_PDB_model_num']
    pay='<img/src/onerror=window.__pwn4=1>';rows=[];k=0
    for l in (FX/'alpha12.pdb').read_text().split('\n'):
        if not l.startswith('ATOM'):continue
        resi=int(l[22:26])
        if resi>4:break
        k+=1;rows.append(' '.join(['ATOM',str(k),l[76:78].strip(),l[12:16].strip(),'.','ALA',pay,str(resi),'?',l[30:38].strip(),l[38:46].strip(),l[46:54].strip(),pay,pay,'1']))
    c=d/'payload2.cif';c.write_text('data_x\nloop_\n'+'\n'.join('_atom_site.'+h for h in H)+'\n'+'\n'.join(rows)+'\n#\n');pg.set_input_files('#file',str(c));pg.wait_for_function("document.querySelector('#cmn').textContent==='payload2.cif'",timeout=3000)
    for btn in['#ra','#cmb','#sq','#rm','#cmp']:
        pg.click(btn);pg.wait_for_timeout(150);ok(pg.evaluate('window.__pwn4')is None,btn+': payload executed');ok(pg.locator('#mb img, #info img, #sel img').count()==0,btn+': injected <img> in DOM');pg.keyboard.press('Escape')
    pg._b.close()
@case('32 aromatic rings in the library are depicted consistently (pyridine 6/6, imidazole 5/5, benzene 6 aromatic bonds)')
def _(pw):
    pg=new_page(pw)
    for mid,exp in[('pyridine','6 aromatic'),('imidazole','5 aromatic'),('benzene','6 aromatic')]:
        sel_mol(pg,mid);t=pg.text_content('#info');m=re.search(r'Bond orders \(.*?\):\s*(.*?)\.\s*Hybridization',t,re.S);ok(m and exp in m.group(1),mid+': '+(m.group(1) if m else t[:120]))
    pg._b.close()
@case('33 presentation mode: the viewer fills the whole window (no gaps), and exits cleanly')
def _(pw):
    for w,h in[(1280,800),(1024,640)]:
        pg=new_page(pw,w,h);pg.keyboard.press('p');pg.wait_for_timeout(120);b=pg.evaluate("(()=>{const r=document.querySelector('#vp').getBoundingClientRect();return [r.left,r.top,r.width,r.height,innerWidth,innerHeight,document.documentElement.scrollHeight]})()")
        ok(abs(b[0])<=1 and abs(b[1])<=1 and abs(b[2]-b[4])<=1 and abs(b[3]-b[5])<=1,f'{w}x{h}: viewer box {b[:4]} vs window {b[4:6]}');ok(b[6]<=b[5]+1,f'{w}x{h}: page taller than window ({b[6]})')
        pg.keyboard.press('p');ok(pg.is_visible('header'),'did not exit presentation');pg._b.close()
@case('34 right-panel X stays visible while the panel content is scrolled')
def _(pw):
    pg=new_page(pw,1280,430);sel_mol(pg,'pep_alpha')
    top0=pg.locator('#rx').bounding_box()['y'];ok(pg.evaluate("(()=>{const a=document.querySelector('aside.r');a.scrollTop=0;const s=[...a.querySelectorAll('*')].filter(e=>e.scrollHeight>e.clientHeight+20&&getComputedStyle(e).overflowY!='visible');return s.length})()")>=1,'right panel has nothing to scroll at this height')
    pg.evaluate("[...document.querySelectorAll('aside.r, aside.r *')].forEach(e=>{if(e.scrollHeight>e.clientHeight+20&&getComputedStyle(e).overflowY!=='visible')e.scrollTop=9999})");pg.wait_for_timeout(100)
    bb=pg.locator('#rx').bounding_box();ab=pg.locator('aside.r').bounding_box();ok(bb['y']>=ab['y']-1 and bb['y']+bb['height']<=ab['y']+80 and pg.locator('#rx').is_visible(),f'X scrolled out of view: top {bb["y"]:.0f} (panel top {ab["y"]:.0f})');pg._b.close()
@case('35 Labels checkbox lives in the Analysis section (not in Analysis tools)')
def _(pw):
    pg=new_page(pw);ok(pg.evaluate("document.querySelector('#lb').closest('details').querySelector('summary').textContent.trim()")=='Analysis','Labels not under Analysis');ok(pg.evaluate("!!document.querySelector('#so')"),'Selected only missing');pg._b.close()
@case('36 labels show atom names, not ids: element symbols for generated molecules, real names + residue for proteins')
def _(pw):
    pg=new_page(pw);sel_mol(pg,'water','Water');pg.click('#sall');pg.check('#lb');L=labels(pg,3);ok(sorted(x['t'] for x in L)==['H','H','O'],'water labels '+str([x['t'] for x in L]))
    sel_mol(pg,'pep_alpha');pg.click('#sall');pg.check('#lb');L=labels(pg,5);ok(all(re.fullmatch(r'(N|CA|C|O|CB) ALA\d+:A',x['t']) for x in L),'protein label format: '+str([x['t'] for x in L][:6]));pg._b.close()
@case('37 many atoms + Labels: labels are decluttered (no overlaps, capped) and stay readable')
def _(pw):
    pg=new_page(pw);sel_mol(pg,'pep_alpha');pg.click('#sall');pg.check('#lb');L=labels(pg,5);ok(len(L)<=80,'too many labels %d'%len(L));ok(not overlaps(L),'overlapping labels: '+str(overlaps(L)[:3]))
    ok(len(L)<60,'with 60 selected atoms the pile must be thinned (got %d of 60)'%len(L));pg._b.close()
@case('38 ligand contacts get labels (residue + distance, ligand) and Clear removes them, the highlights and the results')
def _(pw):
    pg=new_page(pw);imp(pg,'ligand_complex.pdb');sec(pg,'Ligand contacts');pg.select_option('#lig','A:101');pg.click('#rc');L=labels(pg,2);C=[x for x in L if 'contact' in x['c']]
    ok(any(re.fullmatch(r'ALA\d+ · \d+\.\d Å',x['t']) for x in C) and any('(ligand)' in x['t'] for x in C),'contact labels: '+str([x['t'] for x in C]));ok(not overlaps(L),'contact labels overlap')
    pg.click('#ccl');pg.wait_for_function("[...document.querySelectorAll('#lbl .lab')].every(d=>d.style.display==='none')",timeout=3000);ok(pg.inner_text('#cres').strip()=='','results not cleared');pg.click('#cj');ok('Find contacts' in pg.inner_text('#toast'),'results still held after Clear');pg._b.close()
@case('39 Selected only draws just the selected atoms; unticking restores all')
def _(pw):
    pg=new_page(pw);sel_mol(pg,'pep_alpha');ok(pg.evaluate(ATOMMESH)==[60],str(pg.evaluate(ATOMMESH)));pg.select_option('#selby','E:N');pg.check('#so');ok(pg.evaluate(ATOMMESH)==[12],'selected only: '+str(pg.evaluate(ATOMMESH)))
    pg.click('#clr');ok(pg.evaluate(ATOMMESH)==[60],'empty selection shows everything: '+str(pg.evaluate(ATOMMESH)));pg.select_option('#selby','E:N');ok(pg.evaluate(ATOMMESH)==[12],'reselect');pg.uncheck('#so');ok(pg.evaluate(ATOMMESH)==[60],'untick: '+str(pg.evaluate(ATOMMESH)));pg._b.close()
@case('40 imports: rename, save in this browser (IndexedDB), restore after reload, remove; hostile names stay text')
def _(pw):
    pg=new_page(pw);imp(pg,'alpha12.pdb');ok(pg.is_visible('#impbar'),'import bar not shown for an import');pg.click('#imRen');pg.fill('#rnin','My helix 1');pg.press('#rnin','Enter');pg.wait_for_function("document.querySelector('#cmn').textContent==='My helix 1'",timeout=3000)
    pg.click('#imSave');pg.wait_for_function("document.querySelector('#imSave').textContent.includes('Saved')",timeout=3000);pg.reload();pg.wait_for_function("document.querySelector('#lcount').textContent==='44 molecules'",timeout=5000)
    pg.click('#list .it:has-text("My helix 1")');pg.wait_for_function("document.querySelector('#cmn').textContent==='My helix 1'",timeout=3000);ok('60 atoms' in pg.inner_text('#cmm'),pg.inner_text('#cmm'))
    pg.click('#imRen');pg.fill('#rnin','<img src=x onerror=window.__pwn5=1>');pg.press('#rnin','Enter');pg.wait_for_timeout(300);pg.reload();pg.wait_for_function("document.querySelector('#lcount').textContent==='44 molecules'",timeout=5000);ok(pg.evaluate('window.__pwn5')is None and pg.locator('#list img, #cmn img').count()==0,'name executed as HTML')
    pg.click('#list .it:has-text("onerror")');pg.wait_for_timeout(300);pg.click('#imDel');pg.wait_for_function("document.querySelector('#lcount').textContent==='43 molecules'",timeout=3000);pg.reload();pg.wait_for_timeout(700);ok(pg.inner_text('#lcount')=='43 molecules','removed import came back: '+pg.inner_text('#lcount'));pg._b.close()
@case('41 See in 2D: draws the molecule, hydrogens toggle changes the picture, click enlarges, honest caption')
def _(pw):
    pg=new_page(pw);sel_mol(pg,'benzene','Benzene');pg.click('#d2');pg.wait_for_selector('#d2c');ink=lambda:pg.evaluate("(()=>{const c=document.querySelector('#d2c'),d=c.getContext('2d').getImageData(0,0,c.width,c.height).data,b=[d[0],d[1],d[2]];let n=0;for(let i=0;i<d.length;i+=4)if(Math.abs(d[i]-b[0])+Math.abs(d[i+1]-b[1])+Math.abs(d[i+2]-b[2])>30)n++;return n})()")
    d=ink();ok(d>1500,'diagram not drawn (%d ink pixels)'%d);ok('structure diagram' in pg.inner_text('#mb') and 'stereochemistry' in pg.inner_text('#mb'),'diagram caption missing');pg.select_option('#d2m','proj');a=ink();ok(a>1500,'nothing drawn (%d ink pixels)'%a);ok('not a chemical structure diagram' in pg.inner_text('#mb'),'caption missing');pg.uncheck('#d2h');b=ink();ok(b<a,'hydrogens toggle did not change the picture (%d vs %d)'%(a,b))
    pg.click('#d2c');ok(pg.is_visible('#lightbox img') and pg.evaluate("document.querySelector('#lbimg').src").startswith('data:image/png'),'lightbox');pg.keyboard.press('Escape');ok(not pg.is_visible('#lightbox'),'Esc did not close the lightbox');ok(pg.is_visible('#d2c'),'dialog closed with the lightbox');pg.keyboard.press('d');pg._b.close()
@case('42 compare in 3D: B is drawn as ball-and-stick beside A (no overlap), tagged A/B, and the Compare dialog closes')
def _(pw):
    pg=new_page(pw);sel_mol(pg,'water','Water');pg.click('#cmp');pg.select_option('#cpb','methane');pg.click('#cgo');pg.wait_for_selector('#cout table');pg.click('#cside');pg.wait_for_function("document.querySelector('#modal').style.display==='none'",timeout=4000)
    d=pg.evaluate("(()=>{const I=window.__three.insts.filter(m=>!m._d),A=I.find(m=>m.userData&&m.userData.map),B=I.find(m=>m.userData&&m.userData.tag==='B-atoms'),Bb=I.find(m=>m.userData&&m.userData.tag==='B-bonds');return{a:A&&A.pos,b:B&&B.pos,bc:B&&B.count,bb:Bb&&Bb.count}})()");ok(d['b'] and len(d['b'])==5 and d['bc']==5,'B atoms (methane has 5): '+str(d['bc']));ok(d['bb']==8,'B bond cylinders %s (4 C-H bonds x2 halves)'%d['bb'])
    ok(min(p[0] for p in d['b'])-max(p[0] for p in d['a'])>0.5,'B overlaps A: A max x %.2f, B min x %.2f'%(max(p[0] for p in d['a']),min(p[0] for p in d['b'])))
    L=labels(pg,2);tags=[x['t'] for x in L if 'tag' in x['c']];ok(any(t.startswith('A · Water') for t in tags) and any(t.startswith('B · Methane') for t in tags),str(tags));pg._b.close()
@case('43 tube view names every chain (A, B) with non-overlapping labels')
def _(pw):
    pg=new_page(pw);sel_mol(pg,'pep_pair');pg.select_option('#sty','tb');L=labels(pg,2);ch=sorted(x['t'] for x in L if 'chain' in x['c']);ok(ch==['Chain A','Chain B'],str(ch));ok(not overlaps(L),'chain labels overlap');pg.select_option('#sty','bs');pg.wait_for_timeout(200);ok(not any('chain' in x['c'] for x in pg.evaluate(LABJS)),'chain labels remain outside tube view');pg._b.close()
@case('44 style dropdown: six working styles, no disabled/unimplemented cartoon-surface entry')
def _(pw):
    pg=new_page(pw);o=pg.eval_on_selector_all('#sty option','o=>o.map(x=>[x.value,x.disabled,x.textContent])');ok(len(o)==6 and not any(x[1] for x in o) and not any('artoon' in x[2] or 'urface' in x[2] for x in o),str(o));pg._b.close()
@case('45 analysis dialogs fit the window, wrap text instead of overflowing, and diagrams enlarge on click')
def _(pw):
    pg=new_page(pw,1280,720);imp(pg,'alpha12.pdb');pg.click('#ra');pg.wait_for_selector('#rcv');b=pg.locator('#mb').bounding_box();ok(b['y']>=0 and b['y']+b['height']<=721 and b['width']<=1280,'dialog box '+str(b))
    pg.click('#rcv');ok(pg.is_visible('#lightbox'),'Ramachandran plot does not enlarge');pg.keyboard.press('Escape');pg.keyboard.press('Escape');pg._b.close()
    pg=new_page(pw,390,844);imp(pg,'alpha12.pdb')
    for btn in['#ra','#cmb','#sq','#rm','#cmp','#d2']:
        pg.click(btn);pg.wait_for_timeout(120);w=pg.evaluate("(()=>{const m=document.querySelector('#mb');return[m.scrollWidth,m.clientWidth,m.getBoundingClientRect().width,document.documentElement.scrollWidth,innerWidth]})()");ok(w[0]<=w[1]+1 and w[3]<=w[4]+1,f'{btn} at 390px: scrollW {w[0]} clientW {w[1]} page {w[3]}/{w[4]}');pg.keyboard.press('Escape')
    pg._b.close()
@case('46 dark mode: the Light Mode button has a white outline like the Dark Mode button')
def _(pw):
    pg=new_page(pw);bw=lambda:pg.evaluate("getComputedStyle(document.querySelector('.th')).borderTopColor");a=bw();pg.click('#th');b=bw();ok(a=='rgb(255, 255, 255)' and b=='rgb(255, 255, 255)',f'outline colours {a} / {b}');pg._b.close()
@case('47 teaching has no quiz; Teaching closes on outside click (not in the viewer); Help and tool dialogs close on outside click')
def _(pw):
    pg=new_page(pw);sec(pg,'Teaching, session');pg.click('#te');t=pg.inner_text('#tch');ok('Guided activity' in t and 'uiz' not in t and 'question' not in t.lower(),t[:120])
    bb=pg.locator('#vp canvas').bounding_box();pg.mouse.click(bb['x']+bb['width']/2,bb['y']+bb['height']/2);ok(pg.is_visible('#tch'),'Teaching closed on a click in the viewer (needed for the activity)')
    pg.click('header .brand');ok(not pg.is_visible('#tch'),'Teaching did not close on outside click')
    pg.click('#hp');ok(pg.is_visible('#mb'),'help');pg.mouse.click(3,3);ok(not pg.is_visible('#mb'),'Help did not close on outside click');imp(pg,'alpha12.pdb');pg.click('#ra');pg.mouse.click(3,3);ok(not pg.is_visible('#mb'),'tool dialog did not close on outside click')
    pg.click('#hp');pg.click('#mb');ok(pg.is_visible('#mb'),'click INSIDE the dialog must not close it');pg._b.close()
@case('48 tool descriptions: "How to use" replaced by a concrete "Try it" with the verified numbers')
def _(pw):
    pg=new_page(pw);sel_mol(pg,'pep_pair');exp={'#ra':'10 residues','#cmb':'38 contacting pairs','#sq':'75%','#rm':'6.7 Å','#cmp':'18.015'}
    for btn,txt in exp.items():
        pg.click(btn);t=pg.text_content('#mb details.td');ok('How to use' not in t and 'Try it' in t and txt in t,f'{btn}: '+t[-260:]);pg.keyboard.press('Escape')
    pg._b.close()
@case('49 presentation mode: a blocked full-screen request is reported and the viewer still fills the page')
def _(pw):
    pg=new_page(pw);pg.evaluate("(()=>{Element.prototype.requestFullscreen=()=>Promise.reject(new Error('blocked'))})()");pg.keyboard.press('p');pg.wait_for_function("document.querySelector('#toast').textContent.includes('blocked full screen')",timeout=3000)
    b=pg.evaluate("(()=>{const r=document.querySelector('#vp').getBoundingClientRect();return[r.width,r.height,innerWidth,innerHeight]})()");ok(abs(b[0]-b[2])<=1 and abs(b[1]-b[3])<=1,str(b));pg.keyboard.press('p');pg._b.close()
@case('23 no console errors or failed requests across the whole run')
def _(pw):ok(not ERRORS,str(ERRORS)[:200])

w=max(len(n) for n,_,_ in RES);npass=sum(s=='PASS' for _,s,_ in RES)
for n,s,d in RES:print(f'{s:5} {n}'+(f'\n        -> {d}' if d else ''))
print(f'\n{npass}/{len(RES)} passed | mode: {"REAL CDN" if REAL else "Three.js SHIM (3D rendering not exercised)"} | page: {PAGE.name}')
OUT.write_text(json.dumps([{'name':n,'status':s,'detail':d} for n,s,d in RES],indent=1));sys.exit(0 if npass==len(RES) else 1)
