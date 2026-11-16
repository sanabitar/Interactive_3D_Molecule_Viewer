import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs';
// Built-in molecule verification. The library lives inside laura-viewer.html, so it is evaluated here (no DOM needed).
// The expected chemistry comes from SMILES strings parsed by an independent mini-parser (implicit hydrogens from standard valences),
// NOT from the code that builds the molecules. Each atom is compared by its local environment: element + sorted neighbour elements.
const html=fs.readFileSync(new URL('../laura-viewer.html',import.meta.url),'utf8');
const a=html.indexOf('/*HYDROGENS-BEGIN*/'),b=html.indexOf('let imported=[];'),c=html.indexOf('const Core=(()=>{'),d=html.indexOf('})();',html.indexOf('return{parsePDB,'))+5,
 cb=html.indexOf('function computeBonds(A)'),ce=html.indexOf('// ---------- measurements ----------'),h0=html.indexOf('const CPK={'),h1=html.indexOf('function toast(');
const {LIB,Core,computeBonds}=new Function(html.slice(h0,h1)+'\n'+html.slice(a,b)+'\n'+html.slice(c,d)+'\n'+html.slice(cb,ce)+'\nreturn {LIB,Core,computeBonds};')();

function parseSmiles(s){const atoms=[],bonds=[],stack=[],ring={},VAL={C:4,N:3,O:2,S:2};let prev=null,ord=null;
 for(const t of s.match(/\[[^\]]+\]|Cl|Br|[BCNOPSFI]|[cnos]|\(|\)|=|#|\d/g)){
  if(t==='('){stack.push(prev);continue}if(t===')'){prev=stack.pop();continue}if(t==='='){ord=2;continue}if(t==='#'){ord=3;continue}
  if(/^\d$/.test(t)){if(ring[t]){bonds.push([ring[t].i,prev,ring[t].o||ord||(atoms[prev].ar&&atoms[ring[t].i].ar?'a':1)]);delete ring[t]}else ring[t]={i:prev,o:ord};ord=null;continue}
  let el,ar,H=null,q=0;if(t[0]==='['){const m=t.match(/^\[([A-Z][a-z]?|[cnos])(H\d?)?([+-])?\]$/);el=m[1];H=m[2]?(+m[2].slice(1)||1):0;q=m[3]==='+'?1:m[3]==='-'?-1:0}else el=t;
  ar=/^[cnos]$/.test(el);el=el.toUpperCase().length&&ar?el.toUpperCase():el;const i=atoms.length;atoms.push({el,ar,H,q});
  if(prev!==null)bonds.push([prev,i,ord||(ar&&atoms[prev].ar?'a':1)]);prev=i;ord=null}
 atoms.forEach((x,i)=>{if(x.H===null){const used=bonds.filter(b=>b[0]===i||b[1]===i).reduce((s,b)=>s+(b[2]==='a'?1:b[2]),0)+(x.ar?1:0);x.H=Math.max(0,VAL[x.el]-used)}});
 const n0=atoms.length;atoms.forEach((x,i)=>{for(let k=0;k<x.H;k++){bonds.push([i,atoms.length,1]);atoms.push({el:'H'})}});return{atoms,bonds:bonds.map(([i,j])=>[i,j])}}
const sig=(els,bonds)=>{const nb=els.map(()=>[]);bonds.forEach(([i,j])=>{nb[i].push(els[j]);nb[j].push(els[i])});return els.map((e,i)=>e+'('+nb[i].sort().join('')+')').sort()};
const SMILES={water:'O',methane:'C',benzene:'c1ccccc1',co2:'O=C=O',ammonia:'N',ethene:'C=C',ethyne:'C#C',ethane:'CC',
 glycine:'NCC(=O)O',alanine:'NC(C)C(=O)O',alanine_h:'[NH3+]C(C)C(=O)[O-]',ethanol:'CCO',acetic:'CC(=O)O',phenol:'c1ccccc1O',toluene:'Cc1ccccc1',aniline:'Nc1ccccc1',
 hcn:'C#N',h2s:'S',co:'[C-]#[O+]',n2:'N#N',o2:'O=O',methanol:'CO',cyclohexane:'C1CCCCC1',urea:'NC(N)=O',acetone:'CC(C)=O',pyridine:'c1ccncc1',
 imidazole:'c1c[nH]cn1',serine:'NC(CO)C(=O)O',cysteine:'NC(CS)C(=O)O',valine:'NC(C(C)C)C(=O)O',formaldehyde:'C=O'};
const NOH=LIB.filter(e=>e.noH).map(e=>e.id),PEPTIDES=LIB.filter(e=>e.cat==='Idealised peptide models');
test('all 43 built-in molecules are verified: SMILES oracle, peptide check, or (B-DNA only) declared schematic with no real atoms to protonate',()=>{
 assert.equal(LIB.length,43);for(const e of LIB)assert.ok(SMILES[e.id]||PEPTIDES.includes(e)||e.noH,e.id+': not verified');
 assert.deepEqual(NOH,['bdna']);assert.equal(PEPTIDES.length,11)});
for(const e of LIB.filter(e=>SMILES[e.id]))test(`${e.id}: hydrogen count, formula and connectivity match ${SMILES[e.id]}`,()=>{
 const m=e.get().atoms,bonds=computeBonds(m),ref=parseSmiles(SMILES[e.id]),els=m.map(x=>x.el);
 assert.equal(m.length,ref.atoms.length,'atom count (H included)');assert.equal(els.filter(x=>x==='H').length,ref.atoms.filter(x=>x.el==='H').length,'hydrogen count');
 assert.equal(bonds.length,ref.bonds.length,'bond count');assert.deepEqual(sig(els,bonds),sig(ref.atoms.map(x=>x.el),ref.bonds),'atom environments');
 // geometry: sensible X-H lengths and no unbonded atoms closer than 1.45 A
 const B=new Set(bonds.map(([i,j])=>i+','+j)),dist=(p,q)=>Math.hypot(m[p].x-m[q].x,m[p].y-m[q].y,m[p].z-m[q].z),RNG={C:[1.04,1.12],N:[.98,1.05],O:[.93,.99],S:[1.30,1.38]};
 for(const [i,j] of bonds){const[x,y]=m[i].el==='H'?[m[j],m[i]]:[m[i],m[j]];if(y.el==='H'&&RNG[x.el]){const r=dist(i,j);assert.ok(r>=RNG[x.el][0]&&r<=RNG[x.el][1],`${x.el}-H ${r.toFixed(3)}`)}}
 for(let i=0;i<m.length;i++)for(let j=i+1;j<m.length;j++)if(!B.has(i+','+j))assert.ok(dist(i,j)>1.45,`non-bonded ${m[i].el}${i} ${m[j].el}${j} only ${dist(i,j).toFixed(2)} A apart`)});
test('hydrogens sit where the valence says: tetrahedral C-H, trigonal NH2 and aromatic C-H, linear HCN',()=>{
 const ang=(m,i,j,k)=>Core.angle(m[i],m[j],m[k]),get=id=>LIB.find(e=>e.id===id).get().atoms;
 const ph=get('phenol');for(const i of[1,2,3,4,5]){const hh=ph.findIndex((x,k)=>x.el==='H'&&Math.hypot(x.x-ph[i].x,x.y-ph[i].y)<1.2);assert.ok(Math.abs(ph[hh].z)<1e-9,'ring H in the ring plane')}
 const an=get('aniline'),nIdx=an.findIndex(x=>x.el==='N'),hs=an.map((x,i)=>x.el==='H'&&Math.hypot(x.x-an[nIdx].x,x.y-an[nIdx].y,x.z-an[nIdx].z)<1.1?i:-1).filter(i=>i>=0);
 assert.equal(hs.length,2);assert.ok(Math.abs(ang(an,hs[0],nIdx,hs[1])-120)<1,'planar NH2');
 const hc=get('hcn');assert.ok(Math.abs(ang(hc,0,1,2)-180)<1e-6);
 const me=get('methane');assert.ok(Math.abs(ang(me,1,0,2)-109.47)<.1)});
test('L-alanine and L-serine keep the L configuration after hydrogens were added (CORN: CO, R, N clockwise from H)',()=>{
 for(const id of['alanine','serine','cysteine','valine']){const m=LIB.find(e=>e.id===id).get().atoms,v=(i,j)=>[m[i].x-m[j].x,m[i].y-m[j].y,m[i].z-m[j].z],
  CA=1,C=2,N=0,R=5,hCA=m.findIndex((x,k)=>x.el==='H'&&Math.hypot(...v(k,CA))<1.15&&k>R),u=v(hCA,CA),co=v(C,CA),r=v(R,CA),n=v(N,CA);
  // looking from H towards CA, CO -> R -> N should run clockwise for L: triple product sign checked against the heavy-atom scaffold used before (>0 L)
  const tp=co[0]*(r[1]*n[2]-r[2]*n[1])-co[1]*(r[0]*n[2]-r[2]*n[0])+co[2]*(r[0]*n[1]-r[1]*n[0]);assert.ok(tp>0,id+' triple product '+tp.toFixed(2));assert.ok(hCA>0&&Math.hypot(...u)>1)}});
test('the inlined Hydrogens block equals src/core/hydrogens.mjs',()=>{const src=fs.readFileSync(new URL('../src/core/hydrogens.mjs',import.meta.url),'utf8').replace(/\r\n/g,'\n').replace(/^export /gm,''),i=html.indexOf('/*HYDROGENS-BEGIN*/'),j=html.indexOf('/*HYDROGENS-END*/');assert.ok(i>0&&j>i&&html.slice(i,j).includes(src))});
import {layout2D,drawnOrdersAndCharges} from '../src/core/layout2d.mjs';
test('2D structure diagram of every all-atom built-in molecule: every H drawn as a bonded atom, no overlaps; hidden H still counted for labels',()=>{
 let n=0;for(const e of LIB.filter(e=>!e.noH&&!PEPTIDES.includes(e))){const m=e.get().atoms.map(x=>({...x,element:x.el,resName:x.res,hetero:false})),bonds=computeBonds(m),ord=Core.assignBondOrders(m,bonds),dr=drawnOrdersAndCharges(m,bonds,ord),ids=m.map((_,i)=>i);
  const hasHeavy=m.some(x=>x.el!=='H'),r=layout2D(m,bonds,dr.orders,ids,{explicitH:true});if(!hasHeavy){assert.equal(r.ok,false);continue}
  assert.ok(r.ok,e.id+': '+r.reason);assert.equal(r.pos.size,m.length,e.id+': every atom positioned');assert.equal(r.overlaps,0,e.id+' overlaps');
  const hid=layout2D(m,bonds,dr.orders,ids.filter(i=>m[i].el!=='H'),{explicitH:true});assert.equal(hid.pos.size,m.filter(x=>x.el!=='H').length);
  for(const[i,c]of hid.hcount)assert.equal(c,bonds.filter(([a,b])=>(a===i&&m[b].el==='H')||(b===i&&m[a].el==='H')).length,e.id+' hidden H count of atom '+i);n++}
 assert.ok(n>=25)});
test('functional-group labels when hydrogens are hidden: phenol OH, aniline NH2, cysteine SH, glycine NH2 + C=O + OH, zwitterion NH3+ and carboxylate O-',()=>{
 // each heteroatom written as symbol + implicit H count + charge, e.g. N2 = NH2, N3+ = NH3+, O-1 = O-
 const lab=id=>{const m=LIB.find(e=>e.id===id).get().atoms.map(x=>({...x,element:x.el,resName:x.res,hetero:false})),bonds=computeBonds(m),dr=drawnOrdersAndCharges(m,bonds,Core.assignBondOrders(m,bonds)),r=layout2D(m,bonds,dr.orders,m.map((_,i)=>i).filter(i=>m[i].el!=='H'),{explicitH:true});
  return[...r.pos.keys()].filter(i=>m[i].el!=='C').map(i=>m[i].el+(r.hcount.get(i)?'H'+r.hcount.get(i):'')+(dr.charge.get(i)?(dr.charge.get(i)>0?'+':'-'):'')).sort().join(' ')};
 const EXP={phenol:'OH1',aniline:'NH2',ethanol:'OH1',water:'OH2',ammonia:'NH3',cysteine:'NH2 O OH1 SH1',glycine:'NH2 O OH1',alanine_h:'NH3+ O O-',urea:'NH2 NH2 O',imidazole:'N NH1',acetic:'O OH1'};
 for(const[id,x]of Object.entries(EXP))assert.equal(lab(id),x,id)});

// Peptide models: expected formula from the residue list (each residue as free amino acid Gly C2H5NO2 or Ala-like C3H7NO2, minus one H2O per peptide bond),
// then the connectivity of every atom from the bond graph: degrees by atom role.
for(const e of PEPTIDES)test(`${e.id}: free-peptide hydrogens: formula, per-atom connectivity, bond lengths, no heavy-atom clashes`,()=>{
 const m=e.get().atoms,els=m.map(a=>({...a,element:a.el})),bonds=computeBonds(els),deg=m.map(()=>[]);bonds.forEach(([i,j])=>{deg[i].push(j);deg[j].push(i)});
 const res=new Map();m.forEach((a,i)=>{const k=a.chain+'|'+a.resi;if(!res.has(k))res.set(k,{name:a.res,chain:a.chain,at:{}});res.get(k).at[a.name]=i});
 const chains=new Map();for(const r of res.values())chains.set(r.chain,(chains.get(r.chain)||0)+1);
 let C=0,H=0,N=0,O=0;for(const r of res.values()){const g=r.name==='GLY';C+=g?2:3;H+=g?5:7;N+=1;O+=2}for(const n of chains.values()){H-=2*(n-1);O-=n-1}
 const cnt=x=>m.filter(a=>a.el===x).length;assert.deepEqual([cnt('C'),cnt('H'),cnt('N'),cnt('O')],[C,H,N,O],'formula');assert.equal(m.length,C+H+N+O,'no other elements');
 const first=new Set(),last=new Map();for(const[c]of chains){const rs=[...res.values()].filter(r=>r.chain===c);first.add(rs[0]);last.set(rs.at(-1),1)}
 for(const r of res.values()){const A=r.at,f=first.has(r),l=last.has(r),g=r.name==='GLY',nb=i=>deg[i].map(j=>m[j].el).sort().join('');
  assert.equal(nb(A.N),f?'CHH':'CCH','N of '+r.name+r.chain);assert.equal(nb(A.CA),g?'CHHN':'CCHN','CA');
  assert.equal(nb(A.O),'C');assert.equal(nb(A.C),l?'COO':'CNO','C');if(!g)assert.equal(nb(A.CB),'CHHH','CB');if(l)assert.equal(nb(A.OXT),'CH','OXT')}
 const B=new Set(bonds.map(([i,j])=>i+','+j)),dist=(p,q)=>Math.hypot(m[p].x-m[q].x,m[p].y-m[q].y,m[p].z-m[q].z),RNG={C:[1.04,1.12],N:[.98,1.05],O:[.93,.99]};
 for(const[i,j]of bonds){if(m[i].el==='H'||m[j].el==='H'){const x=m[i].el==='H'?m[j]:m[i];const r=dist(i,j);assert.ok(r>=RNG[x.el][0]&&r<=RNG[x.el][1],x.el+'-H '+r.toFixed(3));assert.equal(deg[m[i].el==='H'?i:j].length,1,'every H has exactly one parent')}}
 for(let i=0;i<m.length;i++)for(let j=i+1;j<m.length;j++)if(!B.has(i+','+j)){const hh=m[i].el==='H'&&m[j].el==='H';assert.ok(dist(i,j)>(hh?1.25:1.7),`non-bonded ${m[i].name}${m[i].resi}-${m[j].name}${m[j].resi} ${dist(i,j).toFixed(2)}`)}});
