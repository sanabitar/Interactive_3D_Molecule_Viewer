import test from 'node:test';import assert from 'node:assert/strict';
import {lonePairCount,computeLonePairs} from '../src/core/lonepairs.mjs';import {drawnOrdersAndCharges} from '../src/core/layout2d.mjs';
import {env,viewerAtoms} from './helpers/library-env.mjs';
const {LIB,Core,computeBonds}=env;
// Everything is analysed exactly as the viewer does it: library entry -> atom records -> bonds from geometry -> perceived orders -> orders/charges as drawn -> lone pairs.
function analyse(id,{dropH=false}={}){let atoms=viewerAtoms(LIB.find(e=>e.id===id));if(dropH)atoms=atoms.filter(a=>a.el!=='H');
 const bonds=computeBonds(atoms),od=Core.assignBondOrdersDetailed(atoms,bonds),dr=drawnOrdersAndCharges(atoms,bonds,od.orders),lp=computeLonePairs(atoms,bonds,dr.orders,{charges:dr.charge});return{atoms,bonds,lp,dr}}
const sub=(a,b)=>[a[0]-b[0],a[1]-b[1],a[2]-b[2]],dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2],len=a=>Math.hypot(...a),unit=a=>{const l=len(a);return a.map(x=>x/l)},cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const ang=(a,b)=>Math.acos(Math.max(-1,Math.min(1,dot(unit(a),unit(b)))))*180/Math.PI,pos=a=>[a.x,a.y,a.z];
const nbrs=(R,i)=>R.bonds.filter(([a,b])=>a===i||b===i).map(([a,b])=>a===i?b:a),toward=(R,i,j)=>sub(pos(R.atoms[j]),pos(R.atoms[i]));
const lps=(R,i)=>R.lp.filter(l=>l.atom===i);

test('lone-pair count from valence: water O 2, ammonia N 1, amine N 1, halides 3, ions with their charge; incomplete bonding gives none',()=>{
 const T=[['O',2,0,2],['O',1,-1,3],['O',3,1,1],/* H3O+ */['N',3,0,1],['N',4,1,0],['N',4,0,0],['N',3.5,0,0],['S',2,0,2],['S',6,0,0],['Cl',1,0,3],['Cl',0,-1,4],['F',1,0,3],['P',3,0,1],['P',5,0,0],['C',0,0,0],['H',1,0,0],['Zn',0,0,0],
  ['O',1,0,0],['N',1,0,0],['N',2,0,0],['O',0,0,0]];   // the last four: bonding looks incomplete (hydrogens missing), so nothing is guessed
 for(const[el,sum,q,want]of T)assert.equal(lonePairCount(el,sum,q),want,`${el} bond-order sum ${sum} charge ${q}`)});

// Expected lone pairs per heavy atom, written from chemistry (Lewis structures), NOT derived from the code under test.
const EXPECT={water:{O:[2]},methane:{},benzene:{},co2:{O:[2,2]},ammonia:{N:[1]},ethene:{},ethyne:{},ethane:{},glycine:{N:[1],O:[2,2]},alanine:{N:[1],O:[2,2]},
 alanine_h:{O:[2,3]},/* NH3+ has none; C=O 2; O- 3 */ethanol:{O:[2]},acetic:{O:[2,2]},phenol:{O:[2]},toluene:{},aniline:{N:[1]},hcn:{N:[1]},h2s:{S:[2]},co:{C:[1],O:[1]},n2:{N:[1,1]},o2:{O:[2,2]},
 methanol:{O:[2]},cyclohexane:{},urea:{N:[1,1],O:[2]},acetone:{O:[2]},pyridine:{N:[1]},imidazole:{N:[1]}/* N-H lone pair is in the aromatic pi system: not shown */,serine:{N:[1],O:[2,2,2]},cysteine:{N:[1],O:[2,2],S:[2]},valine:{N:[1],O:[2,2]},formaldehyde:{O:[2]}};
for(const[id,want]of Object.entries(EXPECT))test(`${id}: lone pairs per atom match the Lewis structure`,()=>{
 const R=analyse(id),per=new Map();for(const l of R.lp)per.set(l.atom,(per.get(l.atom)||0)+1);const got={};for(const[i,n]of per)(got[R.atoms[i].el]||(got[R.atoms[i].el]=[])).push(n);for(const k in got)got[k].sort();
 assert.deepEqual(got,Object.fromEntries(Object.entries(want).map(([k,v])=>[k,[...v].sort()])))});
test('peptide models: every amide / N-terminal N has 1, every carbonyl O 2, the C-terminal acid O and OH 2 each, nothing else: 3n+2 per chain',()=>{
 for(const e of LIB.filter(x=>x.cat==='Idealised peptide models')){const R=analyse(e.id),per=new Map();for(const l of R.lp)per.set(l.atom,(per.get(l.atom)||0)+1);
  const chains=new Set(R.atoms.map(a=>a.chain)),n=new Set(R.atoms.map(a=>a.chain+'|'+a.resi)).size;assert.equal(R.lp.length,3*n+2*chains.size,e.id);
  for(const[i,c]of per){const a=R.atoms[i];assert.ok(a.el==='N'?c===1:a.el==='O'?c===2:false,`${e.id} ${a.name}${a.resi} has ${c}`)}
  assert.equal([...per.keys()].filter(i=>R.atoms[i].el==='N').length,n,e.id+': one N per residue')}});
test('B-DNA is schematic (P and X pseudo-atoms): no real atoms, so nothing is derived for it by the viewer',()=>{assert.ok(LIB.find(e=>e.id==='bdna').noH)});

test('water: two lone pairs, tetrahedral with the two O-H bonds (LP-O-H ~109 deg, LP-O-LP ~109 deg), opposite the H atoms',()=>{
 const R=analyse('water'),L=lps(R,0),[h1,h2]=nbrs(R,0);assert.equal(L.length,2);assert.ok(L.every(l=>l.type==='sp3'));
 for(const l of L)for(const h of[h1,h2]){const a=ang(l.dir,toward(R,0,h));assert.ok(a>105&&a<118,'LP-O-H '+a)}
 const aa=ang(L[0].dir,L[1].dir);assert.ok(aa>105&&aa<115,'LP-O-LP '+aa);
 assert.ok(dot(unit([L[0].dir[0]+L[1].dir[0],L[0].dir[1]+L[1].dir[1],L[0].dir[2]+L[1].dir[2]]),unit([...toward(R,0,h1)].map((x,k)=>x+toward(R,0,h2)[k])))<-.99,'the pair points away from the hydrogens')});
test('ammonia: one lone pair on the pyramid axis, opposite the three N-H bonds (LP-N-H ~110 deg)',()=>{
 const R=analyse('ammonia'),L=lps(R,0);assert.equal(L.length,1);const s=nbrs(R,0).map(j=>unit(toward(R,0,j))).reduce((p,q)=>[p[0]+q[0],p[1]+q[1],p[2]+q[2]]);
 assert.ok(dot(L[0].dir,unit(s))<-.999,'opposite the sum of the bonds');for(const j of nbrs(R,0)){const a=ang(L[0].dir,toward(R,0,j));assert.ok(a>105&&a<118,'LP-N-H '+a)}});
test('amine N (glycine) and alcohol O (ethanol, methanol) are tetrahedral: every LP-X-substituent angle ~109 deg',()=>{
 for(const[id,i,n]of[['glycine',0,1],['ethanol',2,2],['methanol',1,2]]){const R=analyse(id),L=lps(R,i);assert.equal(L.length,n,id);for(const l of L){assert.equal(l.type,'sp3');for(const j of nbrs(R,i)){const a=ang(l.dir,toward(R,i,j));assert.ok(a>100&&a<122,`${id} LP-X-subst ${a}`)}}}});
test('carbonyl O (formaldehyde): two lone pairs in the molecular plane, 120 deg from C=O and from each other',()=>{
 const R=analyse('formaldehyde'),O=R.atoms.findIndex(a=>a.el==='O'),C=nbrs(R,O)[0],H=nbrs(R,C).find(j=>R.atoms[j].el==='H'),n=unit(cross(toward(R,C,O),toward(R,C,H))),L=lps(R,O);
 assert.equal(L.length,2);for(const l of L){assert.equal(l.type,'sp2');assert.ok(Math.abs(dot(l.dir,n))<.02,'in the plane');assert.ok(Math.abs(ang(l.dir,toward(R,O,C))-120)<1.5)}assert.ok(Math.abs(ang(L[0].dir,L[1].dir)-120)<2)});
test('linear atoms: N2, HCN nitrogen and CO point straight away from the bond (180 deg)',()=>{
 for(const[id,i]of[['n2',0],['n2',1],['hcn',2],['co',0],['co',1]]){const R=analyse(id),L=lps(R,i);assert.equal(L.length,1,id+i);assert.equal(L[0].type,'sp');assert.ok(Math.abs(ang(L[0].dir,toward(R,i,nbrs(R,i)[0]))-180)<.5)}});
test('lone pair conjugated with a ring (aniline N): ONE lone pair perpendicular to the ring plane (a p lone pair), not tetrahedral',()=>{
 const R=analyse('aniline'),N=R.atoms.findIndex(a=>a.el==='N'),L=lps(R,N);assert.equal(L.length,1);assert.equal(L[0].type,'p');assert.ok(Math.abs(L[0].dir[2])>.99,'perpendicular to the ring plane (z=0): '+L[0].dir)});
test('pyridine N: the lone pair lies in the ring plane pointing straight out of the ring (sp2); phenol O: one in the plane, one perpendicular',()=>{
 let R=analyse('pyridine'),L=lps(R,0);assert.equal(L.length,1);assert.equal(L[0].type,'sp2');assert.ok(Math.abs(L[0].dir[2])<1e-6);assert.ok(dot(L[0].dir,unit(pos(R.atoms[0])))>.99,'radially outward');
 R=analyse('phenol');const O=R.atoms.findIndex(a=>a.el==='O');L=lps(R,O);assert.deepEqual(L.map(l=>l.type).sort(),['p','sp2']);for(const l of L)assert.ok(l.type==='p'?Math.abs(l.dir[2])>.99:Math.abs(l.dir[2])<1e-6)});
test('carboxylate (alanine zwitterion): O- shows three (two in the plane + one p), C=O two; NH3+ shows none; neutral COOH shows two on each O',()=>{
 let R=analyse('alanine_h');const os=R.atoms.map((a,i)=>a.el==='O'?i:-1).filter(i=>i>=0).map(i=>lps(R,i).map(l=>l.type).sort().join('+'));assert.deepEqual(os.sort(),['p+sp2+sp2','sp2+sp2']);
 assert.equal(lps(R,R.atoms.findIndex(a=>a.el==='N')).length,0,'ammonium N');
 R=analyse('acetic');for(const i of R.atoms.map((a,k)=>a.el==='O'?k:-1).filter(k=>k>=0))assert.equal(lps(R,i).length,2)});
test('lone pairs are not atoms: the result holds only {atom,dir,type}, directions are unit vectors, the atom array is untouched',()=>{
 const atoms=viewerAtoms(LIB.find(e=>e.id==='water')),before=JSON.stringify(atoms),bonds=computeBonds(atoms),orders=Core.assignBondOrdersDetailed(atoms,bonds).orders,lp=computeLonePairs(atoms,bonds,orders);
 assert.equal(JSON.stringify(atoms),before);assert.equal(atoms.length,3);assert.equal(bonds.length,2);for(const l of lp){assert.deepEqual(Object.keys(l).sort(),['atom','dir','type']);assert.ok(Math.abs(len(l.dir)-1)<1e-9)}});
test('hydrogens missing from the file: atoms that look short of bonds get NO lone pair (not guessed); complete atoms still do; delocalised atoms are skipped',()=>{
 const R=analyse('glycine',{dropH:true}),O=R.atoms.map((a,i)=>a.el==='O'?i:-1).filter(i=>i>=0),N=R.atoms.findIndex(a=>a.el==='N');
 assert.equal(lps(R,N).length,0,'N lost its H: unknown');assert.deepEqual(O.map(i=>lps(R,i).length).sort(),[0,2],'C=O keeps 2; the O-H without H gets none');
 assert.equal(analyse('pyridine',{dropH:true}).lp.length,0,'ring N could be pyridine-type or N-H without hydrogens');assert.equal(analyse('pyridine').lp.length,1);
 assert.equal(analyse('water',{dropH:true}).lp.length,0,'a bare O (water without H) is not guessed')});
test('structures from files (no hydrogens at all, residue orders): backbone carbonyl O still gets 2 lone pairs, backbone N and side-chain OH are left alone',()=>{
 const R=analyse('pep_alpha',{dropH:true}),perEl={N:0,O:0};for(const l of R.lp)perEl[R.atoms[l.atom].el]++;assert.equal(perEl.N,0,'N-H nitrogens: hydrogens unknown');assert.equal(perEl.O,11*2+2+3,'11 backbone carbonyls x 2, plus the C-terminus read as a carboxylate (C=O 2, O- 3)')});
test('drawnOrdersAndCharges: a delocalised group with an O-H is a protonated acid (C=O and C-OH); two terminal O is still a carboxylate (C=O and C-O-)',()=>{
 const at=['C','O','O','H'].map(el=>({el})),acid=drawnOrdersAndCharges(at,[[0,1],[0,2],[2,3]],[1.5,1.5,1]);assert.deepEqual(acid.orders.slice(0,2),[2,1]);assert.equal(acid.charge.size,0);
 const ox=drawnOrdersAndCharges(['C','O','O'].map(el=>({el})),[[0,1],[0,2]],[1.5,1.5]);assert.deepEqual(ox.orders,[2,1]);assert.equal(ox.charge.get(2),-1)});
