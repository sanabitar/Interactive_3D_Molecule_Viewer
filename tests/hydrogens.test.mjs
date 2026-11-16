import test from 'node:test';import assert from 'node:assert/strict';
import {missingHydrogens,addHydrogens,xhLength} from '../src/core/hydrogens.mjs';
import {layout2D,atomLabelParts,drawnOrdersAndCharges} from '../src/core/layout2d.mjs';
import {zmatrix} from '../src/core/build.mjs';import {distance,angle} from '../src/core/geometry.mjs';
test('missing hydrogens follow valence: C4 N3 O2 S2, charged N+ and O-, nothing for other elements',()=>{
 assert.equal(missingHydrogens('C',4),0);assert.equal(missingHydrogens('C',1),3);assert.equal(missingHydrogens('N',1),2);assert.equal(missingHydrogens('O',1),1);assert.equal(missingHydrogens('S',1),1);
 assert.equal(missingHydrogens('N',1,1),3,'NH3+');assert.equal(missingHydrogens('O',1,-1),0,'alkoxide / carboxylate O-');assert.equal(missingHydrogens('P',3),0);assert.equal(missingHydrogens('C',5),0,'never negative')});
const A=(els)=>els.map(el=>({el,x:0,y:0,z:0}));
test('methanol: CH3 + OH, standard lengths, tetrahedral, staggered O-H',()=>{
 const m=zmatrix([['C'],['O',0,1.43]]),r=addHydrogens(m,[[0,1]],[1]);assert.equal(r.hydrogens,4);assert.equal(r.atoms.length,6);
 for(const[p,h]of r.bonds.slice(1))assert.ok(Math.abs(distance(r.atoms[p],r.atoms[h])-xhLength(r.atoms[p].el))<1e-9);
 const ch=r.bonds.filter(b=>b[0]===0&&b[1]>1).map(b=>b[1]);for(const h of ch)assert.ok(Math.abs(angle(r.atoms[h],r.atoms[0],r.atoms[1])-109.5)<1e-6);
 assert.ok(Math.abs(angle(r.atoms[ch[0]],r.atoms[0],r.atoms[ch[1]])-109.47)<.05)});
test('sp3 / sp2 / sp decided by Kekule bond orders: ethene =CH2 planar 120, ethyne linear, C=O needs none',()=>{
 let r=addHydrogens(zmatrix([['C'],['C',0,1.33]]),[[0,1]],[2]);assert.equal(r.hydrogens,4);const h0=r.bonds.filter(b=>b[0]===0&&b[1]>1).map(b=>b[1]);
 assert.ok(Math.abs(angle(r.atoms[h0[0]],r.atoms[0],r.atoms[1])-120)<1e-6&&Math.abs(angle(r.atoms[h0[0]],r.atoms[0],r.atoms[h0[1]])-120)<1e-6);
 r=addHydrogens(zmatrix([['C'],['C',0,1.2]]),[[0,1]],[3]);assert.equal(r.hydrogens,2);for(const[p,h]of r.bonds.slice(1))assert.ok(Math.abs(angle(r.atoms[h],r.atoms[p],r.atoms[p?0:1])-180)<1e-6);
 r=addHydrogens(zmatrix([['C'],['O',0,1.21]]),[[0,1]],[2]);assert.equal(r.hydrogens,2,'formaldehyde: CH2=O');
 r=addHydrogens(zmatrix([['O'],['C',0,1.16],['O',1,1.16,0,180]]),[[0,1],[1,2]],[2,2]);assert.equal(r.hydrogens,0,'CO2')});
test('charges and explicit torsions are honoured',()=>{
 const r=addHydrogens(zmatrix([['N'],['C',0,1.49]]),[[0,1]],[1],{charges:{0:1}});assert.equal(r.atoms.filter(a=>a.parent===0).length,3,'NH3+');
 const a=addHydrogens(zmatrix([['C'],['C',0,1.52],['O',1,1.43,0,108]]),[[0,1],[1,2]],[1,1],{torsions:{2:180}}),b=addHydrogens(zmatrix([['C'],['C',0,1.52],['O',1,1.43,0,108]]),[[0,1],[1,2]],[1,1],{torsions:{2:60}});
 const hO=r=>r.atoms.findIndex(x=>x.parent===2);assert.ok(distance(a.atoms[hO(a)],a.atoms[0])>distance(b.atoms[hO(b)],b.atoms[0])+.5,'anti H is farther from the methyl than gauche H')});
const G=(els,edges,ids)=>({atoms:els.map(el=>({el})),bonds:edges.map(e=>[e[0],e[1]]),orders:edges.map(e=>e[2]||1),ids:ids||els.map((_,i)=>i)});
test('explicit H: opt-in, bonded atoms placed 0.8 away in free space, no overlaps; default layout is unchanged',()=>{
 const m=G(['C','C','O','H','H','H','H','H','H'],[[0,1],[1,2],[0,3],[0,4],[0,5],[1,6],[1,7],[2,8]]);
 let r=layout2D(m.atoms,m.bonds,m.orders,m.ids);assert.equal(r.pos.size,3,'default: hydrogens only counted');assert.equal(r.hcount.get(0),3);
 r=layout2D(m.atoms,m.bonds,m.orders,m.ids,{explicitH:true});assert.equal(r.pos.size,9);assert.equal(r.overlaps,0);assert.equal(r.hdrawn.get(0),3);
 for(const[h,p]of[[3,0],[4,0],[5,0],[6,1],[7,1],[8,2]])assert.ok(Math.abs(Math.hypot(r.pos.get(h)[0]-r.pos.get(p)[0],r.pos.get(h)[1]-r.pos.get(p)[1])-.8)<1e-9);
 const ids=[0,1,2,3];r=layout2D(m.atoms,m.bonds,m.orders,ids,{explicitH:true});assert.equal(r.pos.size,4,'only listed hydrogens are drawn');assert.equal(r.hcount.get(0),3,'the rest are still counted');assert.equal(r.hdrawn.get(0)||0,1)});
test('explicit H on lone atoms: water is bent, ammonia a tripod, methane a cross; HCN hydrogen opposite the nitrile',()=>{
 let r=layout2D(...Object.values(G(['O','H','H'],[[0,1],[0,2]])).slice(0,4),{explicitH:true}),a=r.pos.get(1),b=r.pos.get(2),o=r.pos.get(0);
 const deg=Math.acos(((a[0]-o[0])*(b[0]-o[0])+(a[1]-o[1])*(b[1]-o[1]))/.64)*180/Math.PI;assert.ok(Math.abs(deg-104.5)<.1,'H-O-H '+deg);
 const m=G(['N','H','H','H'],[[0,1],[0,2],[0,3]]);r=layout2D(m.atoms,m.bonds,m.orders,m.ids,{explicitH:true});assert.equal(r.overlaps,0);
 const c=G(['C','H','H','H','H'],[[0,1],[0,2],[0,3],[0,4]]);r=layout2D(c.atoms,c.bonds,c.orders,c.ids,{explicitH:true});assert.equal(r.overlaps,0);
 const n=G(['C','N','H'],[[0,1,3],[0,2]]);r=layout2D(n.atoms,n.bonds,n.orders,n.ids,{explicitH:true});const p=[0,1,2].map(i=>r.pos.get(i));assert.ok(Math.abs((p[1][0]-p[0][0])*(p[2][1]-p[0][1])-(p[1][1]-p[0][1])*(p[2][0]-p[0][0]))<1e-9,'H-C#N is straight'); assert.ok((p[1][0]-p[0][0])*(p[2][0]-p[0][0])+(p[1][1]-p[0][1])*(p[2][1]-p[0][1])<0,'H on the opposite side of C from N')});
test('ring and chain hydrogens do not collide: cyclohexane, toluene and a pentane chain drawn with every H',()=>{
 const ring=(n)=>Array.from({length:n},(_,i)=>[i,(i+1)%n]);
 for(const[name,els0,edges0]of[['cyclohexane',Array(6).fill('C'),ring(6)],['toluene',Array(7).fill('C'),[...ring(6),[0,6]]],['pentane',Array(5).fill('C'),[[0,1],[1,2],[2,3],[3,4]]]]){
  const els=[...els0],edges=edges0.map(e=>[...e]);const heavy=els.length,deg=els.map(()=>0);edges.forEach(([a,b])=>{deg[a]++;deg[b]++});
  const arom=new Set(name==='toluene'?[0,1,2,3,4,5]:[]);for(let i=0;i<heavy;i++)for(let k=0;k<(arom.has(i)?3:4)-deg[i];k++){edges.push([i,els.length]);els.push('H')}
  const m=G(els,edges),r=layout2D(m.atoms,m.bonds,m.orders,m.ids,{explicitH:true});assert.ok(r.ok&&r.overlaps===0,name+' overlaps '+r.overlaps)}});
test('labels: OH / NH2 / SH after the symbol, HO / H2N before it when bonds leave to the right, charges as superscripts',()=>{
 assert.deepEqual(atomLabelParts('O',1),{pre:'',sym:'O',post:'H'});assert.deepEqual(atomLabelParts('N',2),{pre:'',sym:'N',post:'H₂'});assert.deepEqual(atomLabelParts('S',1,0,true),{pre:'H',sym:'S',post:''});
 assert.deepEqual(atomLabelParts('N',2,0,true),{pre:'H₂',sym:'N',post:''});assert.deepEqual(atomLabelParts('N',3,1),{pre:'',sym:'N',post:'H₃⁺'});assert.deepEqual(atomLabelParts('O',0,-1),{pre:'',sym:'O',post:'⁻'});
 assert.deepEqual(atomLabelParts('C',3),{pre:'',sym:'C',post:'H₃'});assert.deepEqual(atomLabelParts('O',0),{pre:'',sym:'O',post:''});assert.deepEqual(atomLabelParts('Fe',0,2),{pre:'',sym:'Fe',post:'²⁺'})});
test('drawn orders/charges: carboxylate shown as C=O plus C-O-, 4-coordinate N is N+, nothing else is guessed',()=>{
 const at=['C','O','O','N','H','H','H','C'].map(el=>({el})),bonds=[[0,1],[0,2],[7,3],[3,4],[3,5],[3,6],[0,7]],r=drawnOrdersAndCharges(at,bonds,[1.5,1.5,1,1,1,1,1]);
 assert.deepEqual(r.orders.slice(0,2),[2,1]);assert.equal(r.charge.get(2),-1);assert.equal(r.charge.get(3),1);assert.equal(r.charge.size,2);
 const lone=drawnOrdersAndCharges([{el:'C'},{el:'O'}],[[0,1]],[1]);assert.equal(lone.charge.size,0,'a bare O may be an OH whose H is missing')});
import {validateSession} from '../src/core/session.mjs';
test('session settings: hydrogens default to hidden and element labels to shown; old sessions stay readable; bad values are rejected',()=>{
 const base={laura:'0.3',structure:{id:'water',name:'Water',builtin:true}},V=o=>validateSession(JSON.stringify(o));
 let r=V(base);assert.equal(r.ok,true);assert.equal(r.session.settings.hydrogens,false);assert.equal(r.session.settings.elements,true);
 r=V({...base,settings:{hydrogens:true,elements:false}});assert.equal(r.session.settings.hydrogens,true);assert.equal(r.session.settings.elements,false);
 r=V({...base,settings:{elements:'yes'}});assert.equal(r.ok,false);assert.match(r.error,/elements/)});
