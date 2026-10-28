import test from 'node:test';import assert from 'node:assert/strict';
import {assignBondOrders,hybridization,geometryCheck} from '../src/core/bonds.mjs';import {zmatrix,peptide} from '../src/core/build.mjs';
const D=Math.PI/180,at=(el,x,y,z)=>({element:el,x,y,z});
const ord=(A,B)=>assignBondOrders(A,B),ang=(r,deg)=>[r*Math.cos(deg*D),r*Math.sin(deg*D)];
test('bond orders from textbook lengths',()=>{let A=[at('C',0,0,0),at('C',1.54,0,0)];assert.deepEqual(ord(A,[[0,1]]),[1],'ethane C-C');A[1].x=1.33;assert.deepEqual(ord(A,[[0,1]]),[2],'C=C');A[1].x=1.20;assert.deepEqual(ord(A,[[0,1]]),[3],'C#C');A[1].x=1.397;assert.deepEqual(ord(A,[[0,1]]),[1.5],'aromatic');A[1].x=1.46;assert.deepEqual(ord(A,[[0,1]]),[1],'conjugated single');
 const o=(e,d)=>ord([at('C',0,0,0),at(e,d,0,0)],[[0,1]])[0];assert.equal(o('O',1.21),2);assert.equal(o('O',1.128),3,'CO');assert.equal(o('O',1.43),1);assert.equal(o('N',1.156),3,'HCN');assert.equal(o('N',1.329),1,'peptide C-N drawn single');
 assert.equal(ord([at('N',0,0,0),at('N',1.098,0,0)],[[0,1]])[0],3);assert.equal(ord([at('O',0,0,0),at('O',1.208,0,0)],[[0,1]])[0],2);assert.equal(ord([at('C',0,0,0),at('H',1.09,0,0)],[[0,1]])[0],1)});
test('CO2 keeps two double bonds; carboxylate becomes delocalised; carboxylic acid keeps one double',()=>{let A=[at('O',-1.16,0,0),at('C',0,0,0),at('O',1.16,0,0)];assert.deepEqual(ord(A,[[0,1],[1,2]]),[2,2]);
 const Z=(d2)=>zmatrix([['C'],['C',0,1.52],['O',1,1.25,0,117.5],['O',1,d2,0,117.5,2,180]]);let B=Z(1.25);assert.deepEqual(ord(B,[[0,1],[1,2],[1,3]]),[1,1.5,1.5],'carboxylate');B=Z(1.36);assert.deepEqual(ord(B,[[0,1],[1,2],[1,3]]),[1,2,1],'acid')});
test('hybridization rules',()=>{const H=(A,B)=>hybridization(A,B,ord(A,B));
 let A=[at('C',0,0,0),at('C',1.54,0,0)];assert.deepEqual(H(A,[[0,1]]),['sp3','sp3']);A[1].x=1.33;assert.deepEqual(H(A,[[0,1]]),['sp2','sp2']);A[1].x=1.2;assert.deepEqual(H(A,[[0,1]]),['sp','sp']);
 A=[at('O',-1.16,0,0),at('C',0,0,0),at('O',1.16,0,0)];assert.deepEqual(H(A,[[0,1],[1,2]]),['sp2','sp','sp2'],'CO2');
 assert.deepEqual(H([at('C',0,0,0),at('O',1.128,0,0)],[[0,1]]),['sp','sp'],'carbon monoxide: triple-bonded O is sp');
 A=[at('N',0,0,0),at('H',1.01,0,0),at('H',0,1.01,0),at('H',0,0,1.01)];assert.equal(H(A,[[0,1],[0,2],[0,3]])[0],'sp3','ammonia N');
 const P=peptide(4,{phi:-57,psi:-47}).map(a=>({...a,element:a.el})),bonds=[];for(let i=0;i<P.length;i++)for(let j=i+1;j<P.length;j++){const d=Math.hypot(P[i].x-P[j].x,P[i].y-P[j].y,P[i].z-P[j].z);if(d<1.8)bonds.push([i,j])}
 const h=hybridization(P,bonds,ord(P,bonds));assert.equal(h[P.findIndex(a=>a.resi===2&&a.name==='N')],'sp2','amide N');assert.equal(h[P.findIndex(a=>a.resi===2&&a.name==='C')],'sp2','carbonyl C');assert.equal(h[P.findIndex(a=>a.resi===2&&a.name==='CA')],'sp3','C-alpha');assert.equal(h[P.findIndex(a=>a.resi===2&&a.name==='O')],'sp2','carbonyl O')});
const geo=(A,B)=>geometryCheck(A,B,hybridization(A,B,ord(A,B)));
test('geometry check accepts ideal geometry and flags inconsistent geometry',()=>{const t=1.09/Math.sqrt(3),S=[[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]];
 let M=[at('C',0,0,0),...S.map(([x,y,z])=>at('H',t*x,t*y,t*z))],B=[[0,1],[0,2],[0,3],[0,4]];let r=geo(M,B);assert.ok(r.maxDeviation<0.05&&r.outliers.length===0,'tetrahedral methane (109.47 vs 109.5 convention)');
 M=[at('C',0,0,0),at('H',1.09,0,0),at('H',-1.09,0,0),at('H',0,1.09,0),at('H',0,-1.09,0)];r=geo(M,B);assert.ok(r.outliers.length>0&&r.maxDeviation>=70,'square-planar "methane" is sp3 with 90/180 deg angles');
 const e=(deg)=>{const[x,y]=ang(1.06,deg);return[at('C',0,0,0),at('C',1.2,0,0),at('H',-1.06,0,0),at('H',1.2+x,y,0)]};r=geo(e(0),[[0,1],[0,2],[1,3]]);assert.equal(r.outliers.length,0,'linear ethyne');r=geo(e(60),[[0,1],[0,2],[1,3]]);assert.ok(r.outliers.length>0,'bent "ethyne" is sp with a 120 deg angle');
 const w=[at('O',0,0,0),at('H',.756,.586,0),at('H',-.756,.586,0)];r=geo(w,[[0,1],[0,2]]);assert.ok(r.maxDeviation>4&&r.outliers.length===0,'water 104.5 vs 109.5: within lone-pair tolerance')});
// ---- ring perception ----
const poly=(els,R,dz=0)=>els.map((e,i)=>({element:e,x:R*Math.cos(2*Math.PI*i/els.length),y:R*Math.sin(2*Math.PI*i/els.length),z:i%2?dz:-dz})),cyc=n=>Array.from({length:n},(_,i)=>[i,(i+1)%n]);
test('ring perception: pyridine and imidazole become fully aromatic; benzene stays aromatic',()=>{
 const py=poly(['N','C','C','C','C','C'],1.38);assert.deepEqual(assignBondOrders(py,cyc(6)),[1.5,1.5,1.5,1.5,1.5,1.5],'pyridine: C-N bonds no longer drawn single');
 assert.deepEqual(hybridization(py,cyc(6),assignBondOrders(py,cyc(6))),Array(6).fill('sp2'),'all ring atoms sp2');
 const im=poly(['N','C','N','C','C'],1.165);assert.ok(Math.hypot(im[0].x-im[1].x,im[0].y-im[1].y)>1.36&&Math.hypot(im[0].x-im[1].x,im[0].y-im[1].y)<1.38);assert.deepEqual(assignBondOrders(im,cyc(5)),[1.5,1.5,1.5,1.5,1.5],'imidazole');
 assert.deepEqual(assignBondOrders(poly(Array(6).fill('C'),1.397),cyc(6)),Array(6).fill(1.5),'benzene')});
test('ring perception does not fire for puckered rings, aliphatic lengths, or non-ring chains',()=>{
 const chair=poly(['N','C','C','C','C','C'],1.35,0.2),o=assignBondOrders(chair,cyc(6));assert.equal(o[0],1,'N-C bond in a puckered (non-planar) ring is NOT upgraded by the ring pass');assert.equal(o[5],1,'C-N bond likewise');
 assert.deepEqual(assignBondOrders(poly(Array(6).fill('C'),1.54),cyc(6)),Array(6).fill(1),'cyclohexane-length ring');
 const flat=poly(Array(6).fill('C'),1.38,0.04);assert.deepEqual(assignBondOrders(flat,cyc(6)),Array(6).fill(1.5),'slightly noisy flat ring still aromatic');
 const chain=[{element:'C',x:0,y:0,z:0},{element:'N',x:1.34,y:0,z:0},{element:'C',x:2.0,y:1.2,z:0}];assert.deepEqual(assignBondOrders(chain,[[0,1],[1,2]]),[1,1],'amide-like C-N in a chain stays single')});
