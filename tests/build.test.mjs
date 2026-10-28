import test from 'node:test';import assert from 'node:assert/strict';
import {zmatrix,peptide} from '../src/core/build.mjs';import {distance,angle,dihedral} from '../src/core/geometry.mjs';import {phiPsi} from '../src/core/analysis.mjs';
const close=(a,b,t,m)=>assert.ok(Math.abs(a-b)<t,`${m}: ${a} vs ${b}`);
test('zmatrix reproduces specified distance, angle and signed dihedral',()=>{for(const tor of[60,-60,180,-122.6,15]){const A=zmatrix([['C'],['C',0,1.5],['C',1,1.4,0,112],['C',2,1.3,1,105,0,tor]]);
 close(distance(A[2],A[1]),1.4,1e-9,'d');close(angle(A[0],A[1],A[2]),112,1e-9,'ang');close(distance(A[3],A[2]),1.3,1e-9,'d3');close(angle(A[1],A[2],A[3]),105,1e-9,'ang3');close(dihedral(A[0],A[1],A[2],A[3]),tor===180?(dihedral(A[0],A[1],A[2],A[3])>0?180:-180):tor,1e-8,'dihedral '+tor)}});
const MAP=a=>({...a,element:a.el,resName:a.res,resSeq:a.resi,iCode:'',hetero:false});
test('peptide models: phi/psi read back by independent code equal the construction values',()=>{const T={alpha:[-57,-47],beta:[-120,130],ppii:[-75,145],'310':[-49,-26],left:[57,47]};
 for(const[k,[phi,psi]]of Object.entries(T)){const A=peptide(12,{phi,psi}).map(MAP),R=phiPsi(A);assert.equal(A.length,60);assert.equal(R.length,10,k);for(const r of R){close(r.phi,phi,1e-6,k+' phi');close(r.psi,psi,1e-6,k+' psi')}}});
test('peptide geometry: bond lengths, trans peptide bond, no clashes between non-neighbours',()=>{const A=peptide(12,{phi:-57,psi:-47}).map(MAP),g=(r,n)=>A.find(a=>a.resSeq===r&&a.name===n);
 close(distance(g(3,'N'),g(3,'CA')),1.458,1e-9,'N-CA');close(distance(g(3,'C'),g(4,'N')),1.329,1e-9,'C-N');close(distance(g(3,'C'),g(3,'O')),1.231,1e-9,'C=O');close(angle(g(3,'N'),g(3,'CA'),g(3,'C')),111.2,1e-9,'N-CA-C');
 assert.ok(Math.abs(Math.abs(dihedral(g(3,'CA'),g(3,'C'),g(4,'N'),g(4,'CA')))-180)<1e-6,'omega');let min=9;for(let i=0;i<A.length;i++)for(let j=i+1;j<A.length;j++){if(Math.abs(A[i].resSeq-A[j].resSeq)>1)min=Math.min(min,distance(A[i],A[j]))}assert.ok(min>2.0,'min non-neighbour distance '+min);
 // alpha helix rise per residue should be ~1.5 A along the axis: end-to-end CA distance over 11 residues ~ 16-17 A
 const d=distance(g(1,'CA'),g(12,'CA'));assert.ok(d>15&&d<18.5,'helix length '+d)});
test('L-amino acid chirality by the CORN rule',()=>{const A=peptide(3,{phi:-57,psi:-47}).map(MAP),g=n=>A.find(a=>a.resSeq===2&&a.name===n),v=n=>[g(n).x-g('CA').x,g(n).y-g('CA').y,g(n).z-g('CA').z],
 c=v('C'),r=v('CB'),n=v('N'),tp=c[0]*(r[1]*n[2]-r[2]*n[1])-c[1]*(r[0]*n[2]-r[2]*n[0])+c[2]*(r[0]*n[1]-r[1]*n[0]);assert.ok(tp>0,'CORN triple product '+tp+' (>0 means L)')});
