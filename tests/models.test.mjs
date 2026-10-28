import test from 'node:test';import assert from 'node:assert/strict';
import {peptide,transformAtoms} from '../src/core/build.mjs';import {phiPsi,alignGlobal,chainSequence,caContactMap,compareStructures} from '../src/core/analysis.mjs';
import {MODEL_SEQ,seq3,alphaHelix,helixLoopHelix,betaHairpin,ramaTestPeptide,antiparallelHelixPair} from '../src/core/models.mjs';import {molProps} from '../src/core/properties.mjs';
const MAP=a=>({...a,element:a.el,resName:a.res,resSeq:a.resi,hetero:false,iCode:''});
const close=(a,b,t,m)=>assert.ok(Math.abs(a-b)<=t,`${m||''} ${a} vs ${b}`);
const d=(a,b)=>Math.hypot(a.x-b.x,a.y-b.y,a.z-b.z),sub=(a,b)=>[a[0]-b[0],a[1]-b[1],a[2]-b[2]],dot=(u,v)=>u[0]*v[0]+u[1]*v[1]+u[2]*v[2],len=u=>Math.hypot(...u),nrm=u=>u.map(x=>x/len(u));
const clearance=A=>{let m=1e9;for(let i=0;i<A.length;i++)for(let j=i+1;j<A.length;j++){const a=A[i],b=A[j];if(a.chain===b.chain&&Math.abs(a.resi-b.resi)<2)continue;m=Math.min(m,d(a,b))}return m};
const axis=(A,lo,hi)=>{const ca=A.filter(a=>a.name==='CA'&&a.resi>=lo&&a.resi<=hi).map(a=>[a.x,a.y,a.z]),c=[];for(let i=0;i+3<ca.length;i++)c.push([0,1,2].map(k=>(ca[i][k]+ca[i+1][k]+ca[i+2][k]+ca[i+3][k])/4));return{dir:nrm(sub(c[c.length-1],c[0])),mid:[0,1,2].map(k=>c.reduce((s,x)=>s+x[k],0)/c.length)}};
const axisLineDistance=(a,b)=>{const v=sub(b.mid,a.mid),k=dot(v,a.dir);return len([v[0]-k*a.dir[0],v[1]-k*a.dir[1],v[2]-k*a.dir[2]])};

test('peptide(): sequences, per-residue torsions, validation, Gly has no CB, backward compatible',()=>{
 const phi=[0,-60,-120,60,-75],psi=[130,-40,130,40,145],P=peptide(5,{phi,psi,seq:['ALA','GLY','SER','PRO','LYS']}),A=P.map(MAP),R=phiPsi(A);
 assert.equal(P.length,5*5-1);assert.deepEqual([...new Set(P.map(a=>a.res))],['ALA','GLY','SER','PRO','LYS']);assert.ok(!P.some(a=>a.res==='GLY'&&a.name==='CB'));assert.ok(P.filter(a=>a.res!=='GLY').every((a,_,all)=>all.filter(x=>x.name==='CB').length===4));
 R.forEach((r,k)=>{close(r.phi,phi[k+1],1e-6,'phi'+k);close(r.psi,psi[k+1],1e-6,'psi'+k)});
 assert.throws(()=>peptide(5,{phi,psi,seq:['ALA']}),RangeError);assert.throws(()=>peptide(5,{phi:[1,2],psi}),RangeError);assert.throws(()=>peptide(5,{phi:NaN,psi:1}),RangeError);
 assert.equal(peptide(12,{phi:-57,psi:-47}).length,60,'poly-Ala default unchanged');assert.equal(peptide(4,{phi:-57,psi:-47,resName:'GLY'}).length,16)});
test('seq3 and transformAtoms',()=>{assert.deepEqual(seq3('AGW'),['ALA','GLY','TRP']);assert.throws(()=>seq3('AXB'),RangeError);
 const P=peptide(3,{phi:-57,psi:-47}),R=[[0,-1,0],[1,0,0],[0,0,1]],Q=transformAtoms(P,{R,t:[1,2,3],chain:'Z',serialStart:100});
 close(Q[0].x,-P[0].y+1,1e-12);close(Q[0].y,P[0].x+2,1e-12);close(Q[0].z,P[0].z+3,1e-12);assert.equal(Q[0].chain,'Z');assert.equal(Q[0].serial,100);close(d(Q[0],Q[5]),d(P[0],P[5]),1e-9,'rigid')});

test('helix-loop-helix: 25 residues, antiparallel packed helices, loop torsions, no clashes',()=>{
 const P=helixLoopHelix(),A=P.map(MAP);assert.equal(P.filter(a=>a.name==='CA').length,25);assert.equal(chainSequence(A,'A'),MODEL_SEQ.hlh);
 const a1=axis(P,3,10),a2=axis(P,17,24),ang=Math.acos(dot(a1.dir,a2.dir))*180/Math.PI,dist=axisLineDistance(a1,a2);
 assert.ok(ang>150,'axes antiparallel: '+ang);assert.ok(dist>6.5&&dist<11.5,'axis separation '+dist);assert.ok(clearance(P)>2.9,'clearance '+clearance(P));
 const R=phiPsi(A);const loop=R.filter(r=>+r.residue.split(':')[1]>=12&&+r.residue.split(':')[1]<=15);assert.equal(loop.length,4);assert.equal(loop[0].resName,'GLY');assert.equal(loop[1].resName,'PRO');
 const M=caContactMap(A,'A',8),i1=0,j2=24;assert.equal(M.n,25);let cross=0;for(let i=0;i<11;i++)for(let j=15;j<25;j++)cross+=M.matrix[i*M.n+j];assert.ok(cross>=15,'helix-helix CA contacts within 8 A: '+cross)});

test('beta-hairpin: cross-strand CA-CA registry, backbone H-bonds, type I\' turn, strand torsions in beta region',()=>{
 const P=betaHairpin(),A=P.map(MAP),g=(r,n)=>P.find(a=>a.resi===r&&a.name===n);assert.equal(P.length,59);assert.equal(chainSequence(A,'A'),MODEL_SEQ.hairpin);
 for(const[i,j,t]of[[5,8,5.4],[4,9,4.5],[3,10,5.4],[2,11,4.5],[1,12,5.4]])close(d(g(i,'CA'),g(j,'CA')),t,.6,`CA${i}-CA${j}`);
 for(const[a,x,b,y]of[[5,'O',8,'N'],[3,'N',10,'O'],[3,'O',10,'N'],[1,'N',12,'O'],[1,'O',12,'N']])assert.ok(d(g(a,x),g(b,y))<3.3,`${a}${x}-${b}${y} ${d(g(a,x),g(b,y))}`);
 const R=phiPsi(A),t6=R.find(r=>r.residue==='A:6'),t7=R.find(r=>r.residue==='A:7');close(t6.phi,66,1e-6);close(t6.psi,16,1e-6);close(t7.phi,95,1e-6);close(t7.psi,-8,1e-6);assert.equal(t7.resName,'GLY');
 for(const r of R.filter(r=>![6,7].includes(+r.residue.split(':')[1])))assert.equal(r.region,'beta',r.residue);assert.ok(clearance(P)>2.9)});

test('ramaTestPeptide geometry: 16 residues, no clashes',()=>{const P=ramaTestPeptide();assert.equal(P.length,78);assert.equal(chainSequence(P.map(MAP),'A'),MODEL_SEQ.ramaTest);assert.ok(clearance(P)>2.9)});

test('EAAAK helix pair for alignment and RMSD: known identity and near-zero RMSD',()=>{
 const W=alphaHelix(MODEL_SEQ.helixWT).map(MAP),V=alphaHelix(MODEL_SEQ.helixVar).map(MAP);assert.equal(W.length,12*5-0);assert.equal(chainSequence(W,'A'),'AEAAAKEAAAKA');assert.equal(chainSequence(V,'A'),'GSAEAAAREAAGKA');
 const al=alignGlobal(chainSequence(W,'A'),chainSequence(V,'A'));assert.equal(al.identical,10);assert.equal(al.gaps,2);assert.equal(al.alignedPairs,12);close(al.coverageA,1,1e-12);close(al.coverageB,12/14,1e-12);
 const r=compareStructures(W,V,{mode:'ca'});assert.equal(r.nPairs,12);assert.ok(r.rmsd<1e-6,'a register shift of an ideal helix superposes exactly: '+r.rmsd)});

test('antiparallel helix pair: two chains, antiparallel axes ~10 A apart, no clashes, B differs at 2 positions',()=>{
 const P=antiparallelHelixPair(),A=P.map(MAP);assert.equal(P.length,59+60);assert.deepEqual([...new Set(P.map(a=>a.chain))],['A','B']);assert.equal(new Set(P.map(a=>a.serial)).size,P.length,'unique serials');
 const sa=chainSequence(A,'A'),sb=chainSequence(A,'B');assert.equal(sa,MODEL_SEQ.dimerA);assert.equal(sb,MODEL_SEQ.dimerB);assert.equal([...sa].filter((c,i)=>c!==sb[i]).length,2);
 const aA=axis(P.filter(a=>a.chain==='A'),1,12),aB=axis(P.filter(a=>a.chain==='B'),1,12);close(dot(aA.dir,aB.dir),-1,1e-6,'antiparallel');close(axisLineDistance(aA,aB),10,1e-6,'separation');assert.ok(clearance(P)>3.0,'clearance '+clearance(P));
 const r=compareStructures(A,A,{mode:'ca',chainA:'A',chainB:'B'});assert.equal(r.nPairs,12);assert.ok(r.rmsd<1e-6,'chain B is a rigid copy of chain A: '+r.rmsd)});
test('every model has finite coordinates and a sensible formula',()=>{for(const P of[helixLoopHelix(),betaHairpin(),ramaTestPeptide(),antiparallelHelixPair(),alphaHelix(MODEL_SEQ.helixVar)]){assert.ok(P.every(a=>[a.x,a.y,a.z].every(Number.isFinite)));const p=molProps(P.map(MAP));assert.ok(p.atomCount===P.length&&/^C/.test(p.formula))}});
