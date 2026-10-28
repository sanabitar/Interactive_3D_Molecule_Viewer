import test from 'node:test';import assert from 'node:assert/strict';
import {dihedral} from '../src/core/geometry.mjs';import {peptide} from '../src/core/build.mjs';
import {phiPsi,ramaRegion,ramaSummary,RAMA_REGIONS} from '../src/core/analysis.mjs';import {ramaTestPeptide,antiparallelHelixPair} from '../src/core/models.mjs';
const MAP=a=>({...a,element:a.el,resName:a.res,resSeq:a.resi,hetero:false,iCode:''});
const close=(a,b,t,m)=>assert.ok(Math.abs(a-b)<=t,`${m||''} ${a} vs ${b}`);
const cross=(u,v)=>[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]],dot=(u,v)=>u[0]*v[0]+u[1]*v[1]+u[2]*v[2];

// --- Correctness of the angles themselves, against definitions that do not use the code under test ---
test('dihedral sign follows the IUPAC rule: clockwise looking from atom 2 to atom 3 is positive',()=>{
 // B at origin, C on +z, viewer at B looking along +z. With x to the right (y is then down), clockwise on screen = +x turning toward +y.
 const B={x:0,y:0,z:0},C={x:0,y:0,z:1},A={x:1,y:0,z:0};
 for(const t of[60,-60,120,-120,180,0]){const D={x:Math.cos(t*Math.PI/180),y:Math.sin(t*Math.PI/180),z:1};close(((dihedral(A,B,C,D)-t+540)%360)-180,0,1e-6,'t='+t)}
 });
test('phi=-57, psi=-47 builds a RIGHT-handed helix and phi=+57, psi=+47 a LEFT-handed one (real-world handedness check)',()=>{
 const hand=(phi,psi)=>{const ca=peptide(12,{phi,psi}).filter(a=>a.name==='CA'),v=i=>[ca[i+1].x-ca[i].x,ca[i+1].y-ca[i].y,ca[i+1].z-ca[i].z];let s=0;for(let i=0;i<ca.length-3;i++)s+=dot(cross(v(i),v(i+1)),v(i+2));return s};
 assert.ok(hand(-57,-47)>0,'alpha helix must be right-handed');assert.ok(hand(57,47)<0,'phi,psi inverted must be left-handed')});
test('alpha-helix i->i+4 C=O...N distance is a hydrogen bond (~3.1 A) at phi=-57, psi=-47',()=>{
 const P=peptide(12,{phi:-57,psi:-47}),g=(r,n)=>P.find(a=>a.resi===r&&a.name===n),d=(a,b)=>Math.hypot(a.x-b.x,a.y-b.y,a.z-b.z);
 for(let i=1;i<=7;i++)close(d(g(i,'O'),g(i+4,'N')),3.1,.2,'O'+i)});

// --- Region classification ---
test('ramaRegion: canonical conformations land in the expected schematic region',()=>{
 const c=(f,p,r='ALA')=>ramaRegion(f,p,r);
 assert.equal(c(-57,-47),'alphaR');assert.equal(c(-49,-26),'alphaR');assert.equal(c(-120,130),'beta');assert.equal(c(-75,145),'beta');assert.equal(c(57,47),'alphaL');
 assert.equal(c(-120,-175),'beta','psi wrap-around near -180 is still the beta region');assert.equal(c(-135,180),'beta');assert.equal(c(180,180),'beta','+180 and -180 are the same angle');
 assert.equal(c(60,-120),'other','positive phi with strongly negative psi is strained for a non-Gly residue');assert.equal(c(-90,60),'other','bridge region lies between the shaded cores');
 assert.equal(c(60,-120,'GLY'),'gly','Gly is never classified against the non-Gly boxes');assert.equal(c(-57,-47,'PRO'),'alphaR');
 assert.equal(RAMA_REGIONS.length,3)});
test('phiPsi attaches region to each point and ramaSummary counts them',()=>{
 const A=peptide(12,{phi:-57,psi:-47}).map(MAP),R=phiPsi(A);assert.equal(R.length,10);assert.ok(R.every(r=>r.region==='alphaR'));
 const S=ramaSummary(R);assert.equal(S.total,10);assert.equal(S.counts.alphaR,10);assert.equal(S.outside.length,0)});

// --- Test peptide visits every category, and the tool reports what is known to be in it ---
test('Ramachandran test peptide: 14 points covering all regions, two Gly, two Pro, and the known strained residue',()=>{
 const A=ramaTestPeptide().map(MAP),R=phiPsi(A),S=ramaSummary(R);assert.equal(R.length,14,'termini have no phi/psi pair');
 assert.equal(S.counts.gly,2);assert.ok(S.counts.alphaR>=3&&S.counts.beta>=4&&S.counts.alphaL===1,JSON.stringify(S.counts));
 assert.deepEqual(S.outside.map(r=>r.residue).sort(),['A:14','A:5'],'strained Ala 5 (60,-120) and bridge-region Ala 14 (-90,60)');
 assert.equal(R.filter(r=>r.resName==='PRO').length,2);
 const r5=R.find(r=>r.residue==='A:5');close(r5.phi,60,1e-6);close(r5.psi,-120,1e-6);
 const g11=R.find(r=>r.residue==='A:11');assert.equal(g11.region,'gly');close(g11.phi,80,1e-6);close(g11.psi,-170,1e-6)});

// --- Chain handling ---
test('phiPsi: chain filter, blank chain, chain breaks, hetero atoms',()=>{
 const D=antiparallelHelixPair().map(MAP);assert.equal(phiPsi(D).length,20);assert.equal(phiPsi(D,'A').length,10);assert.equal(phiPsi(D,'B').length,10);assert.equal(phiPsi(D,'Z').length,0);
 const noChain=peptide(6,{phi:-57,psi:-47}).map(a=>({...MAP(a),chain:''}));assert.equal(phiPsi(noChain,'').length,4);
 const het=peptide(6,{phi:-57,psi:-47}).map(a=>({...MAP(a),hetero:true}));assert.equal(phiPsi(het).length,0,'HETATM residues (other than MSE) are skipped');
 const broken=peptide(8,{phi:-57,psi:-47}).map(MAP);for(const a of broken)if(a.resSeq>4)a.x+=30;assert.deepEqual(phiPsi(broken).map(r=>r.residue),['A:2','A:3','A:6','A:7'],'C4-N5 > 2 A is a chain break: residues 4 and 5 lose a neighbour, 2-3 and 6-7 keep both')});
