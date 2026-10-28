// Idealised multi-residue test models for the analysis tools. Constructed, NOT experimental structures.
// Heavy atoms only: backbone N, CA, C, O plus CB for every non-Gly residue. Side chains beyond CB are omitted and Pro has no ring.
import {peptide,transformAtoms} from './build.mjs';
const A3={A:'ALA',R:'ARG',N:'ASN',D:'ASP',C:'CYS',Q:'GLN',E:'GLU',G:'GLY',H:'HIS',I:'ILE',L:'LEU',K:'LYS',M:'MET',F:'PHE',P:'PRO',S:'SER',T:'THR',W:'TRP',Y:'TYR',V:'VAL'};
export const seq3=s=>[...s].map(c=>{if(!A3[c])throw new RangeError('Unknown one-letter residue code "'+c+'".');return A3[c]});
export const MODEL_SEQ={helixWT:'AEAAAKEAAAKA',helixVar:'GSAEAAAREAAGKA',hlh:'SEELLKKALEEGPDKEELLKKALEE',hairpin:'SWTWENGKWTWK',dimerA:'AEAAAKEAAAKA',dimerB:'AEAAAREAAGKA',ramaTest:'APAAAAAPAAGGAAAA'};
// ideal right-handed alpha-helix, any sequence
export const alphaHelix=seq=>peptide(seq.length,{phi:-57,psi:-47,seq:seq3(seq)});
// helix (res 1-11) - 4-residue loop G-P-D-K (res 12-15) - helix (res 16-25); loop torsions found by numerical search so the helices pack antiparallel without clashes
export function helixLoopHelix(){const n=25,phi=new Array(n).fill(-57),psi=new Array(n).fill(-47);[[79,16],[-63,150],[-59,-43],[-98,7]].forEach(([f,p],k)=>{phi[11+k]=f;psi[11+k]=p});return peptide(n,{phi,psi,seq:seq3(MODEL_SEQ.hlh)})}
// 12-residue beta-hairpin: strands 1-5 and 8-12, type I' turn at Asn6-Gly7. Torsions optimised for cross-strand CA-CA registry and N-H...O=C distances.
export function betaHairpin(){const phi=[-130,-126,-127,-120,-148,66,95,-86,-113,-127,-110,-124],psi=[155,139,131,127,135,16,-8,135,127,131,132,117];return peptide(12,{phi,psi,seq:seq3(MODEL_SEQ.hairpin)})}
// 16-residue peptide whose 14 phi/psi pairs visit every region of the plot, including Gly-only and strained positions
export function ramaTestPeptide(){const phi=[0,-65,-75,-57,60,57,-120,-65,-135,-60,80,75,-110,-90,-63,-120],psi=[130,145,145,-47,-120,45,130,-35,140,-45,-170,25,135,60,-41,130];return peptide(16,{phi,psi,seq:seq3(MODEL_SEQ.ramaTest)})}
// two helical chains A and B: B is a rigid antiparallel copy of an alpha-helix, axis `separation` A away, with its own sequence
const vs=(a,b)=>[a[0]-b[0],a[1]-b[1],a[2]-b[2]],vc=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],vn=a=>{const l=Math.hypot(...a);return a.map(x=>x/l)};
export function antiparallelHelixPair({seqA=MODEL_SEQ.dimerA,seqB=MODEL_SEQ.dimerB,separation=10}={}){if(seqA.length!==seqB.length)throw new RangeError('Chains must have the same length (same backbone).');
 const A=alphaHelix(seqA),B0=alphaHelix(seqB),ca=A.filter(a=>a.name==='CA').map(a=>[a.x,a.y,a.z]),c=[];
 for(let i=0;i+3<ca.length;i++)c.push([0,1,2].map(k=>(ca[i][k]+ca[i+1][k]+ca[i+2][k]+ca[i+3][k])/4)); // 4-CA centroids cancel the helical periodicity
 const dir=vn(vs(c[c.length-1],c[0])),m=[0,1,2].map(k=>c.reduce((s,x)=>s+x[k],0)/c.length),e=[[1,0,0],[0,1,0],[0,0,1]].sort((p,q)=>Math.abs(p[0]*dir[0]+p[1]*dir[1]+p[2]*dir[2])-Math.abs(q[0]*dir[0]+q[1]*dir[1]+q[2]*dir[2]))[0];
 const w=vn(vc(dir,e)),u=vc(dir,w),R=[0,1,2].map(i=>[0,1,2].map(j=>2*u[i]*u[j]-(i===j?1:0))), // pi rotation about u (perpendicular to the axis) flips the axis direction
  Rm=R.map(r=>r[0]*m[0]+r[1]*m[1]+r[2]*m[2]),t=[0,1,2].map(i=>m[i]-Rm[i]+separation*w[i]);
 return[...A,...transformAtoms(B0,{R,t,chain:'B',serialStart:A.length+1})]}
