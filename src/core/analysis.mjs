// Sequence alignment, phi/psi, CA contact map. Pure functions. Angles in degrees, distances in Angstrom.
import {dihedral,distance} from './geometry.mjs';import {superpose} from './superpose.mjs';import {residueKey} from './structure.mjs';
const AA={ALA:'A',ARG:'R',ASN:'N',ASP:'D',CYS:'C',GLN:'Q',GLU:'E',GLY:'G',HIS:'H',ILE:'I',LEU:'L',LYS:'K',MET:'M',PHE:'F',PRO:'P',SER:'S',THR:'T',TRP:'W',TYR:'Y',VAL:'V',MSE:'M'};
export function residues(atoms,chain){const m=new Map();for(const a of atoms){if(a.hetero&&a.resName!=='MSE')continue;if(chain!==undefined&&a.chain!==chain)continue;const k=residueKey(a);let r=m.get(k);if(!r){r={key:k,resName:a.resName,chain:a.chain,atoms:{}};m.set(k,r)}if(!(a.name in r.atoms))r.atoms[a.name]=a}return[...m.values()]}
export const chainSequence=(atoms,chain)=>residues(atoms,chain).map(r=>AA[r.resName]||'X').join('');
// Needleman-Wunsch global alignment, linear gap penalty, simple match/mismatch scoring (no substitution matrix).
export function alignGlobal(a,b,{match=2,mismatch=-1,gap=-2}={}){a=String(a).replace(/\s/g,'').toUpperCase();b=String(b).replace(/\s/g,'').toUpperCase();
 if(!a||!b)throw new RangeError('Both sequences must be non-empty.');if(!/^[A-Z]+$/.test(a+b))throw new RangeError('Sequences may contain letters only.');
 const n=a.length,m=b.length,W=m+1;if(n*m>4e6)throw new RangeError('Sequences too long for in-browser alignment (limit n*m = 4,000,000).');
 const S=new Float64Array((n+1)*W);for(let i=1;i<=n;i++)S[i*W]=i*gap;for(let j=1;j<=m;j++)S[j]=j*gap;
 for(let i=1;i<=n;i++)for(let j=1;j<=m;j++)S[i*W+j]=Math.max(S[(i-1)*W+j-1]+(a[i-1]===b[j-1]?match:mismatch),S[(i-1)*W+j]+gap,S[i*W+j-1]+gap);
 let i=n,j=m,x='',y='';while(i>0||j>0){if(i>0&&j>0&&S[i*W+j]===S[(i-1)*W+j-1]+(a[i-1]===b[j-1]?match:mismatch)){x=a[--i]+x;y=b[--j]+y}else if(i>0&&S[i*W+j]===S[(i-1)*W+j]+gap){x=a[--i]+x;y='-'+y}else{x='-'+x;y=b[--j]+y}}
 let identical=0,pairs=0;for(let k=0;k<x.length;k++)if(x[k]!=='-'&&y[k]!=='-'){pairs++;if(x[k]===y[k])identical++}
 return{a:x,b:y,score:S[n*W+m],length:x.length,identical,alignedPairs:pairs,gaps:x.length-pairs,identityPct:100*identical/x.length,coverageA:pairs/n,coverageB:pairs/m,params:{algorithm:'Needleman-Wunsch (global)',match,mismatch,gap,matrix:'none'}}}
// Schematic reference regions for a Ramachandran plot (degrees): generous boxes [phiMin,phiMax,psiMin,psiMax] around the commonly populated areas.
// They are NOT validation contours (those are density-derived from large reference sets) and the true allowed area is wider (e.g. the bridge region psi ~ 20-90 deg).
export const RAMA_REGIONS=[{id:'alphaR',label:'α-helix (right-handed)',boxes:[[-150,-35,-75,20]]},{id:'beta',label:'β-strand / PPII',boxes:[[-180,-50,90,180],[-180,-60,-180,-165]]},{id:'alphaL',label:'left-handed helix',boxes:[[35,90,-10,80]]}];
const wrap180=v=>v>180-1e-9?v-360:v; // dihedral() returns (-180,180]; +180 and -180 are the same angle
// Glycine is not classified ('gly'): it lacks a side chain, so its allowed area is much larger and roughly symmetric. Anything else outside the boxes is 'other' (loops, turns, bridge region, or strained).
export function ramaRegion(phi,psi,resName){if(resName==='GLY')return'gly';phi=wrap180(phi);psi=wrap180(psi);for(const r of RAMA_REGIONS)for(const[a,b,c,d]of r.boxes)if(phi>=a&&phi<=b&&psi>=c&&psi<=d)return r.id;return'other'}
export function ramaSummary(points){const counts={alphaR:0,beta:0,alphaL:0,other:0,gly:0},outside=[];for(const p of points){counts[p.region]++;if(p.region==='other')outside.push(p)}return{total:points.length,counts,outside}}
export function phiPsi(atoms,chain){const R=residues(atoms,chain),out=[];for(let i=0;i<R.length;i++){const r=R[i],p=R[i-1],n=R[i+1],a=r.atoms;if(!a.N||!a.CA||!a.C)continue;let phi=null,psi=null;
 if(p&&p.chain===r.chain&&p.atoms.C&&distance(p.atoms.C,a.N)<2.0)phi=dihedral(p.atoms.C,a.N,a.CA,a.C);
 if(n&&n.chain===r.chain&&n.atoms.N&&distance(a.C,n.atoms.N)<2.0)psi=dihedral(a.N,a.CA,a.C,n.atoms.N);
 if(phi!==null&&psi!==null)out.push({residue:r.key,resName:r.resName,phi,psi,region:ramaRegion(phi,psi,r.resName)})}return out}
export function caContactMap(atoms,chain,cutoff=8){if(!(cutoff>0))throw new RangeError('cutoff must be > 0 Angstrom');const R=residues(atoms,chain).filter(r=>r.atoms.CA),n=R.length;
 if(!n)throw new RangeError('No CA atoms found for this chain.');if(n>2000)throw new RangeError('Too many residues (limit 2000).');const M=new Uint8Array(n*n);
 for(let i=0;i<n;i++)for(let j=i;j<n;j++)if(distance(R[i].atoms.CA,R[j].atoms.CA)<=cutoff)M[i*n+j]=M[j*n+i]=1;return{keys:R.map(r=>r.key),n,matrix:M,cutoff}}
// Pair atoms between structures A (reference) and B (moved), superpose B onto A, report RMSD over the paired atoms.
export function compareStructures(A,B,{mode='ca',chainA,chainB,identicalOnly=false}={}){let pa,pb,aln=null,method;
 if(mode==='ca'){const ra=residues(A,chainA).filter(r=>r.atoms.CA),rb=residues(B,chainB).filter(r=>r.atoms.CA);if(ra.length<3||rb.length<3)throw new RangeError('Each chain needs at least 3 residues with CA atoms.');
  const sq=R=>R.map(r=>AA[r.resName]||'X').join('');aln=alignGlobal(sq(ra),sq(rb));pa=[];pb=[];let i=0,j=0;
  for(let k=0;k<aln.length;k++){const x=aln.a[k],y=aln.b[k];if(x!=='-'&&y!=='-'&&(!identicalOnly||x===y)){pa.push(ra[i].atoms.CA);pb.push(rb[j].atoms.CA)}if(x!=='-')i++;if(y!=='-')j++}
  method='CA atoms paired by sequence alignment'+(identicalOnly?' (identical residues only)':'')}
 else{pa=A.filter(a=>a.element!=='H');pb=B.filter(a=>a.element!=='H');if(pa.length!==pb.length||pa.some((a,i)=>a.element!==pb[i].element))throw new RangeError('Atom-order pairing needs identical heavy-atom element sequences; otherwise the comparison is not valid.');method='Heavy atoms paired by file order'}
 const s=superpose(pb,pa),key=a=>residueKey(a)+':'+a.name;
 return{method,params:{mode,chainA,chainB,identicalOnly,fit:'Horn/Kabsch least squares, rotation+translation, equal weights, no outlier rejection'},nPairs:pa.length,rmsd:s.rmsd,maxDeviation:s.deviations.reduce((m,d)=>d>m?d:m,0),rotation:s.rotation,
  alignment:aln&&{length:aln.length,identityPct:aln.identityPct,coverageA:aln.coverageA,coverageB:aln.coverageB},pairs:pa.map((a,i)=>({a:key(a),b:key(pb[i]),d:s.deviations[i]})),pa,pb,moved:pb.map(s.apply)}}
