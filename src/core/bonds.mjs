// Bond-order and hybridization inference from bond lengths + neighbours (heuristic; PDB/mmCIF coordinates do not store bond orders),
// and a consistency check of measured angles against the ideal angle of each atom's hybridization.
import {distance,angle} from './geometry.mjs';
const elOf=a=>a.element||a.el;
function lenOrder(x,y,d){switch([x,y].sort().join('-')){case'C-C':return d<1.25?3:d<1.37?2:d<1.43?1.5:1;case'C-N':return d<1.22?3:d<1.30?2:1;case'C-O':return d<1.14?3:d<1.30?2:1;case'N-N':return d<1.15?3:d<1.30?2:1;case'O-O':return d<1.30?2:1;default:return 1}}
// Ring perception (length path only): a 5- or 6-membered ring of C/N atoms whose bonds ALL have aromatic-range lengths (1.30-1.43 A) and that is planar
// (every atom within 0.12 A of the mean plane) is treated as aromatic: all its bonds become 1.5. Fixes mixed depictions such as pyridine (C-N drawn single, C-C aromatic).
// It is a geometric heuristic: it cannot recognise non-aromatic planar rings with aromatic-range lengths, nor heteroatoms O/S in rings.
const RING_LO=1.30,RING_HI=1.43,PLANAR_TOL=0.12;
function planarRing(atoms,ring){const P=ring.map(i=>atoms[i]),n=P.length,c=['x','y','z'].map(k=>P.reduce((s,a)=>s+a[k],0)/n);let nx=0,ny=0,nz=0;
 for(let i=0;i<n;i++){const a=P[i],b=P[(i+1)%n];nx+=(a.y-b.y)*(a.z+b.z);ny+=(a.z-b.z)*(a.x+b.x);nz+=(a.x-b.x)*(a.y+b.y)}
 const L=Math.hypot(nx,ny,nz);if(!(L>1e-9))return false;return P.every(a=>Math.abs(((a.x-c[0])*nx+(a.y-c[1])*ny+(a.z-c[2])*nz)/L)<=PLANAR_TOL)}
function ringPass(atoms,bonds,ord){const cand=new Map(),add=(u,v,k)=>{if(!cand.has(u))cand.set(u,[]);cand.get(u).push([v,k])};
 bonds.forEach(([i,j],k)=>{const ea=elOf(atoms[i]),eb=elOf(atoms[j]);if((ea!=='C'&&ea!=='N')||(eb!=='C'&&eb!=='N'))return;const d=distance(atoms[i],atoms[j]);if(d<RING_LO||d>RING_HI)return;add(i,j,k);add(j,i,k)});
 const seen=new Set();for(const s of cand.keys()){const path=[s],bk=[];
  const dfs=u=>{for(const[v,k]of cand.get(u)){if(v===s&&path.length>=5&&path.length<=6){const id=[...path].sort((x,y)=>x-y).join(',');if(!seen.has(id)){seen.add(id);if(planarRing(atoms,path))[...bk,k].forEach(b=>{ord[b]=1.5})}continue}
   if(v<s||path.length>=6||path.includes(v))continue;path.push(v);bk.push(k);dfs(v);path.pop();bk.pop()}};dfs(s)}}
// orders: 1, 2, 3, or 1.5 (aromatic / delocalised). Aligned with `bonds` ([i,j] atom-index pairs).
export function assignBondOrders(atoms,bonds){const adj=atoms.map(()=>[]),key=new Map();bonds.forEach(([i,j],k)=>{adj[i].push(j);adj[j].push(i);key.set(i+','+j,k);key.set(j+','+i,k)});
 const ord=bonds.map(([i,j])=>lenOrder(elOf(atoms[i]),elOf(atoms[j]),distance(atoms[i],atoms[j])));
 ringPass(atoms,bonds,ord);
 atoms.forEach((a,c)=>{if(elOf(a)!=='C'||adj[c].length<3)return;const t=adj[c].filter(o=>elOf(atoms[o])==='O'&&adj[o].length===1&&ord[key.get(c+','+o)]===2);if(t.length>=2)t.forEach(o=>{ord[key.get(c+','+o)]=1.5})}); // carboxylate-type resonance
 return ord}
// Hybridization for C, N, O only (others null). N/O attached to an unsaturated carbon are treated as sp2 (conjugation, e.g. amide N, phenol O).
export function hybridization(atoms,bonds,orders){const adj=atoms.map(()=>[]);bonds.forEach(([i,j],k)=>{adj[i].push([j,orders[k]]);adj[j].push([i,orders[k]])});
 const unsatC=i=>elOf(atoms[i])==='C'&&adj[i].some(([,o])=>o>1);
 return atoms.map((a,i)=>{const e=elOf(a),nb=adj[i];if(!nb.length)return null;const dbl=nb.filter(([,o])=>o===2).length,trp=nb.some(([,o])=>o===3),aro=nb.some(([,o])=>o===1.5);
  if(e==='C')return trp||dbl>=2?'sp':dbl||aro?'sp2':'sp3';if(e==='N')return trp?'sp':dbl||aro?'sp2':nb.some(([j])=>unsatC(j))?'sp2':'sp3';
  if(e==='O')return trp?'sp':dbl||aro?'sp2':nb.some(([j])=>unsatC(j))?'sp2':'sp3';return null})}
const IDEAL={sp:180,sp2:120,sp3:109.5},TOL={sp:5,sp2:12,sp3:16};
export function geometryCheck(atoms,bonds,hyb){const adj=atoms.map(()=>[]);bonds.forEach(([i,j])=>{adj[i].push(j);adj[j].push(i)});const outliers=[];let max=0,n=0;
 atoms.forEach((a,c)=>{const h=hyb[c];if(!h||adj[c].length<2)return;const nb=adj[c];for(let x=0;x<nb.length;x++)for(let y=x+1;y<nb.length;y++){const ang=angle(atoms[nb[x]],a,atoms[nb[y]]),dev=Math.abs(ang-IDEAL[h]);n++;if(dev>max)max=dev;
  if(dev>TOL[h])outliers.push({atom:c,hybridization:h,neighbors:[nb[x],nb[y]],angle:ang,ideal:IDEAL[h],deviation:dev})}});return{checked:n,maxDeviation:max,outliers,tolerances:TOL}}

// Authoritative-by-template bond orders for the 20 standard amino acids (residue dictionary), length heuristic otherwise.
const STD=new Set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(' '));
const AROM={PHE:['CG','CD1','CD2','CE1','CE2','CZ'],TYR:['CG','CD1','CD2','CE1','CE2','CZ'],HIS:['CG','ND1','CD2','CE1','NE2'],TRP:['CG','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2']};
const DELOC={ASP:['CG','OD1','OD2'],GLU:['CD','OE1','OE2'],ARG:['CZ','NH1','NH2','NE']},DBL={ASN:[['CG','OD1']],GLN:[['CD','OE1']]};
// returns {orders, source}: source[k] is 'template' (standard residue dictionary; delocalised groups = 1.5) or 'length' (inferred from bond length only).
export function assignBondOrdersDetailed(atoms,bonds){const lenOrd=assignBondOrders(atoms,bonds),key=a=>a.chain+'|'+a.resSeq+'|'+(a.iCode||''),res=new Map(),orders=[],source=[];
 const std=a=>!a.hetero&&STD.has(a.resName);atoms.forEach(a=>{if(std(a)){const k=key(a);if(!res.has(k))res.set(k,{n:a.resName,names:new Set()});res.get(k).names.add(a.name)}});
 bonds.forEach(([i,j],k)=>{const a=atoms[i],b=atoms[j];if(std(a)&&std(b)){let o=1;if(key(a)===key(b)){const r=res.get(key(a)),n=r.n,A=a.name,B=b.name,has=(x,y)=>(A===x&&B===y)||(A===y&&B===x);
   if(has('C','O'))o=r.names.has('OXT')?1.5:2;else if(has('C','OXT'))o=1.5;else if(AROM[n]&&AROM[n].includes(A)&&AROM[n].includes(B))o=1.5;else if(DELOC[n]&&DELOC[n].includes(A)&&DELOC[n].includes(B))o=1.5;else if(DBL[n]&&DBL[n].some(([x,y])=>has(x,y)))o=2}
  orders.push(o);source.push('template')}else{orders.push(lenOrd[k]);source.push('length')}});return{orders,source}}
