import {superpose} from './superpose.mjs';
export const reflect = atoms => atoms.map(a=>({...a,x:-a.x}));
// Enumerate graph automorphisms; only proper rotations are allowed in fitting.
export function checkMirror(atoms,bonds,orders=[],{limit=50000,tolerance=0.02}={}) {
 if(atoms.length<3 || atoms.length>80) throw Error('Use a molecule with 3–80 atoms for this teaching tool.');
 if(atoms.some(a=>![a.x,a.y,a.z].every(Number.isFinite))) throw Error('Invalid coordinates.');
 const n=atoms.length,edge=Array.from({length:n},()=>Array(n).fill(0));
 bonds.forEach(([i,j],k)=>edge[i][j]=edge[j][i]=orders[k]||1);
 const el=a=>a.element||a.el;
 const sig=atoms.map((a,i)=>el(a)+'|'+edge[i].map((v,j)=>v?`${el(atoms[j])}:${v}`:'').filter(Boolean).sort().join(','));
 const candidates=sig.map(s=>sig.flatMap((t,j)=>s===t?[j]:[]));
 const sequence=Array.from({length:n},(_,i)=>i).sort((a,b)=>candidates[a].length-candidates[b].length);
 const mirror=reflect(atoms),map=Array(n).fill(-1),used=new Set();
 let best=null,visits=0,complete=true,matched=false;
 function visit(depth){
  if(matched)return;
  if(++visits>limit){complete=false;return;}
  if(depth===n){const fit=superpose(map.map(j=>mirror[j]),atoms);if(!best||fit.rmsd<best.rmsd)best={...fit,mapping:[...map]};if(fit.rmsd<=tolerance)matched=true;return;}
  const i=sequence[depth];
  for(const j of candidates[i]){if(used.has(j))continue;let ok=true;for(let k=0;k<n;k++)if(map[k]>=0&&edge[i][k]!==edge[j][map[k]]){ok=false;break;}if(!ok)continue;
   map[i]=j;used.add(j);visit(depth+1);used.delete(j);map[i]=-1;if(matched||!complete)return;
  }
 }
 visit(0);
 return {matched,complete,best,visits,tolerance,mirror};
}
