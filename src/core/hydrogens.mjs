// Hydrogen bookkeeping for small molecules. Pure functions, no DOM. Units: Angstrom and degrees.
// Input is a heavy-atom skeleton plus KEKULE bond orders (1, 2 or 3; never 1.5), so that the number of hydrogens follows from valence:
// C 4, N 3, O 2, S 2 (charged N+ 4, O- 1, ...). Positions are idealised (standard X-H lengths, tetrahedral / trigonal / linear directions,
// staggered rotors); they are not energy-minimised.
const R=Math.PI/180;
const XH={C:1.09,N:1.01,O:.96,S:1.34};          // X-H bond lengths (Angstrom)
const VALENCE={C:4,N:3,O:2,S:2};
const add=(a,b)=>[a[0]+b[0],a[1]+b[1],a[2]+b[2]],sub=(a,b)=>[a[0]-b[0],a[1]-b[1],a[2]-b[2]],mul=(a,s)=>[a[0]*s,a[1]*s,a[2]*s];
const dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2],cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],len=a=>Math.hypot(a[0],a[1],a[2]);
const unit=a=>{const l=len(a);return l>1e-12?mul(a,1/l):[0,0,0]};
// unit vector perpendicular to u, as close to `hint` as possible when one is given
function perp(u,hint){let q=hint?sub(hint,mul(u,dot(hint,u))):null;if(!q||len(q)<1e-6){const ax=Math.abs(u[0])<.9?[1,0,0]:[0,1,0];q=sub(ax,mul(u,dot(ax,u)))}return unit(q)}

export const xhLength=el=>XH[el]||1.09;
// Hydrogens an atom still needs, given the sum of its bond orders (Kekule) and formal charge. 0 for elements outside C, N, O, S.
export function missingHydrogens(el,orderSum,charge=0){const v=VALENCE[el];if(v===undefined)return 0;const target=el==='C'?v-Math.abs(charge):v+charge,m=Math.round(target-orderSum);return m>0?m:0}

// atoms: [{el,x,y,z}] heavy atoms; bonds: [[i,j]]; orders: Kekule order per bond.
// opts.charges: {atomIndex: formal charge}; opts.torsions: {atomIndex: degrees or [degrees]} for hydrogens on a terminal rotor (OH, SH, NH3, CH3).
// Returns {atoms (heavy atoms first, then hydrogens with .parent), bonds, orders, hydrogens}.
export function addHydrogens(atoms,bonds,orders,{charges={},torsions={}}={}){
 const n=atoms.length,el=i=>atoms[i].el||atoms[i].element,P=i=>[atoms[i].x,atoms[i].y,atoms[i].z],nb=atoms.map(()=>[]);
 bonds.forEach(([i,j],k)=>{const o=orders[k]||1;nb[i].push([j,o]);nb[j].push([i,o])});
 const sum=i=>nb[i].reduce((s,[,o])=>s+o,0),unsat=i=>nb[i].some(([,o])=>o>=2);
 const out=atoms.map(a=>({...a})),B=bonds.map(b=>[b[0],b[1]]),O=[...orders];let added=0;
 function directions(i,m){
  const e=el(i),u=nb[i].map(([j])=>unit(sub(P(j),P(i)))),d=u.length;
  const planar=unsat(i)||((e==='N'||e==='O')&&nb[i].some(([j])=>unsat(j))); // sp2: double bond here, or N/O conjugated with an unsaturated neighbour
  if(d===0)return[[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]].map(unit).slice(0,m);
  if(d===1){
   const a=u[0],j=nb[i][0][0];
   if(nb[i][0][1]===3)return[mul(a,-1)];                           // sp: H-C#N, H-C#C
   const rs=nb[j].filter(([k])=>k!==i),ref=rs.find(([,o])=>o>=2)||rs[0];  // reference atom: prefer one doubly bonded to the neighbour
   const e1=ref?perp(a,sub(P(ref[0]),P(j))):perp(a),e2=cross(a,e1),th=(e==='S'?96:e==='O'?109.5:planar?120:109.5)*R; // hydroxyl C-O-H ~109.5 even when conjugated (phenol, carboxylic acid); planar NH2 / =CH2 120
   let phis=planar?[0,180]:[180,60,-60];if(torsions[i]!==undefined)phis=[].concat(torsions[i]);
   return phis.slice(0,m).map(f=>add(mul(a,Math.cos(th)),mul(add(mul(e1,Math.cos(f*R)),mul(e2,Math.sin(f*R))),Math.sin(th))));
  }
  if(d===2){
   const [a,b]=u,s=add(a,b),bis=len(s)<1e-6?perp(a):mul(unit(s),-1),nrm=len(cross(a,b))<1e-6?perp(bis):unit(cross(a,b));
   if(m===1&&planar)return[bis];                                    // aromatic C-H, amide / imidazole N-H: in the plane
   const al=54.75*R,c=Math.cos(al),sn=Math.sin(al);
   return[add(mul(bis,c),mul(nrm,sn)),add(mul(bis,c),mul(nrm,-sn))].slice(0,m);
  }
  if(d===3&&m===1){const s=add(add(u[0],u[1]),u[2]);return[len(s)<1e-6?unit(cross(u[0],u[1])):mul(unit(s),-1)]}
  return[];
 }
 for(let i=0;i<n;i++){
  const e=el(i),m=missingHydrogens(e,sum(i),charges[i]||0);if(!m)continue;
  const L=xhLength(e),p=P(i);
  for(const dir of directions(i,m)){const q=add(p,mul(dir,L)),k=out.length;out.push({el:'H',x:q[0],y:q[1],z:q[2],parent:i});B.push([i,k]);O.push(1);added++}
 }
 return{atoms:out,bonds:B,orders:O,hydrogens:added};
}

// Hydrogens for the idealised peptide models (atoms N, CA, C, O and CB per residue; chain + residue number group the atoms). Returns a new array:
// the original heavy atoms unchanged and in place, then one OXT per chain, then the hydrogens.
// Chemistry: a free peptide with neutral ends. N-terminus -NH2, C-terminus -COOH (OXT added), amide N-H, C-alpha H (two for Gly), CB as CH3.
// Side chains beyond CB are not built, so CB is capped like alanine; Pro has no ring in these models, so its N carries an H like the other residues.
export function addPeptideHydrogens(atoms){
 const heavy=atoms.map(a=>({...a})),P=a=>[a.x,a.y,a.z],key=a=>a.chain+'|'+a.resi,byRes=new Map(),chains=new Map(),bonds=[],ords=[];
 heavy.forEach((a,i)=>{const k=key(a);if(!byRes.has(k)){byRes.set(k,{});if(!chains.has(a.chain))chains.set(a.chain,[]);chains.get(a.chain).push(byRes.get(k))}byRes.get(k)[a.name]=i});
 const bond=(i,j,o)=>{bonds.push([i,j]);ords.push(o)};
 for(const rs of chains.values())rs.forEach((r,n)=>{bond(r.N,r.CA,1);bond(r.CA,r.C,1);bond(r.C,r.O,2);if(r.CB!==undefined)bond(r.CA,r.CB,1);if(n)bond(rs[n-1].C,r.N,1)});
 for(const rs of chains.values()){const r=rs[rs.length-1],C=heavy[r.C],a=unit(sub(P(heavy[r.CA]),P(C))),b=unit(sub(P(heavy[r.O]),P(C))),p=add(P(C),mul(unit(add(a,b)),-1.3)),k=heavy.length; // trigonal: the free C-terminal O lies opposite CA and O
  heavy.push({el:'O',x:p[0],y:p[1],z:p[2],name:'OXT',res:C.res,resi:C.resi,chain:C.chain,serial:k+1});bond(r.C,k,1);r.OXT=k}
 const r=addHydrogens(heavy,bonds,ords),total=new Map(),seen=new Map(),firstN=new Set([...chains.values()].map(rs=>rs[0].N));
 for(let k=heavy.length;k<r.atoms.length;k++)total.set(r.atoms[k].parent,(total.get(r.atoms[k].parent)||0)+1);
 return r.atoms.map((a,k)=>{if(k<heavy.length)return a;const p=heavy[a.parent],n=(seen.get(a.parent)||0)+1,t=total.get(a.parent);seen.set(a.parent,n);
  const name=p.name==='N'?(firstN.has(a.parent)?'H'+n:'H'):p.name==='CA'?(t>1?'HA'+(n+1):'HA'):p.name==='CB'?'HB'+n:'HXT'; // N-terminal H1 H2, Gly HA2 HA3, methyl HB1-3, acid HXT
  return{el:'H',x:a.x,y:a.y,z:a.z,name,res:p.res,resi:p.resi,chain:p.chain,serial:k+1}})}
