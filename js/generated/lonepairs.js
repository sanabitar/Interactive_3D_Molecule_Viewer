// GENERATED FILE - do not edit. Source: src/core/lonepairs.mjs
// Rebuild with: node tools/build_viewer.mjs
const LonePairs=(()=>{
// Lone electron pairs: how many each atom has (from valence) and in which directions to draw them (from the real geometry around the atom).
// Pure functions, no DOM. A lone pair is NOT an atom: this module only returns counts and unit directions, never atoms.
//
// Count:  lone pairs = floor((valence electrons - formal charge - sum of bond orders) / 2), for N, O, S, P, Se and the halogens.
//         It is only counted when the atom's bonding is COMPLETE (bond-order sum >= that of the neutral atom: N 3, O 2, S 2, halogen 1, ...).
//         An atom that looks short of bonds (usually because hydrogens are missing from the file) gets none: a missing N-H must not look like a lone pair.
//         Delocalised (1.5) bonds count 1.5, so an aromatic pyrrole-type N-H or furan O (lone pair in the pi system) correctly gets none.
// Place:  the directions complete the geometry that is really there around the atom, so they follow the electron-region geometry:
//         tetrahedral (sp3: water, ammonia, amines, alcohols, ethers), trigonal in the molecular plane (sp2: C=O oxygen, pyridine N),
//         perpendicular to the plane for a lone pair conjugated with a pi system (aniline / amide N, phenol O: the "p" lone pair), linear (sp: N2, HCN, CO).
const R=Math.PI/180;
const VALENCE_E={N:5,O:6,S:6,P:5,Se:6,F:7,Cl:7,Br:7,I:7};      // valence electrons
const NORMAL_SUM={N:3,O:2,S:2,P:3,Se:2,F:1,Cl:1,Br:1,I:1};     // bond-order sum of the neutral atom
const add=(a,b)=>[a[0]+b[0],a[1]+b[1],a[2]+b[2]],sub=(a,b)=>[a[0]-b[0],a[1]-b[1],a[2]-b[2]],mul=(a,s)=>[a[0]*s,a[1]*s,a[2]*s];
const dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2],cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],len=a=>Math.hypot(a[0],a[1],a[2]);
const unit=a=>{const l=len(a);return l>1e-12?mul(a,1/l):[0,0,0]};
function perp(u,hint){let q=hint?sub(hint,mul(u,dot(hint,u))):null;if(!q||len(q)<1e-6){const ax=Math.abs(u[0])<.9?[1,0,0]:[0,1,0];q=sub(ax,mul(u,dot(ax,u)))}return unit(q)}
const TETRA=[[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]].map(unit);

// Lone pairs of one atom from its element, bond-order sum and formal charge. 0 for other elements and for atoms whose bonding looks incomplete.
function lonePairCount(el,orderSum,charge=0){
 const v=VALENCE_E[el];if(v===undefined)return 0;
 if(orderSum<NORMAL_SUM[el]+charge-1e-9)return 0;
 return Math.max(0,Math.floor((v-charge-orderSum)/2+1e-9))}

// atoms: [{el|element,x,y,z}]; bonds: [[i,j]]; orders: bond order per bond (1, 2, 3, 1.5 delocalised).
// opts.charges: Map(atomIndex -> formal charge) for the charges that can be read off the structure (N+ in NH3+, carboxylate O-).
// opts.hasH: does the structure contain hydrogens at all (default: look). Without hydrogens, delocalised atoms are skipped: a ring N could be pyridine-type or N-H.
// Returns [{atom, dir:[x,y,z] unit vector from the atom, type:'sp3'|'sp2'|'p'|'sp'}], one entry per lone pair.
function computeLonePairs(atoms,bonds,orders,opts={}){
 const el=i=>atoms[i].el||atoms[i].element,P=i=>[atoms[i].x,atoms[i].y,atoms[i].z],charges=opts.charges||new Map(),hasH=opts.hasH!==undefined?opts.hasH:atoms.some((_,i)=>el(i)==='H');
 const nb=atoms.map(()=>[]);bonds.forEach(([i,j],k)=>{const o=orders[k]||1;nb[i].push([j,o]);nb[j].push([i,o])});
 const sum=i=>nb[i].reduce((s,[,o])=>s+o,0),unsat=i=>nb[i].some(([,o])=>o>1);   // unsaturated: a double, triple or delocalised bond
 const out=[];
 for(let i=0;i<atoms.length;i++){
  const e=el(i),d=nb[i].length;let m;
  if(e==='C')m=d===1&&nb[i][0][1]===3&&(el(nb[i][0][0])==='O'||el(nb[i][0][0])==='N')?1:0;   // carbon monoxide / isocyanide carbon (formally C-)
  else m=lonePairCount(e,sum(i),charges.get(i)||0);
  if(!m)continue;
  if(!hasH&&nb[i].some(([,o])=>o%1!==0))continue;
  const p=P(i),u=nb[i].map(([j])=>unit(sub(P(j),p)));let dirs=[];   // dirs: [direction, type]
  const conj=(e==='N'||e==='O'||e==='S')&&nb[i].some(([j])=>nb[j].some(([k,o])=>k!==i&&o>1));   // lone pair next to a pi system
  const planar=unsat(i)||conj;
  if(d===0)dirs=TETRA.slice(0,m).map(t=>[t,'sp3']);
  else if(d===1){
   const a=u[0],j=nb[i][0][0],o=nb[i][0][1];
   if(o===3)dirs=[[mul(a,-1),'sp']];                                   // N2, HCN, CO
   else{
    const rs=nb[j].filter(([k])=>k!==i),ref=rs.find(([,q])=>q>1)||rs[0],e1=ref?perp(a,sub(P(ref[0]),P(j))):perp(a),e2=cross(a,e1);
    const at=(deg,ph)=>add(mul(a,Math.cos(deg*R)),mul(add(mul(e1,Math.cos(ph*R)),mul(e2,Math.sin(ph*R))),Math.sin(deg*R)));
    if(planar){dirs=[[at(120,0),'sp2'],[at(120,180),'sp2']].slice(0,m);if(m>=3)dirs.push([e2,'p'])}   // carbonyl O: two in the plane; carboxylate O-: two in the plane + one p
    else dirs=[180,60,-60].slice(0,m).map(ph=>[at(109.5,ph),'sp3']);   // halide, alkoxide-like: three staggered
   }}
  else if(d===2){
   const[a,b]=u,s=add(a,b),bis=len(s)<1e-6?perp(a):mul(unit(s),-1),nrm=len(cross(a,b))<1e-6?perp(bis):unit(cross(a,b)),al=54.75*R;
   const tet=[add(mul(bis,Math.cos(al)),mul(nrm,Math.sin(al))),add(mul(bis,Math.cos(al)),mul(nrm,-Math.sin(al)))];
   if(planar)dirs=m===1?[[bis,'sp2']]:[[bis,'sp2'],[nrm,'p']];             // pyridine N: in the plane; conjugated O (phenol, acid OH): one in the plane, one p
   else dirs=tet.slice(0,m).map(t=>[t,'sp3']);                              // water, alcohols, ethers, thiols
  }
  else if(d===3&&m===1){
   const s=add(add(u[0],u[1]),u[2]);
   dirs=len(s)<.3?[[len(cross(u[0],u[1]))<1e-6?perp(u[0]):unit(cross(u[0],u[1])),'p']]:[[mul(unit(s),-1),'sp3']];   // pyramidal (NH3, amines) or planar-conjugated (aniline, amide)
  }
  for(const[dir,type]of dirs)out.push({atom:i,dir,type});
 }
 return out}

return{lonePairCount,computeLonePairs}})();
