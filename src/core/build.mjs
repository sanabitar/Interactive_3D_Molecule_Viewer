// Constructive geometry: Z-matrix (NeRF) builder and idealised peptide builder. Angstrom / degrees. Output is idealised, not experimental.
const DEG=Math.PI/180,vsub=(a,b)=>[a[0]-b[0],a[1]-b[1],a[2]-b[2]],vadd=(a,b)=>[a[0]+b[0],a[1]+b[1],a[2]+b[2]],vmul=(a,s)=>[a[0]*s,a[1]*s,a[2]*s],vlen=a=>Math.hypot(...a),vnorm=a=>vmul(a,1/vlen(a)),
 vcross=(u,v)=>[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]];
// spec rows: [el] | [el,i,d] | [el,i,d,j,angle] | [el,i,d,j,angle,k,torsion]  (new atom bonded to i; angle(new,i,j); dihedral(new,i,j,k))
export function zmatrix(spec){const P=[];for(const[el,i,d,j,ang,k,tor]of spec){let p;
  if(i===undefined)p=[0,0,0];else if(j===undefined)p=vadd(P[i],[0,0,d]);
  else if(k===undefined){const u=vnorm(vsub(P[j],P[i]));let q=vsub([1,0,0],vmul(u,u[0]));if(vlen(q)<1e-8)q=vsub([0,1,0],vmul(u,u[1]));q=vnorm(q);p=vadd(P[i],vadd(vmul(u,d*Math.cos(ang*DEG)),vmul(q,d*Math.sin(ang*DEG))))}
  else{const a=P[k],b=P[j],c=P[i],bc=vnorm(vsub(c,b)),n=vnorm(vcross(vsub(b,a),bc)),m=vcross(n,bc),t=ang*DEG,w=tor*DEG;p=vadd(c,vadd(vmul(bc,-d*Math.cos(t)),vadd(vmul(m,d*Math.sin(t)*Math.cos(w)),vmul(n,d*Math.sin(t)*Math.sin(w)))))}
  P.push(p)}return P.map((p,n)=>({el:spec[n][0],x:p[0],y:p[1],z:p[2]}))}
export const CB_TORSION=122.6; // dihedral(N,C,CA,CB); sign chosen so that the residue is L (verified in tests by the CORN rule)
// Idealised peptide, heavy atoms only (N,CA,C,O and CB for every non-Gly residue); typical backbone geometry, trans peptide bonds by default.
// phi/psi/omega: one number for all residues, or an array with one entry per residue (phi[0] is unused; psi[n-1] only orients the last C=O).
// seq: optional array of 3-letter residue names (length n). Side chains beyond CB are NOT built; Pro has no ring.
export function peptide(n,{phi,psi,omega=180,resName='ALA',seq}={}){if(!(n>=2))throw new RangeError('At least 2 residues.');
 const names=seq?seq.map(x=>String(x).toUpperCase()):Array(n).fill(resName);if(names.length!==n)throw new RangeError('seq must have one residue name per residue (n = '+n+').');
 const chk=(v,k)=>{if(Array.isArray(v)?(v.length!==n||v.some(x=>!Number.isFinite(x))):!Number.isFinite(v))throw new RangeError(k+' must be a finite number or an array of n numbers.')};chk(phi,'phi');chk(psi,'psi');chk(omega,'omega');
 const at=(v,i)=>Array.isArray(v)?v[i]:v;const S=[],T=[],ix=[];
 for(let r=0;r<n;r++){const prev=ix[r-1],o={};
  if(!prev){o.N=S.push(['N'])-1;o.CA=S.push(['C',o.N,1.458])-1;o.C=S.push(['C',o.CA,1.525,o.N,111.2])-1}
  else{o.N=S.push(['N',prev.C,1.329,prev.CA,116.2,prev.N,at(psi,r-1)])-1;o.CA=S.push(['C',o.N,1.458,prev.C,121.7,prev.CA,at(omega,r-1)])-1;o.C=S.push(['C',o.CA,1.525,o.N,111.2,prev.C,at(phi,r)])-1}
  o.O=S.push(['O',o.C,1.231,o.CA,120.8,o.N,at(psi,r)+180])-1;if(names[r]!=='GLY')o.CB=S.push(['C',o.CA,1.53,o.C,110.1,o.N,CB_TORSION])-1;ix.push(o);
  for(const nm of['N','CA','C','O','CB'])if(o[nm]!==undefined)T[o[nm]]={name:nm,resi:r+1,res:names[r]}}
 return zmatrix(S).map((a,k)=>({el:a.el,x:a.x,y:a.y,z:a.z,name:T[k].name,res:T[k].res,resi:T[k].resi,chain:'A',serial:k+1}))}
// Rigid-body copy of an atom list: rotate (3x3 row-major R), translate by t, relabel chain, renumber serials from serialStart.
export function transformAtoms(atoms,{R=[[1,0,0],[0,1,0],[0,0,1]],t=[0,0,0],chain,serialStart=1}={}){return atoms.map((a,k)=>({...a,x:R[0][0]*a.x+R[0][1]*a.y+R[0][2]*a.z+t[0],y:R[1][0]*a.x+R[1][1]*a.y+R[1][2]*a.z+t[1],z:R[2][0]*a.x+R[2][1]*a.y+R[2][2]*a.z+t[2],chain:chain??a.chain,serial:serialStart+k}))}
