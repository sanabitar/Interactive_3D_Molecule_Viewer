// Descriptive properties from coordinates. Not a similarity measure. Angstrom, g/mol.
// Abridged conventional standard atomic weights (IUPAC), g/mol.
export const ATOMIC_MASS={H:1.008,He:4.0026,Li:6.94,Be:9.0122,B:10.81,C:12.011,N:14.007,O:15.999,F:18.998,Na:22.990,Mg:24.305,Al:26.982,Si:28.085,P:30.974,S:32.06,Cl:35.45,K:39.098,Ca:40.078,Mn:54.938,Fe:55.845,Co:58.933,Ni:58.693,Cu:63.546,Zn:65.38,Se:78.971,Br:79.904,I:126.904};
// Hill system: C first, then H, then alphabetical; if no carbon, everything alphabetical.
export function hillFormula(atoms){const c={};for(const a of atoms){const e=a.element||a.el;if(!e)continue;c[e]=(c[e]||0)+1}
 const ks=Object.keys(c).sort(),order=c.C?['C',...(c.H?['H']:[]),...ks.filter(k=>k!=='C'&&k!=='H')]:ks;return{counts:c,formula:order.map(k=>k+(c[k]>1?c[k]:'')).join('')}}
export function molProps(atoms){if(!atoms||!atoms.length)throw new RangeError('No atoms to describe.');const{counts,formula}=hillFormula(atoms),unknown=Object.keys(counts).filter(e=>!(e in ATOMIC_MASS)),complete=!unknown.length;
 const mw=complete?Object.entries(counts).reduce((s,[e,n])=>s+ATOMIC_MASS[e]*n,0):null,w=a=>complete?ATOMIC_MASS[a.element||a.el]:1,W=atoms.reduce((s,a)=>s+w(a),0);
 const cm=['x','y','z'].map(k=>atoms.reduce((s,a)=>s+w(a)*a[k],0)/W),rg=Math.sqrt(atoms.reduce((s,a)=>s+w(a)*((a.x-cm[0])**2+(a.y-cm[1])**2+(a.z-cm[2])**2),0)/W);
 const ext=k=>{let lo=Infinity,hi=-Infinity;for(const a of atoms){const v=a[k];if(v<lo)lo=v;if(v>hi)hi=v}return hi-lo},H=counts.H||0;
 return{atomCount:atoms.length,hydrogens:H,heavyAtoms:atoms.length-H,formula,counts,molecularWeight:mw,unknownElements:unknown,radiusOfGyration:rg,rgMassWeighted:complete,extent:{x:ext('x'),y:ext('y'),z:ext('z')},
  hydrogensLikelyMissing:H===0&&(counts.C||0)>0,residues:new Set(atoms.filter(a=>!a.hetero&&a.resName).map(a=>`${a.chain}:${a.resSeq}${a.iCode||''}`)).size,chains:new Set(atoms.filter(a=>a.chain).map(a=>a.chain)).size}}
