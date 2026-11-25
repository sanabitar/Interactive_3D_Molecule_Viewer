// js/foundation/core.js: Core: parsers, bond perception, analysis, models and sessions, plus parseFile() and computeBonds().
// A hand-maintained copy of the pure modules in src/core; tests/viewer_parity.test.mjs guards against drift.
// ---------- parsers & validation ----------
const Core=(()=>{
// Structure parsing + validation. Pure functions, no DOM. Units: Angstrom.
class ParseError extends Error{constructor(m,line){super(line?`Line ${line}: ${m}`:m);this.name='ParseError';this.line=line}}
const cap=e=>e?e[0].toUpperCase()+e.slice(1).toLowerCase():'';
const residueKey=a=>`${a.chain}:${a.resSeq}${a.iCode||''}`; // chain + number + insertion code
const num=(s,line,what)=>{if(s===undefined||String(s).trim()===''||!Number.isFinite(Number(s)))throw new ParseError(`invalid ${what} "${String(s).trim()}"`,line);return Number(s)};
function finish(atoms,warnings,format){if(!atoms.length)throw new ParseError('No ATOM/HETATM records found.');
 return{format,atoms,warnings,chains:[...new Set(atoms.map(a=>a.chain))],residues:new Set(atoms.map(a=>residueKey(a)+a.resName)).size}}
const TWO_LETTER=new Set(['ZN','FE','MG','CA','NA','CL','MN','CU','CO','NI','CD','HG','BR','SE','LI','AL','SI','BE','BA','SR','CS','RB','AU','AG','PT','PB','SN','MO','CR','TI']);
function parsePDB(text){const atoms=[],warnings=[],seen=new Set(),L=text.split(/\r?\n/);let alt=null,skippedAlt=0,inferred=0,models=0,done=false,dup=0;
 for(let n=0;n<L.length;n++){const l=L[n];if(l.startsWith('MODEL'))models++;if(l.startsWith('ENDMDL')){done=true;continue}if(done||!/^(ATOM  |HETATM)/.test(l))continue;
  if(l.length<54)throw new ParseError('Atom record is shorter than 54 columns, so its coordinates are incomplete (the file may be truncated).',n+1);
  const al=l[16]||' ';if(al!==' '){if(alt===null)alt=al;if(al!==alt){skippedAlt++;continue}}
  let el=l.slice(76,78).trim();if(!el){const f=l.slice(12,16);el=l[0]==='H'&&/^[A-Za-z]{2}/.test(f)&&TWO_LETTER.has(f.slice(0,2).toUpperCase())?f.slice(0,2):(f.match(/[A-Za-z]/)||['X'])[0];inferred++}
  const a={serial:num(l.slice(6,11),n+1,'serial number'),name:l.slice(12,16).trim(),altLoc:al.trim(),resName:l.slice(17,20).trim(),chain:(l[21]||' ').trim(),resSeq:num(l.slice(22,26),n+1,'residue number'),iCode:(l[26]||' ').trim(),
   x:num(l.slice(30,38),n+1,'x'),y:num(l.slice(38,46),n+1,'y'),z:num(l.slice(46,54),n+1,'z'),element:cap(el),hetero:l.startsWith('HETATM')};
  const k=`${a.chain}|${a.resSeq}|${a.iCode}|${a.name}|${a.altLoc}`;if(seen.has(k))dup++;else seen.add(k);atoms.push(a)}
 if(models>1)warnings.push(`${models-1} additional model(s) were ignored; only model 1 is read.`);
 if(skippedAlt)warnings.push(`${skippedAlt} atoms in alternate conformation(s) other than "${alt}" were skipped.`);
 if(inferred)warnings.push(`Element column missing for ${inferred} atoms; elements were inferred from atom names (heuristic).`);
 if(dup)warnings.push(`${dup} duplicate atom record(s) (same chain, residue, atom name and alternate location) were kept; check the source file.`);return finish(atoms,warnings,'pdb')}
const tok=l=>l.match(/'[^']*'|"[^"]*"|\S+/g)||[];
function parseMmCIF(text){const L=text.split(/\r?\n/),atoms=[],warnings=[];let i=0;
 while(i<L.length&&!(L[i].startsWith('_atom_site.')))i++;if(i>=L.length)throw new ParseError('No _atom_site loop found (not a coordinate mmCIF).');
 const H=[];while(i<L.length&&L[i].startsWith('_atom_site.'))H.push(L[i++].trim().slice(11));
 const c=k=>H.indexOf(k);for(const k of['Cartn_x','Cartn_y','Cartn_z'])if(c(k)<0)throw new ParseError('Missing _atom_site.'+k);
 const idx=(a,b)=>c(a)>=0?c(a):c(b);const ch=idx('auth_asym_id','label_asym_id'),sq=idx('auth_seq_id','label_seq_id'),rn=idx('auth_comp_id','label_comp_id'),an=idx('auth_atom_id','label_atom_id');
 let alt=null,skipped=0,model=null,buf=[];const models=new Set();
 for(;i<L.length;i++){const l=L[i];if(!l.trim())continue;if(/^(#|loop_|_)/.test(l))break;buf=buf.concat(tok(l));if(buf.length<H.length)continue;if(buf.length>H.length)throw new ParseError('Row has more fields than the loop header (check quoting).',i+1);const t=buf;buf=[];
  const m=c('pdbx_PDB_model_num')>=0?t[c('pdbx_PDB_model_num')]:'1';models.add(m);if(model===null)model=m;if(m!==model)continue;
  const al=c('label_alt_id')>=0?t[c('label_alt_id')]:'.';if(al!=='.'&&al!=='?'){if(alt===null)alt=al;if(al!==alt){skipped++;continue}}
  const ic=c('pdbx_PDB_ins_code')>=0?t[c('pdbx_PDB_ins_code')]:'?';
  atoms.push({serial:c('id')>=0?Number(t[c('id')]):atoms.length+1,name:t[an].replace(/^"|"$/g,''),altLoc:al==='.'||al==='?'?'':al,resName:t[rn],chain:t[ch],resSeq:num(t[sq],i+1,'residue number'),
   iCode:ic==='?'||ic==='.'?'':ic,x:num(t[c('Cartn_x')],i+1,'x'),y:num(t[c('Cartn_y')],i+1,'y'),z:num(t[c('Cartn_z')],i+1,'z'),element:cap(t[idx('type_symbol','type_symbol')]),hetero:c('group_PDB')>=0?t[c('group_PDB')]==='HETATM':false})}
 if(buf.length)throw new ParseError('Truncated mmCIF: the last atom_site row is incomplete.');
 if(models.size>1)warnings.push(`${models.size-1} additional model(s) were ignored; only model ${model} is read.`);
 if(skipped)warnings.push(`${skipped} atoms in alternate conformation(s) other than "${alt}" were skipped.`);
 warnings.push('mmCIF: author (auth_*) chain/residue identifiers used where present; only the first model is read.');return finish(atoms,warnings,'mmcif')}
function parseStructure(name,text,{maxBytes=60e6}={}){if(text.length>maxBytes)throw new ParseError('File exceeds the size limit.');
 if(/[\x00-\x08]/.test(text.slice(0,2000)))throw new ParseError('Binary data detected; expected a text structure file.');
 const ext=name.split('.').pop().toLowerCase();if(ext==='pdb'||ext==='ent')return parsePDB(text);if(ext==='cif'||ext==='mmcif')return parseMmCIF(text);
 throw new ParseError(`Unsupported format ".${ext}" (supported here: pdb, ent, cif, mmcif).`)}

// Geometry. Inputs {x,y,z} in Angstrom; angles returned in degrees.
const sub=(a,b)=>[b.x-a.x,b.y-a.y,b.z-a.z],dot=(u,v)=>u[0]*v[0]+u[1]*v[1]+u[2]*v[2],cross=(u,v)=>[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]],len=u=>Math.hypot(...u);
const distance=(a,b)=>len(sub(a,b));
function angle(a,b,c){const u=sub(b,a),v=sub(b,c),m=len(u)*len(v);if(m===0)throw new RangeError('Angle undefined: coincident atoms.');return Math.acos(Math.max(-1,Math.min(1,dot(u,v)/m)))*180/Math.PI}
// IUPAC sign convention; range (-180, 180]
function dihedral(a,b,c,d){const b1=sub(a,b),b2=sub(b,c),b3=sub(c,d);if(!len(b2)||!len(cross(b1,b2))||!len(cross(b2,b3)))throw new RangeError('Dihedral undefined: collinear or coincident atoms.');
 return Math.atan2(len(b2)*dot(b1,cross(b2,b3)),dot(cross(b1,b2),cross(b2,b3)))*180/Math.PI}

// Method B (ligand-associated site): purely geometric. Contact != catalytic. See caveat in output.

function ligandContacts(atoms,{ligandKey,cutoff=4.0,annotated=new Map()}={}){
 if(!(cutoff>0))throw new RangeError('cutoff must be > 0 Angstrom');
 const lig=atoms.filter(a=>a.hetero&&residueKey(a)===ligandKey);if(!lig.length)throw new Error(`Ligand ${ligandKey} not found among HETATM records.`);
 const prot=atoms.filter(a=>!a.hetero),best=new Map();
 for(const p of prot)for(const q of lig){const d=distance(p,q);if(d<=cutoff){const k=residueKey(p),b=best.get(k);if(!b||d<b.minDistance)best.set(k,{residue:k,resName:p.resName,minDistance:d,proteinAtom:p.name,ligandAtom:q.name})}}
 const contacts=[...best.values()].sort((a,b)=>a.minDistance-b.minDistance).map(c=>({...c,evidence:annotated.has(c.residue)?['geometric_contact','annotated_catalytic']:['geometric_contact'],annotationSource:annotated.get(c.residue)||null}));
 const present=new Set(prot.map(residueKey));
 return{method:'Heavy/any-atom minimum distance, ligand vs polymer ATOM records',params:{cutoffAngstrom:cutoff,ligand:ligandKey,ligandResName:lig[0].resName,altConfsHandled:'as parsed'},contacts,
  annotatedNotInContact:[...annotated.keys()].filter(k=>!best.has(k)).map(k=>({residue:k,presentInStructure:present.has(k),annotationSource:annotated.get(k)})),
  caveat:'Contact residues are geometric neighbours of a bound ligand. They are not established catalytic residues, and the ligand may not be a substrate or biologically relevant.'}}

// Least-squares rigid superposition (Horn 1987 quaternion method == Kabsch optimum). Equal weights, rotation+translation only (no reflection).

function eig4(N){const A=N.map(r=>[...r]),V=[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]];
 for(let sw=0;sw<80;sw++){let off=0;for(let i=0;i<4;i++)for(let j=i+1;j<4;j++)off+=A[i][j]**2;if(off<1e-26)break;
  for(let p=0;p<3;p++)for(let q=p+1;q<4;q++){if(Math.abs(A[p][q])<1e-300)continue;const th=(A[q][q]-A[p][p])/(2*A[p][q]),t=(th>=0?1:-1)/(Math.abs(th)+Math.sqrt(th*th+1)),c=1/Math.sqrt(t*t+1),s=t*c;
   for(let k=0;k<4;k++){const x=A[k][p],y=A[k][q];A[k][p]=c*x-s*y;A[k][q]=s*x+c*y}
   for(let k=0;k<4;k++){const x=A[p][k],y=A[q][k];A[p][k]=c*x-s*y;A[q][k]=s*x+c*y}
   for(let k=0;k<4;k++){const x=V[k][p],y=V[k][q];V[k][p]=c*x-s*y;V[k][q]=s*x+c*y}}}
 return{values:[0,1,2,3].map(i=>A[i][i]),vectors:V}}
// Finds the rotation+translation that moves P (moving) onto Q (reference). P[i] pairs with Q[i].
function superpose(P,Q){const n=P.length;if(n!==Q.length)throw new RangeError('Coordinate sets differ in length.');if(n<3)throw new RangeError('At least 3 coordinate pairs are required.');
 const cen=X=>({x:X.reduce((s,a)=>s+a.x,0)/n,y:X.reduce((s,a)=>s+a.y,0)/n,z:X.reduce((s,a)=>s+a.z,0)/n}),cp=cen(P),cq=cen(Q);
 let Sxx=0,Sxy=0,Sxz=0,Syx=0,Syy=0,Syz=0,Szx=0,Szy=0,Szz=0,G=0;
 for(let i=0;i<n;i++){const a=[P[i].x-cp.x,P[i].y-cp.y,P[i].z-cp.z],b=[Q[i].x-cq.x,Q[i].y-cq.y,Q[i].z-cq.z];
  Sxx+=a[0]*b[0];Sxy+=a[0]*b[1];Sxz+=a[0]*b[2];Syx+=a[1]*b[0];Syy+=a[1]*b[1];Syz+=a[1]*b[2];Szx+=a[2]*b[0];Szy+=a[2]*b[1];Szz+=a[2]*b[2];G+=a[0]**2+a[1]**2+a[2]**2+b[0]**2+b[1]**2+b[2]**2}
 const N=[[Sxx+Syy+Szz,Syz-Szy,Szx-Sxz,Sxy-Syx],[Syz-Szy,Sxx-Syy-Szz,Sxy+Syx,Szx+Sxz],[Szx-Sxz,Sxy+Syx,-Sxx+Syy-Szz,Syz+Szy],[Sxy-Syx,Szx+Sxz,Syz+Szy,-Sxx-Syy+Szz]];
 const{values,vectors}=eig4(N);let m=0;for(let k=1;k<4;k++)if(values[k]>values[m])m=k;
 let[w,x,y,z]=[0,1,2,3].map(r=>vectors[r][m]);const L=Math.hypot(w,x,y,z);w/=L;x/=L;y/=L;z/=L;
 const R=[[w*w+x*x-y*y-z*z,2*(x*y-w*z),2*(x*z+w*y)],[2*(x*y+w*z),w*w-x*x+y*y-z*z,2*(y*z-w*x)],[2*(x*z-w*y),2*(y*z+w*x),w*w-x*x-y*y+z*z]];
 const apply=p=>{const a=[p.x-cp.x,p.y-cp.y,p.z-cp.z];return{x:R[0][0]*a[0]+R[0][1]*a[1]+R[0][2]*a[2]+cq.x,y:R[1][0]*a[0]+R[1][1]*a[1]+R[1][2]*a[2]+cq.y,z:R[2][0]*a[0]+R[2][1]*a[1]+R[2][2]*a[2]+cq.z}};
 const deviations=P.map((p,i)=>distance(apply(p),Q[i]));
 return{rmsd:Math.sqrt(deviations.reduce((s,d)=>s+d*d,0)/n),rmsdFromEigenvalue:Math.sqrt(Math.max(0,(G-2*values[m])/n)),rotation:R,apply,deviations,n}}

// Sequence alignment, phi/psi, CA contact map. Pure functions. Angles in degrees, distances in Angstrom.

const AA={ALA:'A',ARG:'R',ASN:'N',ASP:'D',CYS:'C',GLN:'Q',GLU:'E',GLY:'G',HIS:'H',ILE:'I',LEU:'L',LYS:'K',MET:'M',PHE:'F',PRO:'P',SER:'S',THR:'T',TRP:'W',TYR:'Y',VAL:'V',MSE:'M'};
function residues(atoms,chain){const m=new Map();for(const a of atoms){if(a.hetero&&a.resName!=='MSE')continue;if(chain!==undefined&&a.chain!==chain)continue;const k=residueKey(a);let r=m.get(k);if(!r){r={key:k,resName:a.resName,chain:a.chain,atoms:{}};m.set(k,r)}if(!(a.name in r.atoms))r.atoms[a.name]=a}return[...m.values()]}
const chainSequence=(atoms,chain)=>residues(atoms,chain).map(r=>AA[r.resName]||'X').join('');
// Needleman-Wunsch global alignment, linear gap penalty, simple match/mismatch scoring (no substitution matrix).
function alignGlobal(a,b,{match=2,mismatch=-1,gap=-2}={}){a=String(a).replace(/\s/g,'').toUpperCase();b=String(b).replace(/\s/g,'').toUpperCase();
 if(!a||!b)throw new RangeError('Both sequences must be non-empty.');if(!/^[A-Z]+$/.test(a+b))throw new RangeError('Sequences may contain letters only.');
 const n=a.length,m=b.length,W=m+1;if(n*m>4e6)throw new RangeError('Sequences too long for in-browser alignment (limit n*m = 4,000,000).');
 const S=new Float64Array((n+1)*W);for(let i=1;i<=n;i++)S[i*W]=i*gap;for(let j=1;j<=m;j++)S[j]=j*gap;
 for(let i=1;i<=n;i++)for(let j=1;j<=m;j++)S[i*W+j]=Math.max(S[(i-1)*W+j-1]+(a[i-1]===b[j-1]?match:mismatch),S[(i-1)*W+j]+gap,S[i*W+j-1]+gap);
 let i=n,j=m,x='',y='';while(i>0||j>0){if(i>0&&j>0&&S[i*W+j]===S[(i-1)*W+j-1]+(a[i-1]===b[j-1]?match:mismatch)){x=a[--i]+x;y=b[--j]+y}else if(i>0&&S[i*W+j]===S[(i-1)*W+j]+gap){x=a[--i]+x;y='-'+y}else{x='-'+x;y=b[--j]+y}}
 let identical=0,pairs=0;for(let k=0;k<x.length;k++)if(x[k]!=='-'&&y[k]!=='-'){pairs++;if(x[k]===y[k])identical++}
 return{a:x,b:y,score:S[n*W+m],length:x.length,identical,alignedPairs:pairs,gaps:x.length-pairs,identityPct:100*identical/x.length,coverageA:pairs/n,coverageB:pairs/m,params:{algorithm:'Needleman-Wunsch (global)',match,mismatch,gap,matrix:'none'}}}
// Schematic reference regions for a Ramachandran plot (degrees): generous boxes [phiMin,phiMax,psiMin,psiMax] around the commonly populated areas.
// They are NOT validation contours (those are density-derived from large reference sets) and the true allowed area is wider (e.g. the bridge region psi ~ 20-90 deg).
const RAMA_REGIONS=[{id:'alphaR',label:'α-helix (right-handed)',boxes:[[-150,-35,-75,20]]},{id:'beta',label:'β-strand / PPII',boxes:[[-180,-50,90,180],[-180,-60,-180,-165]]},{id:'alphaL',label:'left-handed helix',boxes:[[35,90,-10,80]]}];
const wrap180=v=>v>180-1e-9?v-360:v; // dihedral() returns (-180,180]; +180 and -180 are the same angle
// Glycine is not classified ('gly'): it lacks a side chain, so its allowed area is much larger and roughly symmetric. Anything else outside the boxes is 'other' (loops, turns, bridge region, or strained).
function ramaRegion(phi,psi,resName){if(resName==='GLY')return'gly';phi=wrap180(phi);psi=wrap180(psi);for(const r of RAMA_REGIONS)for(const[a,b,c,d]of r.boxes)if(phi>=a&&phi<=b&&psi>=c&&psi<=d)return r.id;return'other'}
function ramaSummary(points){const counts={alphaR:0,beta:0,alphaL:0,other:0,gly:0},outside=[];for(const p of points){counts[p.region]++;if(p.region==='other')outside.push(p)}return{total:points.length,counts,outside}}
function phiPsi(atoms,chain){const R=residues(atoms,chain),out=[];for(let i=0;i<R.length;i++){const r=R[i],p=R[i-1],n=R[i+1],a=r.atoms;if(!a.N||!a.CA||!a.C)continue;let phi=null,psi=null;
 if(p&&p.chain===r.chain&&p.atoms.C&&distance(p.atoms.C,a.N)<2.0)phi=dihedral(p.atoms.C,a.N,a.CA,a.C);
 if(n&&n.chain===r.chain&&n.atoms.N&&distance(a.C,n.atoms.N)<2.0)psi=dihedral(a.N,a.CA,a.C,n.atoms.N);
 if(phi!==null&&psi!==null)out.push({residue:r.key,resName:r.resName,phi,psi,region:ramaRegion(phi,psi,r.resName)})}return out}
function caContactMap(atoms,chain,cutoff=8){if(!(cutoff>0))throw new RangeError('cutoff must be > 0 Angstrom');const R=residues(atoms,chain).filter(r=>r.atoms.CA),n=R.length;
 if(!n)throw new RangeError('No CA atoms found for this chain.');if(n>2000)throw new RangeError('Too many residues (limit 2000).');const M=new Uint8Array(n*n);
 for(let i=0;i<n;i++)for(let j=i;j<n;j++)if(distance(R[i].atoms.CA,R[j].atoms.CA)<=cutoff)M[i*n+j]=M[j*n+i]=1;return{keys:R.map(r=>r.key),n,matrix:M,cutoff}}
// Pair atoms between structures A (reference) and B (moved), superpose B onto A, report RMSD over the paired atoms.
function compareStructures(A,B,{mode='ca',chainA,chainB,identicalOnly=false}={}){let pa,pb,aln=null,method;
 if(mode==='ca'){const ra=residues(A,chainA).filter(r=>r.atoms.CA),rb=residues(B,chainB).filter(r=>r.atoms.CA);if(ra.length<3||rb.length<3)throw new RangeError('Each chain needs at least 3 residues with CA atoms.');
  const sq=R=>R.map(r=>AA[r.resName]||'X').join('');aln=alignGlobal(sq(ra),sq(rb));pa=[];pb=[];let i=0,j=0;
  for(let k=0;k<aln.length;k++){const x=aln.a[k],y=aln.b[k];if(x!=='-'&&y!=='-'&&(!identicalOnly||x===y)){pa.push(ra[i].atoms.CA);pb.push(rb[j].atoms.CA)}if(x!=='-')i++;if(y!=='-')j++}
  method='CA atoms paired by sequence alignment'+(identicalOnly?' (identical residues only)':'')}
 else{pa=A.filter(a=>a.element!=='H');pb=B.filter(a=>a.element!=='H');if(pa.length!==pb.length||pa.some((a,i)=>a.element!==pb[i].element))throw new RangeError('Atom-order pairing needs identical heavy-atom element sequences; otherwise the comparison is not valid.');method='Heavy atoms paired by file order'}
 const s=superpose(pb,pa),key=a=>residueKey(a)+':'+a.name;
 return{method,params:{mode,chainA,chainB,identicalOnly,fit:'Horn/Kabsch least squares, rotation+translation, equal weights, no outlier rejection'},nPairs:pa.length,rmsd:s.rmsd,maxDeviation:s.deviations.reduce((m,d)=>d>m?d:m,0),rotation:s.rotation,
  alignment:aln&&{length:aln.length,identityPct:aln.identityPct,coverageA:aln.coverageA,coverageB:aln.coverageB},pairs:pa.map((a,i)=>({a:key(a),b:key(pb[i]),d:s.deviations[i]})),pa,pb,moved:pb.map(s.apply)}}

// Constructive geometry: Z-matrix (NeRF) builder and idealised peptide builder. Angstrom / degrees. Output is idealised, not experimental.
const DEG=Math.PI/180,vsub=(a,b)=>[a[0]-b[0],a[1]-b[1],a[2]-b[2]],vadd=(a,b)=>[a[0]+b[0],a[1]+b[1],a[2]+b[2]],vmul=(a,s)=>[a[0]*s,a[1]*s,a[2]*s],vlen=a=>Math.hypot(...a),vnorm=a=>vmul(a,1/vlen(a)),
 vcross=(u,v)=>[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]];
// spec rows: [el] | [el,i,d] | [el,i,d,j,angle] | [el,i,d,j,angle,k,torsion]  (new atom bonded to i; angle(new,i,j); dihedral(new,i,j,k))
function zmatrix(spec){const P=[];for(const[el,i,d,j,ang,k,tor]of spec){let p;
  if(i===undefined)p=[0,0,0];else if(j===undefined)p=vadd(P[i],[0,0,d]);
  else if(k===undefined){const u=vnorm(vsub(P[j],P[i]));let q=vsub([1,0,0],vmul(u,u[0]));if(vlen(q)<1e-8)q=vsub([0,1,0],vmul(u,u[1]));q=vnorm(q);p=vadd(P[i],vadd(vmul(u,d*Math.cos(ang*DEG)),vmul(q,d*Math.sin(ang*DEG))))}
  else{const a=P[k],b=P[j],c=P[i],bc=vnorm(vsub(c,b)),n=vnorm(vcross(vsub(b,a),bc)),m=vcross(n,bc),t=ang*DEG,w=tor*DEG;p=vadd(c,vadd(vmul(bc,-d*Math.cos(t)),vadd(vmul(m,d*Math.sin(t)*Math.cos(w)),vmul(n,d*Math.sin(t)*Math.sin(w)))))}
  P.push(p)}return P.map((p,n)=>({el:spec[n][0],x:p[0],y:p[1],z:p[2]}))}
const CB_TORSION=122.6; // dihedral(N,C,CA,CB); sign chosen so that the residue is L (verified in tests by the CORN rule)
// Idealised peptide, heavy atoms only (N,CA,C,O and CB for every non-Gly residue); typical backbone geometry, trans peptide bonds by default.
// phi/psi/omega: one number for all residues, or an array with one entry per residue (phi[0] is unused; psi[n-1] only orients the last C=O).
// seq: optional array of 3-letter residue names (length n). Side chains beyond CB are NOT built; Pro has no ring.
function peptide(n,{phi,psi,omega=180,resName='ALA',seq}={}){if(!(n>=2))throw new RangeError('At least 2 residues.');
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
function transformAtoms(atoms,{R=[[1,0,0],[0,1,0],[0,0,1]],t=[0,0,0],chain,serialStart=1}={}){return atoms.map((a,k)=>({...a,x:R[0][0]*a.x+R[0][1]*a.y+R[0][2]*a.z+t[0],y:R[1][0]*a.x+R[1][1]*a.y+R[1][2]*a.z+t[1],z:R[2][0]*a.x+R[2][1]*a.y+R[2][2]*a.z+t[2],chain:chain??a.chain,serial:serialStart+k}))}

// Idealised multi-residue test models for the analysis tools. Constructed, NOT experimental structures.
// Heavy atoms only: backbone N, CA, C, O plus CB for every non-Gly residue. Side chains beyond CB are omitted and Pro has no ring.

const A3={A:'ALA',R:'ARG',N:'ASN',D:'ASP',C:'CYS',Q:'GLN',E:'GLU',G:'GLY',H:'HIS',I:'ILE',L:'LEU',K:'LYS',M:'MET',F:'PHE',P:'PRO',S:'SER',T:'THR',W:'TRP',Y:'TYR',V:'VAL'};
const seq3=s=>[...s].map(c=>{if(!A3[c])throw new RangeError('Unknown one-letter residue code "'+c+'".');return A3[c]});
const MODEL_SEQ={helixWT:'AEAAAKEAAAKA',helixVar:'GSAEAAAREAAGKA',hlh:'SEELLKKALEEGPDKEELLKKALEE',hairpin:'SWTWENGKWTWK',dimerA:'AEAAAKEAAAKA',dimerB:'AEAAAREAAGKA',ramaTest:'APAAAAAPAAGGAAAA'};
// ideal right-handed alpha-helix, any sequence
const alphaHelix=seq=>peptide(seq.length,{phi:-57,psi:-47,seq:seq3(seq)});
// helix (res 1-11) - 4-residue loop G-P-D-K (res 12-15) - helix (res 16-25); loop torsions found by numerical search so the helices pack antiparallel without clashes
function helixLoopHelix(){const n=25,phi=new Array(n).fill(-57),psi=new Array(n).fill(-47);[[79,16],[-63,150],[-59,-43],[-98,7]].forEach(([f,p],k)=>{phi[11+k]=f;psi[11+k]=p});return peptide(n,{phi,psi,seq:seq3(MODEL_SEQ.hlh)})}
// 12-residue beta-hairpin: strands 1-5 and 8-12, type I' turn at Asn6-Gly7. Torsions optimised for cross-strand CA-CA registry and N-H...O=C distances.
function betaHairpin(){const phi=[-130,-126,-127,-120,-148,66,95,-86,-113,-127,-110,-124],psi=[155,139,131,127,135,16,-8,135,127,131,132,117];return peptide(12,{phi,psi,seq:seq3(MODEL_SEQ.hairpin)})}
// 16-residue peptide whose 14 phi/psi pairs visit every region of the plot, including Gly-only and strained positions
function ramaTestPeptide(){const phi=[0,-65,-75,-57,60,57,-120,-65,-135,-60,80,75,-110,-90,-63,-120],psi=[130,145,145,-47,-120,45,130,-35,140,-45,-170,25,135,60,-41,130];return peptide(16,{phi,psi,seq:seq3(MODEL_SEQ.ramaTest)})}
// two helical chains A and B: B is a rigid antiparallel copy of an alpha-helix, axis `separation` A away, with its own sequence
const vs=(a,b)=>[a[0]-b[0],a[1]-b[1],a[2]-b[2]],vc=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],vn=a=>{const l=Math.hypot(...a);return a.map(x=>x/l)};
function antiparallelHelixPair({seqA=MODEL_SEQ.dimerA,seqB=MODEL_SEQ.dimerB,separation=10}={}){if(seqA.length!==seqB.length)throw new RangeError('Chains must have the same length (same backbone).');
 const A=alphaHelix(seqA),B0=alphaHelix(seqB),ca=A.filter(a=>a.name==='CA').map(a=>[a.x,a.y,a.z]),c=[];
 for(let i=0;i+3<ca.length;i++)c.push([0,1,2].map(k=>(ca[i][k]+ca[i+1][k]+ca[i+2][k]+ca[i+3][k])/4)); // 4-CA centroids cancel the helical periodicity
 const dir=vn(vs(c[c.length-1],c[0])),m=[0,1,2].map(k=>c.reduce((s,x)=>s+x[k],0)/c.length),e=[[1,0,0],[0,1,0],[0,0,1]].sort((p,q)=>Math.abs(p[0]*dir[0]+p[1]*dir[1]+p[2]*dir[2])-Math.abs(q[0]*dir[0]+q[1]*dir[1]+q[2]*dir[2]))[0];
 const w=vn(vc(dir,e)),u=vc(dir,w),R=[0,1,2].map(i=>[0,1,2].map(j=>2*u[i]*u[j]-(i===j?1:0))), // pi rotation about u (perpendicular to the axis) flips the axis direction
  Rm=R.map(r=>r[0]*m[0]+r[1]*m[1]+r[2]*m[2]),t=[0,1,2].map(i=>m[i]-Rm[i]+separation*w[i]);
 return[...A,...transformAtoms(B0,{R,t,chain:'B',serialStart:A.length+1})]}

// Descriptive properties from coordinates. Not a similarity measure. Angstrom, g/mol.
// Abridged conventional standard atomic weights (IUPAC), g/mol.
const ATOMIC_MASS={H:1.008,He:4.0026,Li:6.94,Be:9.0122,B:10.81,C:12.011,N:14.007,O:15.999,F:18.998,Na:22.990,Mg:24.305,Al:26.982,Si:28.085,P:30.974,S:32.06,Cl:35.45,K:39.098,Ca:40.078,Mn:54.938,Fe:55.845,Co:58.933,Ni:58.693,Cu:63.546,Zn:65.38,Se:78.971,Br:79.904,I:126.904};
// Hill system: C first, then H, then alphabetical; if no carbon, everything alphabetical.
function hillFormula(atoms){const c={};for(const a of atoms){const e=a.element||a.el;if(!e)continue;c[e]=(c[e]||0)+1}
 const ks=Object.keys(c).sort(),order=c.C?['C',...(c.H?['H']:[]),...ks.filter(k=>k!=='C'&&k!=='H')]:ks;return{counts:c,formula:order.map(k=>k+(c[k]>1?c[k]:'')).join('')}}
function molProps(atoms){if(!atoms||!atoms.length)throw new RangeError('No atoms to describe.');const{counts,formula}=hillFormula(atoms),unknown=Object.keys(counts).filter(e=>!(e in ATOMIC_MASS)),complete=!unknown.length;
 const mw=complete?Object.entries(counts).reduce((s,[e,n])=>s+ATOMIC_MASS[e]*n,0):null,w=a=>complete?ATOMIC_MASS[a.element||a.el]:1,W=atoms.reduce((s,a)=>s+w(a),0);
 const cm=['x','y','z'].map(k=>atoms.reduce((s,a)=>s+w(a)*a[k],0)/W),rg=Math.sqrt(atoms.reduce((s,a)=>s+w(a)*((a.x-cm[0])**2+(a.y-cm[1])**2+(a.z-cm[2])**2),0)/W);
 const ext=k=>{let lo=Infinity,hi=-Infinity;for(const a of atoms){const v=a[k];if(v<lo)lo=v;if(v>hi)hi=v}return hi-lo},H=counts.H||0;
 return{atomCount:atoms.length,hydrogens:H,heavyAtoms:atoms.length-H,formula,counts,molecularWeight:mw,unknownElements:unknown,radiusOfGyration:rg,rgMassWeighted:complete,extent:{x:ext('x'),y:ext('y'),z:ext('z')},
  hydrogensLikelyMissing:H===0&&(counts.C||0)>0,residues:new Set(atoms.filter(a=>!a.hetero&&a.resName).map(a=>`${a.chain}:${a.resSeq}${a.iCode||''}`)).size,chains:new Set(atoms.filter(a=>a.chain).map(a=>a.chain)).size}}

// Bond-order and hybridization inference from bond lengths + neighbours (heuristic; PDB/mmCIF coordinates do not store bond orders),
// and a consistency check of measured angles against the ideal angle of each atom's hybridization.

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
function assignBondOrders(atoms,bonds){const adj=atoms.map(()=>[]),key=new Map();bonds.forEach(([i,j],k)=>{adj[i].push(j);adj[j].push(i);key.set(i+','+j,k);key.set(j+','+i,k)});
 const ord=bonds.map(([i,j])=>lenOrder(elOf(atoms[i]),elOf(atoms[j]),distance(atoms[i],atoms[j])));
 ringPass(atoms,bonds,ord);
 atoms.forEach((a,c)=>{if(elOf(a)!=='C'||adj[c].length<3)return;const t=adj[c].filter(o=>elOf(atoms[o])==='O'&&adj[o].length===1&&ord[key.get(c+','+o)]===2);if(t.length>=2)t.forEach(o=>{ord[key.get(c+','+o)]=1.5})}); // carboxylate-type resonance
 return ord}
// Hybridization for C, N, O only (others null). N/O attached to an unsaturated carbon are treated as sp2 (conjugation, e.g. amide N, phenol O).
function hybridization(atoms,bonds,orders){const adj=atoms.map(()=>[]);bonds.forEach(([i,j],k)=>{adj[i].push([j,orders[k]]);adj[j].push([i,orders[k]])});
 const unsatC=i=>elOf(atoms[i])==='C'&&adj[i].some(([,o])=>o>1);
 return atoms.map((a,i)=>{const e=elOf(a),nb=adj[i];if(!nb.length)return null;const dbl=nb.filter(([,o])=>o===2).length,trp=nb.some(([,o])=>o===3),aro=nb.some(([,o])=>o===1.5);
  if(e==='C')return trp||dbl>=2?'sp':dbl||aro?'sp2':'sp3';if(e==='N')return trp?'sp':dbl||aro?'sp2':nb.some(([j])=>unsatC(j))?'sp2':'sp3';
  if(e==='O')return trp?'sp':dbl||aro?'sp2':nb.some(([j])=>unsatC(j))?'sp2':'sp3';return null})}
const IDEAL={sp:180,sp2:120,sp3:109.5},TOL={sp:5,sp2:12,sp3:16};
function geometryCheck(atoms,bonds,hyb){const adj=atoms.map(()=>[]);bonds.forEach(([i,j])=>{adj[i].push(j);adj[j].push(i)});const outliers=[];let max=0,n=0;
 atoms.forEach((a,c)=>{const h=hyb[c];if(!h||adj[c].length<2)return;const nb=adj[c];for(let x=0;x<nb.length;x++)for(let y=x+1;y<nb.length;y++){const ang=angle(atoms[nb[x]],a,atoms[nb[y]]),dev=Math.abs(ang-IDEAL[h]);n++;if(dev>max)max=dev;
  if(dev>TOL[h])outliers.push({atom:c,hybridization:h,neighbors:[nb[x],nb[y]],angle:ang,ideal:IDEAL[h],deviation:dev})}});return{checked:n,maxDeviation:max,outliers,tolerances:TOL}}

// Authoritative-by-template bond orders for the 20 standard amino acids (residue dictionary), length heuristic otherwise.
const STD=new Set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(' '));
const AROM={PHE:['CG','CD1','CD2','CE1','CE2','CZ'],TYR:['CG','CD1','CD2','CE1','CE2','CZ'],HIS:['CG','ND1','CD2','CE1','NE2'],TRP:['CG','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2']};
const DELOC={ASP:['CG','OD1','OD2'],GLU:['CD','OE1','OE2'],ARG:['CZ','NH1','NH2','NE']},DBL={ASN:[['CG','OD1']],GLN:[['CD','OE1']]};
// returns {orders, source}: source[k] is 'template' (standard residue dictionary; delocalised groups = 1.5) or 'length' (inferred from bond length only).
function assignBondOrdersDetailed(atoms,bonds){const lenOrd=assignBondOrders(atoms,bonds),key=a=>a.chain+'|'+a.resSeq+'|'+(a.iCode||''),res=new Map(),orders=[],source=[];
 const std=a=>!a.hetero&&STD.has(a.resName);atoms.forEach(a=>{if(std(a)){const k=key(a);if(!res.has(k))res.set(k,{n:a.resName,names:new Set()});res.get(k).names.add(a.name)}});
 bonds.forEach(([i,j],k)=>{const a=atoms[i],b=atoms[j];if(std(a)&&std(b)){let o=1;if(key(a)===key(b)){const r=res.get(key(a)),n=r.n,A=a.name,B=b.name,has=(x,y)=>(A===x&&B===y)||(A===y&&B===x);
   if(has('C','O'))o=r.names.has('OXT')?1.5:2;else if(has('C','OXT'))o=1.5;else if(AROM[n]&&AROM[n].includes(A)&&AROM[n].includes(B))o=1.5;else if(DELOC[n]&&DELOC[n].includes(A)&&DELOC[n].includes(B))o=1.5;else if(DBL[n]&&DBL[n].some(([x,y])=>has(x,y)))o=2}
  orders.push(o);source.push('template')}else{orders.push(lenOrd[k]);source.push('length')}});return{orders,source}}

// Session-file validation. Treats input as untrusted data: never evaluated, size-checked before parsing, strictly typed.
const SESSION_MAX_CHARS=6_000_000; // 5 MB structure limit + JSON overhead
const STYLES=['bs','sf','st','wf','bb','tb'],fin=v=>typeof v==='number'&&Number.isFinite(v),vec3=a=>Array.isArray(a)&&a.length===3&&a.every(fin);
const bad=m=>({ok:false,error:'Not a valid LAURA session: '+m});
function validateSession(text,{maxChars=SESSION_MAX_CHARS}={}){if(typeof text!=='string')return bad('the file could not be read as text.');
 if(text.length>maxChars)return bad(`the file is ${(text.length/1e6).toFixed(1)} MB; the limit is ${(maxChars/1e6).toFixed(1)} MB.`);
 let S;try{S=JSON.parse(text)}catch(e){return bad('it is not valid JSON (the file may be truncated or is not a session file).')}
 if(!S||typeof S!=='object'||Array.isArray(S))return bad('the top level must be an object.');
 if(typeof S.laura!=='string'||!/^\d+\.\d+$/.test(S.laura))return bad('missing or malformed "laura" version field; this does not look like a LAURA session.');
 const[major,minor]=S.laura.split('.').map(Number);if(major!==0||minor<3||minor>4)return bad(`it was written by LAURA session version ${S.laura}, which this version cannot read (supported: 0.3–0.4).`);
 const st=S.structure;if(!st||typeof st!=='object')return bad('the "structure" section is missing.');
 if(typeof st.id!=='string'||typeof st.name!=='string')return bad('structure.id and structure.name must be text.');
 if(st.builtin!==true&&typeof st.raw!=='string')return bad(st.rawOmitted?'the structure file was too large to package; re-import the original file.':'the structure file is not packaged and it is not a built-in molecule.');
 const set={style:'bs',hydrogens:false,solvent:true,labels:false,elements:true,lonepairs:false,...(S.settings||{})};if(!STYLES.includes(set.style))return bad(`unknown display style "${String(set.style).slice(0,20)}".`);
 for(const k of['hydrogens','solvent','labels','elements','lonepairs'])if(typeof set[k]!=='boolean')return bad(`settings.${k} must be true or false.`);
 const sel=S.selection===undefined?[]:S.selection;if(!Array.isArray(sel)||!sel.every(x=>Array.isArray(x)&&x.length===3&&typeof x[0]==='string'&&fin(x[1])&&typeof x[2]==='string'))return bad('selection entries must be [chain, serial, atom name].');
 const con={ligand:'',cutoff:'4.0',annotated:'',source:'',result:null,...(S.contacts||{})};for(const k of['ligand','cutoff','annotated','source'])if(typeof con[k]!=='string')return bad(`contacts.${k} must be text.`);
 if(S.camera!==undefined&&S.camera!==null&&!(S.camera&&vec3(S.camera.p)&&vec3(S.camera.t)))return bad('camera.p and camera.t must be three finite numbers.');
 return{ok:true,session:{...S,structure:st,settings:set,selection:sel,contacts:con,camera:S.camera||null}}}

// Pure UI decision helpers (testable without a DOM).
function shortcutAllowed(e,{modalOpen=false}={}){if(e.ctrlKey||e.metaKey||e.altKey)return false;const k=String(e.key||'').toLowerCase();if(k==='escape')return true;if(modalOpen)return false;
 const t=e.target||{},tag=String(t.tagName||'').toUpperCase();if(t.isContentEditable||tag==='TEXTAREA'||tag==='SELECT')return false;if(tag==='INPUT'&&t.type!=='range')return false;
 if(k===' '&&(tag==='BUTTON'||tag==='SUMMARY'||tag==='A'||(t.getAttribute&&t.getAttribute('role')==='separator')))return false;return true}
// Display-name sanitising for renamed imports: trims, collapses whitespace, removes control characters, caps length. Returns null if nothing is left.
function cleanName(s,max=80){const t=String(s??'').replace(/[\u0000-\u001f\u007f]/g,' ').replace(/\s+/g,' ').trim().slice(0,max).trim();return t.length?t:null}

// Camera framing: centre and distance that fit a set of points inside the narrower of the vertical/horizontal field of view.
function framing(points,{fov=40,aspect=1,minRadius=1.5,margin=1.15}={}){if(!points||!points.length)return null;const lo=[Infinity,Infinity,Infinity],hi=[-Infinity,-Infinity,-Infinity],K=['x','y','z'];
 for(const p of points)K.forEach((k,i)=>{if(!Number.isFinite(p[k]))throw new RangeError('Non-finite coordinate in framing input.');lo[i]=Math.min(lo[i],p[k]);hi[i]=Math.max(hi[i],p[k])});
 const c={x:(lo[0]+hi[0])/2,y:(lo[1]+hi[1])/2,z:(lo[2]+hi[2])/2};let r=minRadius;for(const p of points)r=Math.max(r,Math.hypot(p.x-c.x,p.y-c.y,p.z-c.z));
 const v=fov*Math.PI/360,h=Math.atan(Math.tan(v)*aspect);return{center:c,radius:r,distance:margin*r/Math.sin(Math.min(v,h))}}

// Flat "2D view" of 3D coordinates: orthographic projection onto the best-fit plane (principal-component analysis).
// This is NOT a chemical structure diagram: it preserves distances only for planar molecules; non-planar atoms overlap.
function eigSymN(N){const n=N.length,A=N.map(r=>[...r]),V=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>i===j?1:0));
 for(let sw=0;sw<80;sw++){let off=0;for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)off+=A[i][j]**2;if(off<1e-28)break;
  for(let p=0;p<n-1;p++)for(let q=p+1;q<n;q++){if(Math.abs(A[p][q])<1e-300)continue;const th=(A[q][q]-A[p][p])/(2*A[p][q]),t=(th>=0?1:-1)/(Math.abs(th)+Math.sqrt(th*th+1)),c=1/Math.sqrt(t*t+1),s=t*c;
   for(let k=0;k<n;k++){const x=A[k][p],y=A[k][q];A[k][p]=c*x-s*y;A[k][q]=s*x+c*y}for(let k=0;k<n;k++){const x=A[p][k],y=A[q][k];A[p][k]=c*x-s*y;A[q][k]=s*x+c*y}for(let k=0;k<n;k++){const x=V[k][p],y=V[k][q];V[k][p]=c*x-s*y;V[k][q]=s*x+c*y}}}
 return{values:A.map((r,i)=>r[i]),vectors:V}}
function projectPlane(points){const n=points.length;if(!n)return{xy:[],depth:[],variance:[0,0,0]};const c=['x','y','z'].map(k=>points.reduce((s,p)=>s+p[k],0)/n),C=[[0,0,0],[0,0,0],[0,0,0]];
 for(const p of points){const d=[p.x-c[0],p.y-c[1],p.z-c[2]];for(let i=0;i<3;i++)for(let j=0;j<3;j++)C[i][j]+=d[i]*d[j]/n}
 const{values,vectors}=eigSymN(C),ord=[0,1,2].sort((a,b)=>values[b]-values[a]),ax=ord.map(k=>[vectors[0][k],vectors[1][k],vectors[2][k]]);
 const dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2],xy=[],depth=[];
 for(const p of points){const d=[p.x-c[0],p.y-c[1],p.z-c[2]];xy.push([dot(d,ax[0]),dot(d,ax[1])]);depth.push(dot(d,ax[2]))}
 return{xy,depth,variance:ord.map(k=>Math.max(0,values[k]))}}

// Screen-space label placement. Greedy by priority; a label is dropped if it overlaps one already placed or lies off-screen.
function declutter(items,{width=1000,height=700,pad=3,cap=80}={}){const order=items.map((it,i)=>i).sort((a,b)=>((items[a].prio??0)-(items[b].prio??0))||a-b),placed=[],grid=new Map(),CELL=80;
 const cells=r=>{const out=[];for(let gx=Math.floor(r[0]/CELL);gx<=Math.floor(r[2]/CELL);gx++)for(let gy=Math.floor(r[1]/CELL);gy<=Math.floor(r[3]/CELL);gy++)out.push(gx+','+gy);return out};
 for(const i of order){if(placed.length>=cap)break;const it=items[i];if(!(it.x>=0&&it.x<=width&&it.y>=0&&it.y<=height))continue;
  const r=[it.x-it.w/2-pad,it.y-it.h/2-pad,it.x+it.w/2+pad,it.y+it.h/2+pad],cs=cells(r);let hit=false;
  for(const k of cs){const L=grid.get(k);if(L&&L.some(q=>!(r[2]<=q[0]||r[0]>=q[2]||r[3]<=q[1]||r[1]>=q[3]))){hit=true;break}}
  if(hit)continue;placed.push(i);for(const k of cs){if(!grid.has(k))grid.set(k,[]);grid.get(k).push(r)}}
 return placed}
// Label text: real atom name (+ residue) for imported/protein atoms; element symbol for generated molecules whose "names" are just running numbers (C1, O2...).
function labelText(a){const el=a.element||a.el||'';if(a.autoName)return el;const res=a.resName&&a.resSeq!==undefined&&a.resSeq!==0?` ${a.resName}${a.resSeq}${a.iCode||''}${a.chain?':'+a.chain:''}`:'';return (a.name||el)+res}

return{parsePDB,parseMmCIF,parseStructure,ParseError,residueKey,distance,angle,dihedral,ligandContacts,alignGlobal,chainSequence,phiPsi,caContactMap,superpose,compareStructures,zmatrix,peptide,transformAtoms,RAMA_REGIONS,ramaRegion,ramaSummary,hillFormula,molProps,assignBondOrders,assignBondOrdersDetailed,hybridization,geometryCheck,validateSession,SESSION_MAX_CHARS,shortcutAllowed,framing,declutter,labelText,projectPlane,cleanName,models:{seq3,MODEL_SEQ,alphaHelix,helixLoopHelix,betaHairpin,ramaTestPeptide,antiparallelHelixPair}}})();
function parseXYZ(t){const L=t.split(/\r?\n/),n=parseInt(L[0]);if(!(n>0))throw Error('First line of an XYZ file must be the atom count.');const atoms=[];
for(let i=0;i<n;i++){const p=(L[i+2]||'').trim().split(/\s+/),x=+p[1],y=+p[2],z=+p[3];if(p.length<4||[x,y,z].some(isNaN))throw Error('Bad XYZ line '+(i+3));const el=p[0][0].toUpperCase()+p[0].slice(1).toLowerCase();atoms.push({el,x,y,z,serial:i+1,name:el+(i+1),res:'',resi:0,chain:''})}return{atoms}}
function parseFile(name,text){const ext=name.split('.').pop().toLowerCase();let r;
if(ext==='xyz'){if(text.length>60e6)throw Error('File exceeds the 60 MB limit.');r=parseXYZ(text);r.atoms=r.atoms.map(a=>({...a,element:a.el,resName:a.res,resSeq:0,iCode:'',hetero:false}))}
else{r=Core.parseStructure(name,text);r.atoms=r.atoms.map(a=>({...a,el:a.element,res:a.resName,resi:a.resSeq}))}
r.raw=text;r.ext=ext;r.warnings=r.warnings||[];return r}
function computeBonds(A){const h=new Map(),B=[],c=2.6,k=(x,y,z)=>x+','+y+','+z;A.forEach((a,i)=>{const key=k(Math.floor(a.x/c),Math.floor(a.y/c),Math.floor(a.z/c));(h.get(key)||h.set(key,[]).get(key)).push(i)});
A.forEach((a,i)=>{const cx=Math.floor(a.x/c),cy=Math.floor(a.y/c),cz=Math.floor(a.z/c);for(let dx=-1;dx<=1;dx++)for(let dy=-1;dy<=1;dy++)for(let dz=-1;dz<=1;dz++)for(const j of h.get(k(cx+dx,cy+dy,cz+dz))||[]){if(j<=i)continue;const b=A[j];if(a.el==='H'&&b.el==='H')continue;const d=Math.hypot(a.x-b.x,a.y-b.y,a.z-b.z);if(d>.4&&d<cov(a.el)+cov(b.el)+.4)B.push([i,j])}});return B}
