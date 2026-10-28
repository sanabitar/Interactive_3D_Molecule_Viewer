// Structure parsing + validation. Pure functions, no DOM. Units: Angstrom.
export class ParseError extends Error{constructor(m,line){super(line?`Line ${line}: ${m}`:m);this.name='ParseError';this.line=line}}
const cap=e=>e?e[0].toUpperCase()+e.slice(1).toLowerCase():'';
export const residueKey=a=>`${a.chain}:${a.resSeq}${a.iCode||''}`; // chain + number + insertion code
const num=(s,line,what)=>{if(s===undefined||String(s).trim()===''||!Number.isFinite(Number(s)))throw new ParseError(`invalid ${what} "${String(s).trim()}"`,line);return Number(s)};
function finish(atoms,warnings,format){if(!atoms.length)throw new ParseError('No ATOM/HETATM records found.');
 return{format,atoms,warnings,chains:[...new Set(atoms.map(a=>a.chain))],residues:new Set(atoms.map(a=>residueKey(a)+a.resName)).size}}
const TWO_LETTER=new Set(['ZN','FE','MG','CA','NA','CL','MN','CU','CO','NI','CD','HG','BR','SE','LI','AL','SI','BE','BA','SR','CS','RB','AU','AG','PT','PB','SN','MO','CR','TI']);
export function parsePDB(text){const atoms=[],warnings=[],seen=new Set(),L=text.split(/\r?\n/);let alt=null,skippedAlt=0,inferred=0,models=0,done=false,dup=0;
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
export function parseMmCIF(text){const L=text.split(/\r?\n/),atoms=[],warnings=[];let i=0;
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
export function parseStructure(name,text,{maxBytes=60e6}={}){if(text.length>maxBytes)throw new ParseError('File exceeds the size limit.');
 if(/[\x00-\x08]/.test(text.slice(0,2000)))throw new ParseError('Binary data detected; expected a text structure file.');
 const ext=name.split('.').pop().toLowerCase();if(ext==='pdb'||ext==='ent')return parsePDB(text);if(ext==='cif'||ext==='mmcif')return parseMmCIF(text);
 throw new ParseError(`Unsupported format ".${ext}" (supported here: pdb, ent, cif, mmcif).`)}
