// Method B (ligand-associated site): purely geometric. Contact != catalytic. See caveat in output.
import {distance} from './geometry.mjs';import {residueKey} from './structure.mjs';
export function ligandContacts(atoms,{ligandKey,cutoff=4.0,annotated=new Map()}={}){
 if(!(cutoff>0))throw new RangeError('cutoff must be > 0 Angstrom');
 const lig=atoms.filter(a=>a.hetero&&residueKey(a)===ligandKey);if(!lig.length)throw new Error(`Ligand ${ligandKey} not found among HETATM records.`);
 const prot=atoms.filter(a=>!a.hetero),best=new Map();
 for(const p of prot)for(const q of lig){const d=distance(p,q);if(d<=cutoff){const k=residueKey(p),b=best.get(k);if(!b||d<b.minDistance)best.set(k,{residue:k,resName:p.resName,minDistance:d,proteinAtom:p.name,ligandAtom:q.name})}}
 const contacts=[...best.values()].sort((a,b)=>a.minDistance-b.minDistance).map(c=>({...c,evidence:annotated.has(c.residue)?['geometric_contact','annotated_catalytic']:['geometric_contact'],annotationSource:annotated.get(c.residue)||null}));
 const present=new Set(prot.map(residueKey));
 return{method:'Heavy/any-atom minimum distance, ligand vs polymer ATOM records',params:{cutoffAngstrom:cutoff,ligand:ligandKey,ligandResName:lig[0].resName,altConfsHandled:'as parsed'},contacts,
  annotatedNotInContact:[...annotated.keys()].filter(k=>!best.has(k)).map(k=>({residue:k,presentInStructure:present.has(k),annotationSource:annotated.get(k)})),
  caveat:'Contact residues are geometric neighbours of a bound ligand. They are not established catalytic residues, and the ligand may not be a substrate or biologically relevant.'}}
