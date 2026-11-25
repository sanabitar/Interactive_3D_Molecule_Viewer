// js/app/ligand-contacts.js: Ligand contacts: list ligands, find the residues in contact, show, copy and export them.
// ---------- ligand contacts (Core.ligandContacts) ----------
function popLig(){const L=new Map();mol.atoms.forEach(a=>{if(a.hetero&&a.resName!=='HOH')L.set(Core.residueKey(a),a.resName)});const s=$('#lig');
s.innerHTML=[...L].map(([k,n])=>`<option value="${esc(k)}">${esc(k)} (${esc(n)})</option>`).join('');$('#rc').disabled=!L.size;if(!L.size)s.innerHTML='<option>No ligand (HETATM) records</option>'}
function runContacts(){if(!mol)return;const ann=new Map(),src=$('#asrc').value.trim()||'user-supplied, source not stated';$('#ann').value.split(/[,\s]+/).filter(Boolean).forEach(x=>ann.set(x,src));
try{lastC=Core.ligandContacts(mol.atoms,{ligandKey:$('#lig').value,cutoff:parseFloat($('#cut').value),annotated:ann})}catch(e){return toast(e.message)}
const cs=new Set(lastC.contacts.map(c=>c.residue)),as=new Set(lastC.contacts.filter(c=>c.evidence.length>1).map(c=>c.residue)),k=$('#lig').value;hl=[];
mol.atoms.forEach((a,i)=>{const r=Core.residueKey(a);if(a.hetero&&r===k)hl.push([i,0xff7a00]);else if(!a.hetero&&cs.has(r))hl.push([i,as.has(r)?0xff00ff:0xff7a00])});
$('#cres').innerHTML=`<b>${lastC.contacts.length}</b> residues within ${lastC.params.cutoffAngstrom} Å of ${esc(k)}<table style="width:100%;font-size:12px">`+lastC.contacts.map(c=>`<tr><td>${esc(c.residue)} ${esc(c.resName)}</td><td>${c.minDistance.toFixed(2)} Å</td><td>${c.evidence.length>1?'contact + annotated*':'contact'}</td></tr>`).join('')+'</table>'+
(lastC.annotatedNotInContact.length?`<div>Annotated, not in contact: ${lastC.annotatedNotInContact.map(x=>esc(x.residue)+(x.presentInStructure?'':' (not in structure)')).join(', ')}</div>`:'')+(ann.size?'<div>*Annotations are user-supplied; LAURA does not verify them.</div>':'')+`<div>${esc(lastC.caveat)}</div>`;makeContactLabels(k);drawSel()}
$('#ccl').onclick=()=>{hl=[];cLab=[];lastC=null;$('#cres').innerHTML='';if(mol)drawSel()};
$('#rc').onclick=runContacts;
$('#cj').onclick=()=>lastC?save('ligand_contacts_'+slug()+'.json',JSON.stringify({tool:'LAURA core contacts (browser)',structure:mol.entry.name,source:mol.entry.src,exportedAt:new Date().toISOString(),result:lastC},null,1)):toast('Run "Find contacts" first.');
$('#cc').onclick=()=>{if(!lastC)return toast('Run "Find contacts" first.');const q=v=>'"'+String(v??'').replace(/"/g,'""')+'"';save('ligand_contacts_'+slug()+'.csv',['residue,resName,minDistance_A,proteinAtom,ligandAtom,evidence,annotationSource,cutoff_A,ligand'].concat(lastC.contacts.map(c=>[c.residue,c.resName,c.minDistance.toFixed(3),c.proteinAtom,c.ligandAtom,c.evidence.join('+'),c.annotationSource,lastC.params.cutoffAngstrom,lastC.params.ligand].map(q).join(','))).join('\n'),'text/csv')};
