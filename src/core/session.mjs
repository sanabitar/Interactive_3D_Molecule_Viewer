// Session-file validation. Treats input as untrusted data: never evaluated, size-checked before parsing, strictly typed.
export const SESSION_MAX_CHARS=6_000_000; // 5 MB structure limit + JSON overhead
const STYLES=['bs','sf','st','wf','bb','tb'],fin=v=>typeof v==='number'&&Number.isFinite(v),vec3=a=>Array.isArray(a)&&a.length===3&&a.every(fin);
const bad=m=>({ok:false,error:'Not a valid LAURA session: '+m});
export function validateSession(text,{maxChars=SESSION_MAX_CHARS}={}){if(typeof text!=='string')return bad('the file could not be read as text.');
 if(text.length>maxChars)return bad(`the file is ${(text.length/1e6).toFixed(1)} MB; the limit is ${(maxChars/1e6).toFixed(1)} MB.`);
 let S;try{S=JSON.parse(text)}catch(e){return bad('it is not valid JSON (the file may be truncated or is not a session file).')}
 if(!S||typeof S!=='object'||Array.isArray(S))return bad('the top level must be an object.');
 if(typeof S.laura!=='string'||!/^\d+\.\d+$/.test(S.laura))return bad('missing or malformed "laura" version field; this does not look like a LAURA session.');
 const[major,minor]=S.laura.split('.').map(Number);if(major!==0||minor<3||minor>4)return bad(`it was written by LAURA session version ${S.laura}, which this version cannot read (supported: 0.3–0.4).`);
 const st=S.structure;if(!st||typeof st!=='object')return bad('the "structure" section is missing.');
 if(typeof st.id!=='string'||typeof st.name!=='string')return bad('structure.id and structure.name must be text.');
 if(st.builtin!==true&&typeof st.raw!=='string')return bad(st.rawOmitted?'the structure file was too large to package; re-import the original file.':'the structure file is not packaged and it is not a built-in molecule.');
 const set={style:'bs',hydrogens:true,solvent:true,labels:false,...(S.settings||{})};if(!STYLES.includes(set.style))return bad(`unknown display style "${String(set.style).slice(0,20)}".`);
 for(const k of['hydrogens','solvent','labels'])if(typeof set[k]!=='boolean')return bad(`settings.${k} must be true or false.`);
 const sel=S.selection===undefined?[]:S.selection;if(!Array.isArray(sel)||!sel.every(x=>Array.isArray(x)&&x.length===3&&typeof x[0]==='string'&&fin(x[1])&&typeof x[2]==='string'))return bad('selection entries must be [chain, serial, atom name].');
 const con={ligand:'',cutoff:'4.0',annotated:'',source:'',result:null,...(S.contacts||{})};for(const k of['ligand','cutoff','annotated','source'])if(typeof con[k]!=='string')return bad(`contacts.${k} must be text.`);
 if(S.camera!==undefined&&S.camera!==null&&!(S.camera&&vec3(S.camera.p)&&vec3(S.camera.t)))return bad('camera.p and camera.t must be three finite numbers.');
 return{ok:true,session:{...S,structure:st,settings:set,selection:sel,contacts:con,camera:S.camera||null}}}
