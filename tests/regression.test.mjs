import test from 'node:test';import assert from 'node:assert/strict';import fs from 'fs';
import {parsePDB,parseMmCIF,ParseError} from '../src/core/structure.mjs';import {assignBondOrders,assignBondOrdersDetailed} from '../src/core/bonds.mjs';import {peptide} from '../src/core/build.mjs';
import {validateSession,SESSION_MAX_CHARS} from '../src/core/session.mjs';import {shortcutAllowed} from '../src/core/ui.mjs';import {framing} from '../src/core/camera.mjs';
const FX=f=>fs.readFileSync(new URL('./fixtures/'+f,import.meta.url),'utf8');
test('D7 regression: extra PDB models are reported, not silently dropped',()=>{const r=parsePDB(FX('two_models.pdb'));assert.equal(r.atoms.length,10);assert.ok(r.warnings.some(w=>/1 additional model/.test(w)),JSON.stringify(r.warnings))});
test('D8 regression: a record truncated inside the coordinates is rejected with its line number',()=>{assert.throws(()=>parsePDB(FX('bad_truncated.pdb')),e=>e instanceof ParseError&&/truncated/.test(e.message)&&e.line===60)});
test('duplicate atom records are kept but reported',()=>{const t=FX('alpha12.pdb'),l=t.split('\n')[1],r=parsePDB(t+l+'\n');assert.equal(r.atoms.length,61);assert.ok(r.warnings.some(w=>/duplicate/.test(w)))});
test('D13 regression: HETATM ion without an element column keeps its two-letter element; Calpha does not become calcium',()=>{
 const zn='HETATM  210 ZN    ZN A 102     -10.000 -10.000 -10.000  1.00  0.00',ca='ATOM      2  CA  ALA A   1       0.000   0.000   1.458  1.00  0.00',caIon='HETATM  211 CA    CA A 103       1.000   1.000   1.000  1.00  0.00';
 const r=parsePDB([zn,ca,caIon].join('\n'));assert.deepEqual(r.atoms.map(a=>a.element),['Zn','C','Ca']);assert.ok(r.warnings.some(w=>/inferred/.test(w)))});
const CIF=h=>`data_x\nloop_\n${h.map(x=>'_atom_site.'+x).join('\n')}\n`;const H=['group_PDB','id','type_symbol','label_atom_id','label_alt_id','label_comp_id','label_asym_id','label_seq_id','pdbx_PDB_ins_code','Cartn_x','Cartn_y','Cartn_z','pdbx_PDB_model_num'];
test('mmCIF: wrapped rows are joined; truncated final row and extra models are reported',()=>{
 const r=parseMmCIF(CIF(H)+'ATOM 1 C CA . ALA A 1 ?\n 1.0 2.0 3.0 1\nATOM 2 C CA . ALA A 1 ? 4.0 5.0 6.0 2\n#\n');assert.equal(r.atoms.length,1);assert.equal(r.atoms[0].z,3);assert.ok(r.warnings.some(w=>/additional model/.test(w)));
 assert.throws(()=>parseMmCIF(CIF(H)+'ATOM 1 C CA . ALA A 1 ?\n'),/Truncated mmCIF/);assert.throws(()=>parseMmCIF(CIF(H)+'ATOM 1 C CA . ALA A 1 ? 1 2 3 1 extra\n'),/more fields/)});
function rng(seed){return()=>{seed|=0;seed=seed+0x6D2B79F5|0;let t=Math.imul(seed^seed>>>15,1|seed);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296}}
const gauss=r=>Math.sqrt(-2*Math.log(r()+1e-12))*Math.cos(2*Math.PI*r());
test('D11 quantified: length-only bond orders misclassify noisy protein bonds; residue templates do not',(t)=>{
 const P=peptide(8,{phi:-57,psi:-47}).map(a=>({...a,element:a.el,resName:a.res,resSeq:a.resi,iCode:'',hetero:false})),bonds=[];for(let i=0;i<P.length;i++)for(let j=i+1;j<P.length;j++)if(Math.hypot(P[i].x-P[j].x,P[i].y-P[j].y,P[i].z-P[j].z)<1.8)bonds.push([i,j]);
 const truth=assignBondOrdersDetailed(P,bonds).orders;const r=rng(42),noisy=P.map(a=>({...a,x:a.x+.06*gauss(r),y:a.y+.06*gauss(r),z:a.z+.06*gauss(r)}));
 const lenOnly=assignBondOrders(noisy,bonds),tmpl=assignBondOrdersDetailed(noisy,bonds),bad=(o)=>o.filter((v,k)=>v!==truth[k]).length;
 t.diagnostic(`${bonds.length} bonds, coordinate noise sigma 0.06 A: length-only misclassified ${bad(lenOnly)}, template misclassified ${bad(tmpl.orders)}`);assert.equal(bad(tmpl.orders),0);assert.ok(bad(lenOnly)>0);assert.ok(tmpl.source.every(s=>s==='template'))});
const A=(chain,resSeq,resName,name,x,y,z,hetero=false)=>({chain,resSeq,iCode:'',resName,name,element:name[0]==='O'?'O':name[0]==='N'?'N':'C',x,y,z,hetero});
test('residue-template bond orders: aromatic ring, carboxylate, amide, guanidinium, C-terminus; ligands fall back to length',()=>{
 const ring=['CG','CD1','CE1','CZ','CE2','CD2'].map((n,i)=>A('A',1,'PHE',n,1.4*Math.cos(i*Math.PI/3),1.4*Math.sin(i*Math.PI/3),0)),CB=A('A',1,'PHE','CB',3,0,0),at=[...ring,CB];
 let b=[[0,1],[1,2],[2,3],[3,4],[4,5],[5,0],[0,6]];let d=assignBondOrdersDetailed(at,b);assert.deepEqual(d.orders,[1.5,1.5,1.5,1.5,1.5,1.5,1],'PHE ring 1.5 each, CB-CG single');
 d=assignBondOrdersDetailed([A('A',2,'ASP','CB',0,0,0),A('A',2,'ASP','CG',1.5,0,0),A('A',2,'ASP','OD1',2.2,1,0),A('A',2,'ASP','OD2',2.2,-1,0)],[[0,1],[1,2],[1,3]]);assert.deepEqual(d.orders,[1,1.5,1.5],'ASP carboxylate');
 d=assignBondOrdersDetailed([A('A',3,'ASN','CG',0,0,0),A('A',3,'ASN','OD1',1.2,0,0),A('A',3,'ASN','ND2',-1,1,0)],[[0,1],[0,2]]);assert.deepEqual(d.orders,[2,1],'ASN amide C=O double');
 d=assignBondOrdersDetailed([A('A',4,'ARG','NE',0,0,0),A('A',4,'ARG','CZ',1.3,0,0),A('A',4,'ARG','NH1',2,1,0),A('A',4,'ARG','NH2',2,-1,0),A('A',4,'ARG','CD',-1,1,0)],[[0,1],[1,2],[1,3],[0,4]]);assert.deepEqual(d.orders,[1.5,1.5,1.5,1],'ARG guanidinium delocalised, NE-CD single');
 d=assignBondOrdersDetailed([A('A',5,'GLY','C',0,0,0),A('A',5,'GLY','O',1.23,0,0),A('A',5,'GLY','OXT',-1.2,0.5,0)],[[0,1],[0,2]]);assert.deepEqual(d.orders,[1.5,1.5],'C-terminal carboxylate');
 d=assignBondOrdersDetailed([A('A',6,'ALA','C',0,0,0),A('A',6,'ALA','O',1.23,0,0),A('A',7,'ALA','N',-1.3,0,0)],[[0,1],[0,2]]);assert.deepEqual(d.orders,[2,1],'backbone C=O double, peptide bond single');
 const lig=[A('A',9,'LIG','C1',0,0,0,true),A('A',9,'LIG','C2',1.33,0,0,true)];d=assignBondOrdersDetailed(lig,[[0,1]]);assert.deepEqual([d.orders[0],d.source[0]],[2,'length'],'ligand uses (disclosed) length inference')});
const S0=()=>({laura:'0.3',savedAt:'x',structure:{id:'water',name:'Water',builtin:true,raw:null},settings:{style:'bs',hydrogens:true,solvent:true,labels:false},selection:[['A',1,'CA']],contacts:{ligand:'',cutoff:'4.0',annotated:'',source:'',result:null},camera:{p:[0,0,5],t:[0,0,0]}});
test('D6 regression: session validation gives clear errors and keeps old sessions readable',()=>{
 const E=o=>validateSession(typeof o==='string'?o:JSON.stringify(o));assert.equal(E(S0()).ok,true);
 assert.match(E('{"laura":"0.3","structure":{"id"').error,/not valid JSON/);assert.match(E('hello').error,/not valid JSON/);assert.match(E({laura:'9.9',structure:{id:'a',name:'a',builtin:true}}).error,/cannot read|version 9\.9/);
 assert.match(E({laura:'0.3'}).error,/structure/);assert.match(E({...S0(),settings:{style:'zz'}}).error,/unknown display style/);assert.match(E({...S0(),selection:[['A','x','CA']]}).error,/selection entries/);
 assert.match(E({...S0(),camera:{p:[0,0],t:[0,0,0]}}).error,/camera/);assert.match(E({...S0(),structure:{id:'a',name:'a'}}).error,/not packaged/);assert.match(E({...S0(),structure:{id:'a',name:'a',rawOmitted:true}}).error,/too large/);
 const big='{"laura":"0.3","pad":"'+'x'.repeat(SESSION_MAX_CHARS)+'"}';assert.match(validateSession(big).error,/limit is 6\.0 MB/,'size checked before JSON.parse');
 const old={laura:'0.3',structure:{id:'water',name:'Water',builtin:true}};const r=E(old);assert.equal(r.ok,true);assert.equal(r.session.settings.style,'bs');assert.deepEqual(r.session.selection,[]);assert.equal(r.session.camera,null);
 const inj=E({...S0(),structure:{id:'x',name:'<img onerror=alert(1)>',builtin:true}});assert.equal(inj.ok,true,'text is data; rendering must escape it (checked in browser tests)')});
test('D3/D4 regression: shortcut gating',()=>{const e=(o)=>({key:'t',ctrlKey:false,metaKey:false,altKey:false,target:{tagName:'DIV'},...o});
 assert.equal(shortcutAllowed(e({})),true);assert.equal(shortcutAllowed(e({target:{tagName:'TEXTAREA'}})),false);assert.equal(shortcutAllowed(e({target:{tagName:'INPUT',type:'text'}})),false);assert.equal(shortcutAllowed(e({target:{tagName:'INPUT',type:'range'}})),true);
 assert.equal(shortcutAllowed(e({target:{tagName:'SELECT'}})),false);assert.equal(shortcutAllowed(e({target:{tagName:'DIV',isContentEditable:true}})),false);for(const m of['ctrlKey','metaKey','altKey'])assert.equal(shortcutAllowed(e({[m]:true})),false,m);
 assert.equal(shortcutAllowed(e({key:' ',target:{tagName:'BUTTON'}})),false,'Space on a focused button activates it');assert.equal(shortcutAllowed(e({key:' '})),true);assert.equal(shortcutAllowed(e({key:'p'}),{modalOpen:true}),false);
 assert.equal(shortcutAllowed(e({key:'Escape',target:{tagName:'TEXTAREA'}})),true,'Esc always allowed');assert.equal(shortcutAllowed(e({key:'Escape',ctrlKey:true})),false)});
test('camera framing: analytic distance, aspect ratio, degenerate cases',()=>{const P=[{x:-10,y:0,z:0},{x:10,y:0,z:0}],near=(a,b,m)=>assert.ok(Math.abs(a-b)<1e-9,m+': '+a+' vs '+b);
 let f=framing(P,{fov:40,aspect:1,margin:1});near(f.distance,10/Math.sin(20*Math.PI/180),'square view');near(f.radius,10,'radius');near(f.center.x,0,'centre');
 f=framing(P,{fov:40,aspect:.5,margin:1});near(f.distance,10/Math.sin(Math.atan(Math.tan(20*Math.PI/180)*.5)),'narrow view needs more distance');assert.ok(f.distance>framing(P,{fov:40,aspect:1,margin:1}).distance);
 f=framing([{x:3,y:4,z:5}],{minRadius:2.5,margin:1});near(f.radius,2.5,'single point uses minRadius');assert.deepEqual(f.center,{x:3,y:4,z:5});assert.equal(framing([]),null);assert.throws(()=>framing([{x:NaN,y:0,z:0}]),RangeError)});
import {molProps} from '../src/core/properties.mjs';import {compareStructures} from '../src/core/analysis.mjs';
test('D16 regression: structure statistics and RMSD do not use argument spreading (no stack overflow at 300k atoms)',()=>{
 const N=300000,big=Array.from({length:N},(_,i)=>({element:i%3?'C':'N',x:i*0.01,y:(i%7),z:(i%11),hetero:false,chain:'A',resSeq:i,iCode:'',resName:'ALA'}));
 const p=molProps(big);assert.equal(p.atomCount,N);assert.ok(Math.abs(p.extent.x-(N-1)*0.01)<1e-6);
 const mv=big.map(a=>({...a,x:a.x+1}));const r=compareStructures(big,mv,{mode:'ord'});assert.equal(r.nPairs,N);assert.ok(r.rmsd<1e-6,'pure translation -> RMSD 0');assert.ok(r.maxDeviation<1e-6)});
