import test from 'node:test';import assert from 'node:assert/strict';
import {parsePDB,parseMmCIF,parseStructure,ParseError,residueKey} from '../src/core/structure.mjs';
import {distance,angle,dihedral} from '../src/core/geometry.mjs';import {ligandContacts} from '../src/core/contacts.mjs';
// SYNTHETIC fixtures (hand-built, not experimental data); exact PDB columns via helper.
const P=(rec,ser,name,alt,res,ch,seq,ic,x,y,z,el)=>rec.padEnd(6)+String(ser).padStart(5)+' '+name.padEnd(4)+alt+res.padStart(3).padEnd(4)+ch+String(seq).padStart(4)+ic+'   '+x.toFixed(3).padStart(8)+y.toFixed(3).padStart(8)+z.toFixed(3).padStart(8)+'  1.00  0.00          '+el.padStart(2);
const PDB=[P('ATOM',1,' N  ',' ','SER','A',52,' ',0,0,0,'N'),P('ATOM',2,' CA ','A','SER','A',52,' ',3,4,0,'C'),P('ATOM',3,' CA ','B','SER','A',52,' ',9,9,9,'C'),
 P('ATOM',4,' CA ',' ','GLY','A',52,'A',20,0,0,'C'),P('HETATM',5,' C1 ',' ','LIG','A',301,' ',0,3.5,0,'C'),P('HETATM',6,' O  ',' ','HOH','A',401,' ',50,50,50,'O'),'ENDMDL',P('ATOM',7,' CA ',' ','ALA','A',99,' ',1,1,1,'C')].join('\n');
test('PDB: altloc, insertion codes, hetero, first model only',()=>{const s=parsePDB(PDB);assert.equal(s.atoms.length,5);
 assert.ok(!s.atoms.some(a=>a.x===9),'altloc B skipped');assert.equal(s.warnings.length,1);
 const keys=new Set(s.atoms.filter(a=>!a.hetero).map(residueKey));assert.deepEqual([...keys].sort(),['A:52','A:52A']);assert.equal(s.atoms.filter(a=>a.hetero).length,2);assert.deepEqual(s.chains,['A'])});
test('PDB: malformed / empty input',()=>{assert.throws(()=>parsePDB(PDB.replace('  50.000','  xx.xxx')),ParseError);assert.throws(()=>parsePDB('REMARK only'),/No ATOM/);
 assert.throws(()=>parseStructure('a.pdb','A\0B'),/Binary/);assert.throws(()=>parseStructure('a.mol','x'),/Unsupported/);assert.throws(()=>parsePDB(P('ATOM',1,' CA ',' ','ALA','A',1,' ',0,0,0,'C').replace('   0.000','        ')),/invalid/)});
const CIF=`data_t\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.pdbx_PDB_model_num\nATOM 1 C CA . SER A B 10 ? 1.0 2.0 3.0 1\nATOM 2 C CA . GLY A B 10 A 4.0 5.0 6.0 1\nHETATM 3 O "O1'" . LIG C D 7 ? 7.0 8.0 9.0 1\nATOM 4 C CA . SER A B 10 ? 0 0 0 2\n#\n`;
test('mmCIF: auth ids, insertion codes, quoted names, first model',()=>{const s=parseMmCIF(CIF);assert.equal(s.atoms.length,3);assert.equal(residueKey(s.atoms[1]),'B:10A');assert.equal(s.atoms[2].name,"O1'");assert.equal(s.atoms[2].chain,'D');assert.ok(s.atoms[2].hetero);
 assert.throws(()=>parseMmCIF('data_x\n'),/_atom_site/);assert.throws(()=>parseMmCIF(CIF.replace('4.0 5.0 6.0 1','4.0 5.0')),ParseError)});
const o={x:0,y:0,z:0},pt=(x,y,z)=>({x,y,z});
test('geometry reference cases',()=>{assert.equal(distance(o,pt(3,4,0)),5);assert.ok(Math.abs(angle(pt(1,0,0),o,pt(0,1,0))-90)<1e-9);assert.ok(Math.abs(angle(pt(1,0,0),o,pt(-1,0,0))-180)<1e-9);
 const a=pt(0,1,0),b=o,c=pt(1,0,0);assert.ok(Math.abs(dihedral(a,b,c,pt(1,1,0)))<1e-9,'cis=0');assert.ok(Math.abs(Math.abs(dihedral(a,b,c,pt(1,-1,0)))-180)<1e-9,'trans=180');
 assert.ok(Math.abs(dihedral(a,b,c,pt(1,0,1))-90)<1e-9,'IUPAC +90');assert.ok(Math.abs(dihedral(a,b,c,pt(1,0,-1))+90)<1e-9,'IUPAC -90');
 assert.throws(()=>angle(o,o,pt(1,0,0)),RangeError);assert.throws(()=>dihedral(o,pt(1,0,0),pt(2,0,0),pt(3,1,0)),RangeError)});
test('ligand contacts: geometric vs annotated evidence kept separate',()=>{const s=parsePDB(PDB),ann=new Map([['A:52','TEST-ANNOTATION'],['A:999','TEST-ANNOTATION']]);
 const r=ligandContacts(s.atoms,{ligandKey:'A:301',cutoff:4.0,annotated:ann});assert.equal(r.contacts.length,1);assert.equal(r.contacts[0].residue,'A:52');assert.ok(Math.abs(r.contacts[0].minDistance-Math.hypot(3,0.5))<1e-9,'min over atoms: CA(3,4,0)-C1(0,3.5,0)');
 
 assert.deepEqual(r.contacts[0].evidence,['geometric_contact','annotated_catalytic']);assert.equal(r.annotatedNotInContact[0].residue,'A:999');assert.equal(r.annotatedNotInContact[0].presentInStructure,false);
 const r2=ligandContacts(s.atoms,{ligandKey:'A:301',cutoff:4.0});assert.deepEqual(r2.contacts[0].evidence,['geometric_contact']);assert.match(r2.caveat,/not established catalytic/);assert.throws(()=>ligandContacts(s.atoms,{ligandKey:'Z:1'}),/not found/);assert.throws(()=>ligandContacts(s.atoms,{ligandKey:'A:301',cutoff:0}),RangeError)});
