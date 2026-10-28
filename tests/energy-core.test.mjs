import test from 'node:test';import assert from 'node:assert/strict';
import {ethaneEnergy,ethanePhase} from '../src/core/energy-core.mjs';
import {rotateBond} from '../src/core/conformation-core.mjs';
const atoms=[{x:-.77,y:0,z:0,el:'C'},{x:.77,y:0,z:0,el:'C'}];
for(const side of [0,1])for(let k=0;k<3;k++){const t=(120*k+side*60)*Math.PI/180;atoms.push({x:side?1.13:-1.13,y:Math.cos(t),z:Math.sin(t),el:'H'});}
const bonds=[[0,1],[0,2],[0,3],[0,4],[1,5],[1,6],[1,7]],orders=bonds.map(()=>1);
test('staggered minimum, eclipsed maximum, periodicity and range',()=>{assert.ok(ethaneEnergy(0,60)<1e-10);assert.equal(ethaneEnergy(60,60),12);for(let a=-180;a<=180;a++){const e=ethaneEnergy(a,60);assert.ok(e>=0&&e<=12);assert.ok(Math.abs(e-ethaneEnergy(a+120,60))<1e-10);}});
test('phase follows coordinates and rotation rather than assuming starting conformation',()=>{const phase=ethanePhase(atoms,bonds,orders,0,1);assert.notEqual(phase,null);assert.ok(ethaneEnergy(0,phase)<1e-8);const rotated=rotateBond(atoms,0,1,[1,5,6,7],37),p=ethanePhase(rotated,bonds,orders,0,1);assert.ok(Math.abs(ethaneEnergy(37,phase)-ethaneEnergy(0,p))<1e-9);});
test('reversing axis gives same current energy',()=>{const r=rotateBond(atoms,0,1,[1,5,6,7],23);assert.ok(Math.abs(ethaneEnergy(0,ethanePhase(r,bonds,orders,0,1))-ethaneEnergy(0,ethanePhase(r,bonds,orders,1,0)))<1e-9);});
test('reject missing hydrogens, substituted molecules and distorted projections',()=>{assert.equal(ethanePhase(atoms.slice(0,7),bonds,orders,0,1),null);const changed=atoms.map(a=>({...a}));changed[2].el='F';assert.equal(ethanePhase(changed,bonds,orders,0,1),null);changed[2]={...atoms[2],y:0,z:1};assert.equal(ethanePhase(changed,bonds,orders,0,1),null);});
